# GEMINI.md

Guidance for Gemini CLI when working in this repository.

## Project

**ParMetaD** — GPU-parallelizable metadynamics simulations driven by the [WESTPA](https://westpa.github.io/westpa/overview.html) weighted-ensemble framework. Pluggable MD engine: [GPUMD](https://gpumd.org/), [AMBER](https://ambermd.org/), or [OpenMM](https://openmm.org/). Metadynamics bias is applied through [PLUMED](https://www.plumed.org/).

This is a simulation **harness**, not a library. There is no build step and no test suite. "Running" means submitting `run_WE.sh` to a SLURM scheduler on a GPU cluster.

## Repository Map

- `run_WE.sh` — SLURM submission script. Sources `env.sh`, calls `init.sh`, starts a WESTPA ZMQ master, then SSHs into each allocated node to start workers via `node.sh`. **This is the entry point.**
- `init.sh` — Cleans previous run state (`west.h5`, `traj_segs/`, `seg_logs/`, `istates/`) and runs `w_init` with the basis-state list. Sources `env.sh`. Already invoked from `run_WE.sh` — do not call separately.
- `env.sh` — Cluster modules, conda env activation (`westpa-2.0`), and exports for `PLUMED_KERNEL`, `LD_LIBRARY_PATH`, `GPUMD_PATH`, etc. Sourced by both `init.sh` and `run_WE.sh`.
- `node.sh` — Per-node worker launcher used over SSH by `run_WE.sh`.
- `west.cfg` — WESTPA YAML config: pcoord dimensionality and dtype, bin boundaries, target walkers per bin, `max_total_iterations`, and the executable propagator wiring to `westpa_scripts/`.
- `bstates/` — Basis states. `bstates.txt` is `<index> <weight> <filename>`, one per line. Coordinate file (e.g. `model.xyz`) lives here too.
- `common_files/` — Files shared into every walker by `runseg.sh`: engine inputs for both "new trajectory" and "continued" segments, force-field / potential files, a reference topology for PLUMED `MOLINFO`, and PLUMED inputs (`plumed_init.dat` without `RESTART`, `plumed.dat` with `RESTART`).
- `westpa_scripts/` — WESTPA hook scripts:
  - `runseg.sh` — **The hot path.** Runs one walker segment. Branches on `WEST_CURRENT_SEG_INITPOINT_TYPE` (`SEG_INITPOINT_NEWTRAJ` vs `SEG_INITPOINT_CONTINUES`), launches the MD engine, then writes the progress coordinate to `$WEST_PCOORD_RETURN` (currently column 2 of PLUMED's `COLVAR`).
  - `get_pcoord.sh`, `gen_istate.sh`, `post_iter.sh`, `tar_segs.sh`, `cat_trajectory.py` — supporting hooks.
- `AMBER_files/`, `AMBER_files_new/`, `OpenMM_files/` — Engine-specific templates. Not all are wired up; the currently active path is GPUMD.
- `ParMetaD_FES.ipynb` — Free energy surface analysis (reads `west.h5` and HILLS data).
- `data_extract.py`, `simtime.py`, `run_data.sh` — Post-processing helpers.

## How Pieces Connect

`run_WE.sh` (SLURM) → sources `env.sh` → runs `init.sh` (which also sources `env.sh` and runs `w_init`) → starts `w_run` ZMQ master → SSHs `node.sh` on each node → workers pull tasks → WESTPA invokes `westpa_scripts/runseg.sh` per walker per iteration → `runseg.sh` writes pcoord back → WESTPA bins, resamples, advances iteration.

`west.cfg` is the contract between WESTPA and `westpa_scripts/`: the `executable` block names which script handles `propagator`, `get_pcoord`, `gen_istate`, `post_iteration`.

## Working Conventions

- **Cluster paths in `env.sh` and `run_WE.sh` are hardcoded** to one user's Expanse layout. Treat them as templates; when editing, replace only the paths the user names, don't sweep through "fixing" everything.
- **Two parallel inputs per engine**: one for `SEG_INITPOINT_NEWTRAJ` (cold start, seeds velocities, no PLUMED `RESTART`) and one for `SEG_INITPOINT_CONTINUES` (restart from parent, PLUMED `RESTART`, carries `HILLS` forward). When changing one, check whether the other needs the same change.
- **Progress coordinate plumbing lives in `runseg.sh`**, not in `west.cfg`. If the user changes the CV, update both `plumed*.dat` (the `PRINT ARG=...` line) and the awk column in `runseg.sh` that writes `$WEST_PCOORD_RETURN`, and update bin boundaries in `west.cfg`.
- **`init.sh` is destructive** — it deletes `west.h5`, `traj_segs/`, and `seg_logs/`. Never suggest running it on a simulation the user wants to keep. To extend an existing run use `w_truncate` / resubmit logic, not a re-init.
- **No tests, no linter, no CI.** "Validation" is a short SLURM job that completes a couple of iterations without errors in `seg_logs/`. When asked to verify a change, propose checking `job.err`, `west-*.log`, and a sample `seg_logs/<iter>-<seg>.log`.
- **Don't invent engine support.** If a directory exists (`AMBER_files/`, `OpenMM_files/`) but `runseg.sh` does not reference it, that engine is not currently wired up. Say so rather than assuming it works.

## Common Tasks

- **New CV** → edit `common_files/plumed.dat` and `plumed_init.dat` (CV definition, `METAD ARG=`, `PRINT ARG=`), update the awk column in `westpa_scripts/runseg.sh`, update bin boundaries and `pcoord_ndim` in `west.cfg`.
- **New MD engine** → add engine inputs under `common_files/`, add a branch in `runseg.sh` that copies/links them and invokes the engine, update `env.sh` with the engine's env vars.
- **Different cluster** → edit `module load` lines and exported paths in `env.sh`, then SLURM directives and SSH/CUDA assumptions in `run_WE.sh` and `node.sh`.
- **Analysis** → `ParMetaD_FES.ipynb` and `data_extract.py` read `west.h5` plus per-segment HILLS. `run_data.sh` is the convenience wrapper.

## What Not To Do

- Do not commit large generated outputs (`west.h5`, `traj_segs/`, `seg_logs/`, `istates/`).
- Do not edit `AMBER_files.zip` — it is a snapshot.
- Do not run `init.sh` "to be safe" before a resubmit; `run_WE.sh` already calls it.
- Do not rewrite paths in `env.sh` to anything portable on the user's behalf unless explicitly asked — the hardcoded paths are intentional per-deploy edits.
