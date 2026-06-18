# ParMetaD

GPU-parallelizable metadynamics simulations powered by the [WESTPA](https://westpa.github.io/westpa/overview.html) weighted-ensemble framework. Supports [GPUMD](https://gpumd.org/), [AMBER](https://ambermd.org/), and [OpenMM](https://openmm.org/) as the underlying MD engine, with [PLUMED](https://www.plumed.org/) providing the metadynamics bias.

This repo is not final, and is actively undergoing changes.

## Requirements

- WESTPA 2.0 (a conda env named `westpa-2.0` is assumed by the scripts)
- PLUMED 2.x compiled as a shared library (`libplumedKernel.so`)
- One MD engine: GPUMD, AMBER, or OpenMM
- A SLURM-managed cluster with GPU nodes (the provided `run_WE.sh` is written for SLURM)
- Python 3.8+, NumPy

## Repository Layout

- `bstates/` — basis states: starting structures WESTPA samples walkers from
- `common_files/` — engine inputs, PLUMED inputs, and any other files shared across walkers
- `westpa_scripts/` — WESTPA hook scripts (`runseg.sh`, `get_pcoord.sh`, `gen_istate.sh`, `post_iter.sh`)
- `AMBER_files/`, `AMBER_files_new/`, `OpenMM_files/` — engine-specific templates and examples
- `west.cfg` — WESTPA configuration (progress coordinate, bins, propagator)
- `env.sh`, `init.sh`, `run_WE.sh`, `node.sh` — environment, initialization, and run scripts
- `ParMetaD_FES.ipynb` — free energy surface analysis notebook
- `data_extract.py`, `simtime.py`, `run_data.sh` — post-processing helpers

## Setup

### 1. Edit `env.sh`

Update the cluster-specific paths and the conda env name:

```bash
conda activate westpa-2.0                       # your westpa conda env
export LD_LIBRARY_PATH=/path/to/plumed/lib:$LD_LIBRARY_PATH
export PLUMED_KERNEL=/path/to/plumed/lib/libplumedKernel.so
# Plus any engine-specific variables, e.g.:
# export GPUMD_PATH=/path/to/GPUMD/src/gpumd
# export AMBERHOME=/path/to/amber
```

Also adjust the `module load` lines to match your cluster's module system.

### 2. Place your initial structure(s) in `bstates/`

`bstates/` holds the basis states WESTPA initializes walkers from. Drop the starting coordinate/restart file for your engine here — for example a `.xyz` for GPUMD, a `.rst`/`.ncrst` for AMBER, or a `.pdb`/state file for OpenMM.

Then list each basis state in `bstates/bstates.txt`, one per line, in the format:

```
<index> <weight> <filename>
```

Example:

```
0 1 model.xyz
```

To use multiple basis states, add more rows and update `--segs-per-state` in `init.sh` accordingly.

### 3. Populate `common_files/`

`common_files/` is symlinked or copied into every walker's working directory by `westpa_scripts/runseg.sh`. Put files shared across all walkers here:

- **Engine input files** for both a fresh segment and a continued segment (separate inputs are needed because the first segment has to initialize velocities / state, while later segments restart from the parent walker). Examples: `run_init.in` and `run.in` for GPUMD, `md.in` and a restart-aware variant for AMBER, an OpenMM driver script, etc.
- **Force-field or potential files**: NEP potential for GPUMD, `prmtop` / parameter files for AMBER, system XML or forcefield XML for OpenMM.
- **Topology / reference structure** that PLUMED can read for atom names (e.g. a reference PDB for PLUMED's `MOLINFO`).
- **PLUMED inputs**: typically a `plumed_init.dat` (used on the first segment, no `RESTART`) and a `plumed.dat` (with `RESTART`, used on continued segments). Define your collective variables, `METAD` parameters (`SIGMA`, `HEIGHT`, `BIASFACTOR`, `TEMP`, `PACE`, `GRID_*`), any walls, and a `PRINT` line that writes the CV to `COLVAR`.

### 4. Edit `west.cfg`

Set the progress coordinate dimensionality, bin boundaries, target walkers per bin, and `max_total_iterations` to match your CV and the run length you want.

### 5. Edit `westpa_scripts/runseg.sh`

`runseg.sh` is what each walker actually executes. Confirm:

- The right engine input, potential/parameter files, topology, and PLUMED files are linked or copied for your engine.
- The branch for `SEG_INITPOINT_NEWTRAJ` (fresh start) uses your "init" inputs and links the basis-state coordinates, while `SEG_INITPOINT_CONTINUES` uses the restart-aware inputs and pulls the restart file and HILLS from `$WEST_PARENT_DATA_REF`.
- The MD engine is invoked correctly for your setup.
- The GPU assignment block (`CUDA_VISIBLE_DEVICES`) matches your hardware.
- The progress coordinate is pulled from the correct column of the right output file and written to `$WEST_PCOORD_RETURN`. The shipped version reads column 2 of PLUMED's `COLVAR`:

  ```bash
  tail -n +2 COLVAR | awk '{print $2}' > $WEST_PCOORD_RETURN
  ```

### 6. Edit `run_WE.sh`

Update SLURM directives (account, partition, nodes, GPUs, wall time, mail) and the `num_gpu_per_node` variable to match your allocation.

## Running

`run_WE.sh` already sources `env.sh` and invokes `init.sh` (which cleans previous state, recreates `traj_segs/`, `seg_logs/`, `istates/`, and calls `w_init`), so you only need to submit it:

```bash
sbatch run_WE.sh
```

To re-run from scratch, just resubmit — `init.sh` removes `west.h5`, `traj_segs/`, and `seg_logs/` on each run. Do **not** call `init.sh` manually before submitting unless you want a fresh initialization; running it twice has no extra effect but is unnecessary.

## Analysis

After the run finishes, extract data and plot the free energy surface:

```bash
./run_data.sh           # runs data_extract.py / simtime.py
jupyter notebook ParMetaD_FES.ipynb
```

## License

MIT — see [LICENSE](LICENSE).
