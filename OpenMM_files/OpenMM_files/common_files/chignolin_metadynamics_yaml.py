#!/usr/bin/env python3
"""
OpenMM Metadynamics with YAML Configuration
Reads simulation parameters from YAML file (similar to AMBER .in files)
"""

import os
import sys
import argparse
import yaml
from pathlib import Path
import datetime

from openmm.app import (
    AmberPrmtopFile,
    AmberInpcrdFile,
    PDBFile,
    Simulation,
    HBonds,
    HCT,
    OBC1,
    OBC2,
    GBn,
    GBn2,
    NoCutoff,
    CutoffNonPeriodic,
    PME,
    DCDReporter,
    PDBReporter,
    StateDataReporter,
)
from openmm import Platform, LangevinIntegrator
from openmm.unit import kelvin, picoseconds, amu
from openmmplumed import PlumedForce
import mdtraj.reporters


def load_config(config_file):
    """Load configuration from YAML file"""
    if not os.path.exists(config_file):
        raise FileNotFoundError(f"Configuration file not found: {config_file}")
    
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    
    return config


def setup_run_log(config, config_file):
    """Setup run.log file with general simulation information"""
    run_log_file = "run.log"
    
    with open(run_log_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("OPENMM METADYNAMICS SIMULATION RUN LOG\n")
        f.write("=" * 80 + "\n")
        f.write(f"Start time: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Configuration file: {config_file}\n")
        f.write(f"Working directory: {os.getcwd()}\n")
        f.write("\n")
        
        # Simulation mode and files
        f.write("SIMULATION SETUP:\n")
        f.write("-" * 40 + "\n")
        f.write(f"Mode: {config['mode']}\n")
        
        input_files = config['input_files']
        f.write(f"Topology file: {input_files['prmtop']}\n")
        f.write(f"PLUMED file: {input_files['plumed']}\n")
        
        if config['mode'] == 'init':
            if 'inpcrd' in input_files and input_files['inpcrd']:
                f.write(f"Coordinate file: {input_files['inpcrd']} (inpcrd)\n")
            elif 'pdb' in input_files and input_files['pdb']:
                f.write(f"Coordinate file: {input_files['pdb']} (pdb)\n")
        elif config['mode'] == 'restart':
            if 'checkpoint' in input_files and input_files['checkpoint']:
                f.write(f"Restart file: {input_files['checkpoint']} (checkpoint)\n")
            elif 'nc_restart' in input_files and input_files['nc_restart']:
                f.write(f"Restart file: {input_files['nc_restart']} (netcdf)\n")
        
        f.write("\n")
        
        # MD parameters
        md_params = config['md_parameters']
        f.write("MOLECULAR DYNAMICS PARAMETERS:\n")
        f.write("-" * 40 + "\n")
        f.write(f"Platform: {md_params['platform']}\n")
        f.write(f"Temperature: {md_params['temperature']} K\n")
        f.write(f"Timestep: {md_params['timestep']} ps\n")
        f.write(f"Number of steps: {md_params['nsteps']:,}\n")
        f.write(f"Output frequency: {md_params['output_freq']} steps\n")
        f.write(f"Implicit solvent: {md_params['implicit_solvent']}\n")
        f.write(f"Nonbonded method: {md_params['nonbonded_method']}\n")
        f.write(f"Constraints: {md_params['constraints']}\n")
        f.write(f"Hydrogen mass: {md_params['hydrogen_mass']} amu\n")
        f.write("\n")
        
        # Output files
        output_files = config['output_files']
        f.write("OUTPUT FILES:\n")
        f.write("-" * 40 + "\n")
        f.write(f"DCD trajectory: {output_files['dcd_trajectory']}\n")
        f.write(f"PDB snapshots: {output_files['pdb_snapshots']}\n")
        f.write(f"Simulation log: {output_files['log_file']}\n")
        if 'output_checkpoint' in output_files and output_files['output_checkpoint']:
            f.write(f"Output checkpoint: {output_files['output_checkpoint']}\n")
        if 'nc_restart' in output_files and output_files['nc_restart']:
            f.write(f"NetCDF restart: {output_files['nc_restart']}\n")
        f.write("\n")
        
        # Advanced options
        advanced = config.get('advanced', {})
        if advanced:
            f.write("ADVANCED OPTIONS:\n")
            f.write("-" * 40 + "\n")
            f.write(f"Energy minimization: {advanced.get('minimize_energy', True)}\n")
            f.write(f"Minimization iterations: {advanced.get('min_iterations', 1000)}\n")
            f.write(f"GPU device: {advanced.get('gpu_device', 0)}\n")
            f.write(f"Precision: {advanced.get('precision', 'mixed')}\n")
            f.write("\n")
        
        f.write("=" * 80 + "\n")
        f.write("SIMULATION PROGRESS:\n")
        f.write("=" * 80 + "\n")
    
    return run_log_file


def log_to_run_log(message, run_log_file="run.log"):
    """Append message to run.log file"""
    with open(run_log_file, 'a') as f:
        timestamp = datetime.datetime.now().strftime('%H:%M:%S')
        f.write(f"[{timestamp}] {message}\n")
        f.flush()  # Ensure immediate write


def create_implicit_system(config):
    """Create OpenMM system from configuration"""
    input_files = config['input_files']
    md_params = config['md_parameters']
    
    prmtop_path = input_files['prmtop']
    if not os.path.exists(prmtop_path):
        raise FileNotFoundError(f"prmtop not found: {prmtop_path}")
    
    prmtop = AmberPrmtopFile(prmtop_path)
    
    # Determine coordinate file and type
    mode = config['mode']
    positions = None
    checkpoint_file = None
    
    if mode == "init":
        if 'inpcrd' in input_files and input_files['inpcrd']:
            coord_file = input_files['inpcrd']
            if not os.path.exists(coord_file):
                raise FileNotFoundError(f"inpcrd not found: {coord_file}")
            inpcrd = AmberInpcrdFile(coord_file)
            positions = inpcrd.positions
        elif 'pdb' in input_files and input_files['pdb']:
            coord_file = input_files['pdb']
            if not os.path.exists(coord_file):
                raise FileNotFoundError(f"pdb not found: {coord_file}")
            pdb = PDBFile(coord_file)
            positions = pdb.positions
        else:
            raise ValueError("init mode requires either inpcrd or pdb in input_files")
    
    elif mode == "restart":
        if 'checkpoint' in input_files and input_files['checkpoint']:
            checkpoint_file = input_files['checkpoint']
            if not os.path.exists(checkpoint_file):
                raise FileNotFoundError(f"checkpoint not found: {checkpoint_file}")
        elif 'nc_restart' in input_files and input_files['nc_restart']:
            coord_file = input_files['nc_restart']
            if not os.path.exists(coord_file):
                raise FileNotFoundError(f"nc_restart not found: {coord_file}")
            inpcrd = AmberInpcrdFile(coord_file)
            positions = inpcrd.positions
        else:
            raise ValueError("restart mode requires either checkpoint or nc_restart in input_files")
    
    # Create system with specified parameters
    implicit_solvent_map = {
        'HCT': HCT,
        'OBC1': OBC1,
        'OBC2': OBC2,
        'GBn': GBn,
        'GBn2': GBn2
    }
    
    nonbonded_method_map = {
        'NoCutoff': NoCutoff,
        'CutoffNonPeriodic': CutoffNonPeriodic,
        'PME': PME
    }
    
    constraints_map = {
        'None': None,
        'HBonds': HBonds,
        'AllBonds': 'AllBonds',  # OpenMM constant
        'HAngles': 'HAngles'     # OpenMM constant
    }
    
    system = prmtop.createSystem(
        nonbondedMethod=nonbonded_method_map[md_params['nonbonded_method']],
        constraints=constraints_map[md_params['constraints']],
        implicitSolvent=implicit_solvent_map[md_params['implicit_solvent']],
        hydrogenMass=md_params['hydrogen_mass'] * amu,
    )
    
    topology = prmtop.topology
    return system, topology, positions, checkpoint_file


def add_plumed(system, plumed_path: str) -> PlumedForce:
    """Add PLUMED force to system"""
    if not os.path.exists(plumed_path):
        raise FileNotFoundError(f"PLUMED script not found: {plumed_path}")
    
    with open(plumed_path, "r") as f:
        script = f.read()
    
    pf = PlumedForce(script)
    system.addForce(pf)
    return pf


def setup_simulation(system, topology, positions, config, checkpoint_file=None):
    """Setup simulation from configuration"""
    md_params = config['md_parameters']
    
    # Set platform
    platform = Platform.getPlatformByName(md_params['platform'])
    
    # Create integrator
    temperature = md_params['temperature'] * kelvin
    timestep = md_params['timestep'] * picoseconds
    integrator = LangevinIntegrator(temperature, 1.0 / picoseconds, timestep)
    
    # Create simulation
    sim = Simulation(topology, system, integrator, platform)
    
    if checkpoint_file:
        # Restart from checkpoint (loads positions, velocities, integrator state)
        sim.loadCheckpoint(checkpoint_file)
        print(f"Loaded checkpoint from {checkpoint_file}")
        print("Restarting with velocities from checkpoint")
    else:
        # Fresh start (set positions and generate new velocities)
        sim.context.setPositions(positions)
        sim.context.setVelocitiesToTemperature(temperature)
        print(f"Fresh start: generated velocities at {temperature}")
    
    return sim


def setup_reporters(simulation, config):
    """Setup output reporters from configuration"""
    output_files = config['output_files']
    md_params = config['md_parameters']
    advanced = config.get('advanced', {})
    
    # State data reporter
    simulation.reporters.append(
        StateDataReporter(
            output_files['log_file'],
            md_params['output_freq'],
            step=True,
            time=True,
            potentialEnergy=True,
            kineticEnergy=True,
            totalEnergy=True,
            temperature=True,
        )
    )
    
    # DCD trajectory reporter
    simulation.reporters.append(
        DCDReporter(output_files['dcd_trajectory'], md_params['output_freq'])
    )
    
    # PDB snapshots reporter
    if output_files['final_pdb_only']:
        simulation.reporters.append(
            PDBReporter(output_files['pdb_snapshots'], md_params['nsteps'])
        )
    else:
        simulation.reporters.append(
            PDBReporter(output_files['pdb_snapshots'], md_params['output_freq'] * 10)
        )
    
    # NetCDF restart reporter
    if 'nc_restart' in output_files and output_files['nc_restart']:
        nc_reporter = mdtraj.reporters.NetCDFReporter(
            output_files['nc_restart'], 
            md_params['nsteps'], 
            coordinates=True, 
            time=True, 
            cell=True
        )
        simulation.reporters.append(nc_reporter)


def run_simulation(simulation, config):
    """Run the simulation with configuration"""
    md_params = config['md_parameters']
    advanced = config.get('advanced', {})
    output_files = config['output_files']
    
    # Energy minimization
    if advanced.get('minimize_energy', True):
        print("Minimizing energy...")
        log_to_run_log("Starting energy minimization", "run.log")
        simulation.minimizeEnergy(maxIterations=advanced.get('min_iterations', 1000))
        log_to_run_log("Energy minimization completed", "run.log")
    
    # Run simulation
    print(f"Running simulation for {md_params['nsteps']} steps...")
    log_to_run_log(f"Starting MD run for {md_params['nsteps']:,} steps", "run.log")
    simulation.step(md_params['nsteps'])
    log_to_run_log("MD run completed", "run.log")
    
    # Save checkpoint
    if 'output_checkpoint' in output_files and output_files['output_checkpoint']:
        simulation.saveCheckpoint(output_files['output_checkpoint'])
        print(f"Checkpoint saved to {output_files['output_checkpoint']}")
        log_to_run_log(f"Checkpoint saved to {output_files['output_checkpoint']}", "run.log")




def main():
    parser = argparse.ArgumentParser(
        description="OpenMM Metadynamics with YAML Configuration",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run with YAML config
  python chignolin_metadynamics_yaml.py -c metadynamics_config.yaml
  
  # Override config file
  python chignolin_metadynamics_yaml.py -c config.yaml --nsteps 500000
        """
    )
    
    parser.add_argument('-c', '--config', required=True, help='YAML configuration file')
    parser.add_argument('--nsteps', type=int, help='Override number of steps')
    parser.add_argument('--temperature', type=float, help='Override temperature')
    parser.add_argument('--output-checkpoint', help='Override output checkpoint file')
    
    args = parser.parse_args()
    
    # Load configuration
    config = load_config(args.config)
    
    # Override parameters from command line
    if args.nsteps:
        config['md_parameters']['nsteps'] = args.nsteps
    if args.temperature:
        config['md_parameters']['temperature'] = args.temperature
    if args.output_checkpoint:
        config['output_files']['output_checkpoint'] = args.output_checkpoint
    
    # Setup run log
    run_log_file = setup_run_log(config, args.config)
    
    print("=" * 60)
    print("OpenMM Metadynamics with YAML Configuration")
    print("=" * 60)
    print(f"Mode: {config['mode']}")
    print(f"Platform: {config['md_parameters']['platform']}")
    print(f"Temperature: {config['md_parameters']['temperature']} K")
    print(f"Steps: {config['md_parameters']['nsteps']}")
    print(f"Implicit solvent: {config['md_parameters']['implicit_solvent']}")
    print(f"Run log: {run_log_file}")
    print("=" * 60)
    
    # Log initial setup
    log_to_run_log("Starting simulation setup", run_log_file)
    
    try:
        # Create system
        log_to_run_log("Creating OpenMM system", run_log_file)
        system, topology, positions, checkpoint_file = create_implicit_system(config)
        log_to_run_log(f"System created with {system.getNumParticles()} particles", run_log_file)
        
        # Add PLUMED
        log_to_run_log("Adding PLUMED force", run_log_file)
        add_plumed(system, config['input_files']['plumed'])
        log_to_run_log("PLUMED force added successfully", run_log_file)
        
        # Setup simulation
        log_to_run_log("Setting up simulation", run_log_file)
        simulation = setup_simulation(system, topology, positions, config, checkpoint_file)
        log_to_run_log("Simulation setup complete", run_log_file)
        
        # Setup reporters
        log_to_run_log("Setting up output reporters", run_log_file)
        setup_reporters(simulation, config)
        log_to_run_log("Reporters configured", run_log_file)
        
        # Run simulation
        log_to_run_log(f"Starting MD simulation for {config['md_parameters']['nsteps']:,} steps", run_log_file)
        run_simulation(simulation, config)
        log_to_run_log("MD simulation completed", run_log_file)
        
        print("\nSimulation completed successfully!")
        print("\nOutput files:")
        for key, filename in config['output_files'].items():
            if filename:
                print(f"  {key}: {filename}")
        
        # Log completion
        log_to_run_log("Simulation completed successfully", run_log_file)
        log_to_run_log(f"End time: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", run_log_file)
        
    except Exception as e:
        error_msg = f"Error during simulation: {e}"
        print(error_msg)
        log_to_run_log(f"ERROR: {error_msg}", run_log_file)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
