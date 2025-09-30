#!/bin/bash
#SBATCH --job-name="chignolin_GaMD"
#SBATCH --output="job.out"
#SBATCH --partition=gpu-debug
#SBATCH --nodes=1
#SBATCH --gpus=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=50G
#SBATCH --account=ucd187
#SBATCH --no-requeue
#SBATCH --mail-user=ssiddharth@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH -t 00:30:00


module purge
module load shared
module load gpu/0.15.4
module load slurm
module load openmpi/4.0.4
module load cuda/11.0.2
module load plumed/2.6.1
#module load amber/20
conda init
conda activate westpa-2.0

export PATH=$PATH:$HOME/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH
#source $AMBERHOME/amber.sh
export PYTHONPATH=/home/ssonti/miniconda3/envs/westpa-2.0/bin/python
export PLUMED_KERNEL=/cm/shared/apps/spack/gpu/opt/spack/linux-centos8-skylake_avx512/gcc-8.3.1/plumed-2.6.1-63lfaa2clqpjeif3aa3kdk44ozvlzkac/lib/libplumedKernel.so
echo $PLUMED_KERNEL
#pmemd.cuda -O -i md.in -o md.out -p chignolin.prmtop -c chignolin.rst -r md_cmd.rst -x md.nc
#python3 check_forcefields.py
#python chignolin_metadynamics.py -p chignolin_works.pdb -l plumed_rmsd.dat --final-pdb-only -o final.pdb
#python -u chignolin_metadynamics_implicit.py --prmtop chignolin.prmtop --inpcrd chignolin_box.inpcrd --plumed plumed_rmsd.dat > log.out
#python chignolin_metadynamics_implicit.py \
 #   --prmtop chignolin.prmtop --inpcrd chignolin_box.inpcrd \
 #   --plumed plumed_rmsd.dat -n 500000 -f 500 \
  #  -d traj_implicit.dcd --final-pdb-only -o final.pdb 

python chignolin_metadynamics_yaml.py -c metadynamics_config.yaml
