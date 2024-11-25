#!/bin/bash
#SBATCH --job-name="NEP_WE_run"
#SBATCH --output="job.out"
#SBATCH --error="job.err"
#SBATCH --partition=gpu
#SBATCH --nodes=4
#SBATCH --gpus=16
#SBATCH --ntasks-per-node=1
#SBATCH --mem=50G
#SBATCH --account=ucd192
#SBATCH --no-requeue
#SBATCH --mail-user=ssiddharth@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH -t 24:00:00



set -x
#export WEST_SIM_ROOT=$PWD
cd $SLURM_SUBMIT_DIR
#cd $WEST_SIM_ROOT
source ~/.bashrc
module purge
module load shared
module load gpu/0.15.4
module load slurm
module load openmpi/4.0.4
module load cuda/11.0.2
#module load amber/20-patch15
#module load conda
conda activate westpa-2.0
#source env.sh || exit 1
#env | sort

export LD_LIBRARY_PATH=/expanse/lustre/scratch/ssonti/temp_project/amber_learn/chignolin_tutorial/westpa_tutorials/tutorial7.3-chignolin/plumed/lib:$LD_LIBRARY_PATH
export WEST_SIM_ROOT=$SLURM_SUBMIT_DIR
cd $WEST_SIM_ROOT
export PYTHONPATH=/home/ssonti/miniconda3/envs/westpa-2.0/bin/python
export GPUMD_PATH=/expanse/lustre/scratch/ssonti/temp_project/amber_learn/chignolin_tutorial/westpa_tutorials/tutorial7.3-chignolin/GPUMD-3.9.4/src/gpumd
echo $GPUMD_PATH
export PLUMED_KERNEL=/expanse/lustre/scratch/ssonti/temp_project/amber_learn/chignolin_tutorial/westpa_tutorials/tutorial7.3-chignolin/plumed/lib/libplumedKernel.so
echo $PLUMED_KERNEL
#cp cMD/gamd-restart.dat common_files/gamd-restart.dat
#cp cMD/md_cmd.rst bstates/bstate.rst

./init.sh
echo "init.sh ran"
source env.sh || exit 1
env | sort
SERVER_INFO=$WEST_SIM_ROOT/west_zmq_info.json

#TODO: set num_gpu_per_node
num_gpu_per_node=4
rm -rf nodefilelist.txt
scontrol show hostname $SLURM_JOB_NODELIST > nodefilelist.txt

# start server
#w_truncate -n 11
#rm -rf traj_segs/000011
#rm -rf seg_logs/000011*
w_run --work-manager=zmq --n-workers=0 --zmq-mode=master --zmq-write-host-info=$SERVER_INFO --zmq-comm-mode=tcp &> west-$SLURM_JOBID-local.log &

# wait on host info file up to 1 min
for ((n=0; n<60; n++)); do
    if [ -e $SERVER_INFO ] ; then
        echo "== server info file $SERVER_INFO =="
        cat $SERVER_INFO
        break
    fi
    sleep 1
done

# exit if host info file doesn't appear in one minute
if ! [ -e $SERVER_INFO ] ; then
    echo 'server failed to start'
    exit 1
fi
export CUDA_VISIBLE_DEVICES=0,1,2,3
echo$CUDA_VISIBLE_DEVICES
for node in $(cat nodefilelist.txt); do
    ssh -o StrictHostKeyChecking=no $node $PWD/node.sh $SLURM_SUBMIT_DIR $SLURM_JOBID $node $CUDA_VISIBLE_DEVICES --work-manager=zmq --n-workers=$num_gpu_per_node --zmq-mode=client --zmq-read-host-info=$SERVER_INFO --zmq-comm-mode=tcp &
done
wait
