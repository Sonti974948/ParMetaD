#!/bin/bash
#
# node.sh


set -x
umask g+r
cd $1; shift
source env.sh
export WEST_JOBID=$1; shift
export SLURM_NODENAME=$1; shift
export CUDA_VISIBLE_DEVICES_ALLOCATED=$1; shift
echo "starting WEST client processes on: "; hostname
echo "current directory is $PWD"
echo "environment is: "
env | sort


# Extract individual GPUs from CUDA_VISIBLE_DEVICES_ALLOCATED
#export CUDA_DEVICES=(`echo $CUDA_VISIBLE_DEVICES_ALLOCATED | tr , ' '`)

# Set CUDA_VISIBLE_DEVICES based on worker index (WM_PROCESS_INDEX should be set for each worker)
#export CUDA_VISIBLE_DEVICES=${CUDA_DEVICES[$WM_PROCESS_INDEX]}

#echo "Assigned GPU (CUDA_VISIBLE_DEVICES) = " $CUDA_VISIBLE_DEVICES
#echo "Assigned GPU (CUDA_VISIBLE_DEVICES) for walker $WM_PROCESS_INDEX = " $CUDA_VISIBLE_DEVICES
#w_run "$@" &> west-$SLURM_NODENAME-node.log
#echo "Shutting down. Hopefully this was on purpose?"

echo "CUDA_VISIBLE_DEVICES = " $CUDA_VISIBLE_DEVICES
w_run "$@" &> west-$SLURM_NODENAME-node.log
echo "Shutting down.  Hopefully this was on purpose?"
