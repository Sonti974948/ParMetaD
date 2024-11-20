#!/bin/bash

if [ -n "$SEG_DEBUG" ] ; then
  set -x
  env | sort
fi

cd $WEST_SIM_ROOT
mkdir -pv $WEST_CURRENT_SEG_DATA_REF
cd $WEST_CURRENT_SEG_DATA_REF

ln -sv $WEST_SIM_ROOT/common_files/initial_traj.pdb .
#ln -sv $WEST_SIM_ROOT/common_files/plumed.dat .
ln -sv $WEST_SIM_ROOT/common_files/nep5.txt .

if [ "$WEST_CURRENT_SEG_INITPOINT_TYPE" = "SEG_INITPOINT_CONTINUES" ]; then
  cp $WEST_SIM_ROOT/common_files/run.in run.in
  cp $WEST_SIM_ROOT/common_files/plumed.dat plumed.dat
  ln -sv $WEST_PARENT_DATA_REF/restart.xyz ./model.xyz
  cp $WEST_PARENT_DATA_REF/HILLS HILLS
elif [ "$WEST_CURRENT_SEG_INITPOINT_TYPE" = "SEG_INITPOINT_NEWTRAJ" ]; then
  cp $WEST_SIM_ROOT/common_files/run_init.in run.in
  cp $WEST_SIM_ROOT/common_files/plumed_init.dat plumed.dat
  ln -sv $WEST_PARENT_DATA_REF ./model.xyz
fi

export CUDA_DEVICES=(`echo $CUDA_VISIBLE_DEVICES_ALLOCATED | tr , ' '`)
export CUDA_VISIBLE_DEVICES=${CUDA_DEVICES[$WM_PROCESS_INDEX]}

echo "RUNSEG.SH: CUDA_VISIBLE_DEVICES_ALLOCATED = " $CUDA_VISIBLE_DEVICES_ALLOCATED
echo "RUNSEG.SH: WM_PROCESS_INDEX = " $WM_PROCESS_INDEX
echo "RUNSEG.SH: CUDA_VISIBLE_DEVICES = " $CUDA_VISIBLE_DEVICES

while ! grep -q "Finished executing the commands in run.in" log.out; do
	$GPUMD_PATH > log.out
	min_value=$(awk 'NR==1 {min=$2} NR>1 && $2 < min {min=$2} END {print min}' HILLS)
	max_value=$(awk '$2+0==$2 {if (NR==1 || $2 > max) max=$2} END {print max}' HILLS)
	plumed sum_hills --hills HILLS --mintozero --min "$min_value" --max "$max_value" --bin 100
done

#export CUDA_DEVICES=(`echo $CUDA_VISIBLE_DEVICES_ALLOCATED | tr , ' '`)
#export CUDA_VISIBLE_DEVICES=${CUDA_DEVICES[$WM_PROCESS_INDEX]}

#echo "RUNSEG.SH: CUDA_VISIBLE_DEVICES_ALLOCATED = " $CUDA_VISIBLE_DEVICES_ALLOCATED
#echo "RUNSEG.SH: WM_PROCESS_INDEX = " $WM_PROCESS_INDEX
#echo "RUNSEG.SH: CUDA_VISIBLE_DEVICES = " $CUDA_VISIBLE_DEVICES

#cat $RMSD > rmsd.dat
#cat $RG > rg.dat
#paste <(cat COLVAR | tail -n +2 | awk {'print $1'}) >$WEST_PCOORD_RETURN
tail -n +2 COLVAR | awk '{print $2}' > $WEST_PCOORD_RETURN
#cat $TEMP | tail -n +2 | awk '{print $2}' > $WEST_PCOORD_RETURN
#paste <(cat $TEMP | tail -n 1 | awk {'print $2'}) <(cat $RG | tail -n 1 | awk {'print $2'})>$WEST_PCOORD_RETURN
#cat $TEMP >pcoord.dat
# Clean up
#rm -f $TEMP md.in seg.nfo seg.pdb
