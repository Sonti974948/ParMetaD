#!/bin/bash

if [ -n "$SEG_DEBUG" ] ; then
  set -x
  env | sort
fi


#module purge

#module load shared
#module load gpu/0.15.4
#module load slurm
#module load openmpi/4.0.4
#module load cuda/11.0.2
#module load plumed/2.6.1
#module load amber/20
#conda init
#conda activate westpa-openmm

#export PATH=$PATH:$HOME/bin
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH
#source $AMBERHOME/amber.sh
#export PYTHONPATH=/home/ssonti/miniconda3/envs/westpa-2.0/bin/python
#export PLUMED_KERNEL=/cm/shared/apps/spack/gpu/opt/spack/linux-centos8-skylake_avx512/gcc-8.3.1/plumed-2.6.1-63lfaa2clqpjeif3aa3kdk44ozvlzkac/lib/libplumedKernel.so
#echo $PLUMED_KERNEL


cd $WEST_SIM_ROOT
mkdir -pv $WEST_CURRENT_SEG_DATA_REF
cd $WEST_CURRENT_SEG_DATA_REF

ln -sv $WEST_SIM_ROOT/common_files/chignolin.prmtop .
#ln -sv $WEST_SIM_ROOT/common_files/gamd-restart.dat .

if [ "$WEST_CURRENT_SEG_INITPOINT_TYPE" = "SEG_INITPOINT_CONTINUES" ]; then
  #sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/common_files/md.in > md.in
  cp $WEST_SIM_ROOT/common_files/plumed.dat plumed.dat
  cp $WEST_SIM_ROOT/common_files/chignolin_works.pdb .
  cp $WEST_PARENT_DATA_REF/final.chk init.chk
  cp $WEST_SIM_ROOT/common_files/chignolin_metadynamics_yaml.py .
  cp $WEST_SIM_ROOT/common_files/config_cont.yaml config.yaml
  cp $WEST_PARENT_DATA_REF/HILLS HILLS
elif [ "$WEST_CURRENT_SEG_INITPOINT_TYPE" = "SEG_INITPOINT_NEWTRAJ" ]; then
  #sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/common_files/md.in > md.in
  cp $WEST_SIM_ROOT/common_files/plumed_init.dat plumed.dat
  cp $WEST_SIM_ROOT/common_files/chignolin_works.pdb .
  cp $WEST_SIM_ROOT/common_files/chignolin_box.inpcrd .
  cp $WEST_SIM_ROOT/common_files/chignolin_metadynamics_yaml.py .
  cp $WEST_SIM_ROOT/common_files/config_init.yaml config.yaml
  cp $WEST_SIM_ROOT/COLVAR COLVAR0
fi

export CUDA_DEVICES=(`echo $CUDA_VISIBLE_DEVICES_ALLOCATED | tr , ' '`)
export CUDA_VISIBLE_DEVICES=${CUDA_DEVICES[$WM_PROCESS_INDEX]}

echo "RUNSEG.SH: CUDA_VISIBLE_DEVICES_ALLOCATED = " $CUDA_VISIBLE_DEVICES_ALLOCATED
echo "RUNSEG.SH: WM_PROCESS_INDEX = " $WM_PROCESS_INDEX
echo "RUNSEG.SH: CUDA_VISIBLE_DEVICES = " $CUDA_VISIBLE_DEVICES

while ! grep -q "Start time:" run.log; do
	python chignolin_metadynamics_yaml.py -c config.yaml
done

#RMSD=rmsd.dat
#RG=rg.dat
#COMMAND="         parm chignolin.prmtop\n"
#COMMAND="${COMMAND} trajin $WEST_CURRENT_SEG_DATA_REF/parent.rst\n"
#COMMAND="${COMMAND} trajin $WEST_CURRENT_SEG_DATA_REF/seg.nc\n"
#COMMAND="${COMMAND} reference $WEST_SIM_ROOT/common_files/chignolin_works.pdb\n"
#COMMAND="${COMMAND} rms ca-rmsd @CA reference out $RMSD mass\n"
#COMMAND="${COMMAND} radgyr ca-rg @CA  out $RG  mass\n"
#COMMAND="${COMMAND} go\n"

#echo -e $COMMAND | $CPPTRAJ
#cat $RMSD > rmsd.dat
#cat $RG > rg.dat
#paste <(cat rmsd.dat | tail -n +2 | awk {'print $2'}) <(cat rg.dat | tail -n +2 | awk {'print $2'})>$WEST_PCOORD_RETURN
#if [ "$WEST_CURRENT_SEG_INITPOINT_TYPE" = "SEG_INITPOINT_CONTINUES" ]; then
	#( awk 'NR==2 {print $2}' $WEST_SIM_ROOT/COLVAR; tail -n +2 COLVAR | awk '{print $2}' ) > $WEST_PCOORD_RETURN
#elif [ "$WEST_CURRENT_SEG_INITPOINT_TYPE" = "SEG_INITPOINT_NEWTRAJ" ]; then
	#tail -n +2 COLVAR | awk '{print $2}' > $WEST_PCOORD_RETURN
#fi
#tail -n +2 COLVAR | awk '{print $2}' > $WEST_PCOORD_RETURN
if [ -f COLVAR0 ]; then
    ( awk 'NR==2 {print $2}' COLVAR0; tail -n +2 COLVAR | awk '{print $2}' ) > $WEST_PCOORD_RETURN
else
    tail -n +2 COLVAR | awk '{print $2}' > $WEST_PCOORD_RETURN
fi

#cat $TEMP | tail -n +2 | awk '{print $2}' > $WEST_PCOORD_RETURN
#paste <(cat $TEMP | tail -n 1 | awk {'print $2'}) <(cat $RG | tail -n 1 | awk {'print $2'})>$WEST_PCOORD_RETURN
#cat $TEMP >pcoord.dat
# Clean up
#rm -f $TEMP md.in seg.nfo seg.pdb
