#!/bin/bash

if [ -n "$SEG_DEBUG" ] ; then
  set -x
  env | sort
fi

cd $WEST_SIM_ROOT

tail -n +2 COLVAR | awk '{print $2}' > $WEST_PCOORD_RETURN

echo $WEST_PCOORD_RETURN
#rm $DIST

if [ -n "$SEG_DEBUG" ] ; then
  head -v $WEST_PCOORD_RETURN
fi
