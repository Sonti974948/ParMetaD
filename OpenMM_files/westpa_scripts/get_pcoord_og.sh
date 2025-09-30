#!/bin/bash

if [ -n "$SEG_DEBUG" ] ; then
  set -x
  env | sort
fi

cd $WEST_SIM_ROOT





#cat $DIST | tail -n +2 | awk '{print $2}' > $WEST_PCOORD_RETURN
tail -n +2 COLVAR | awk '{print $2}' > $WEST_PCOORD_RETURN
#rm $DIST

if [ -n "$SEG_DEBUG" ] ; then
  head -v $WEST_PCOORD_RETURN
fi
