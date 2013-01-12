#!/bin/sh

CURRENT=`pwd`
BASENAME=`basename $CURRENT`
v_number=`echo $BASENAME | sed s/[a-zA-Z]*_//`;
#echo "version_number="$v_number;
echo $v_number;

exit;