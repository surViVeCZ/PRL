#!/bin/bash
#pocet cisel bud zadam nebo 10 :)
if [ $# -lt 1 ];then 
    numbers=10;
else
    numbers=$1;
fi;
dd if=/dev/random bs=1 count=$numbers of=numbers