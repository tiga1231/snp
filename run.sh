#!/bin/bash

FILES=./experiment/experiment_1*
for i in $FILES 
do
    #./gen.sh 
    a=${i:24}
    echo $a
    ./gen.sh $a
done
