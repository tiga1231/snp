#!/bin/bash

#ibis -d experiment_1328 -q "select chr,start,stop,type,ref,alt where 1=1 limit 5"

ibis -d experiment/experiment_$1 -q "select chr,start,stop,type,ref,alt where 1=1 order by chr start limit 99999999999999" 2>result0.txt;
awk '
BEGIN{
    i=1;
}

{
    if(i>3 && NF==6){
        gsub("\"", ""); ## remove "
        gsub(" ",  "");## remove space
        print;
    }else{
        i++;
    }
}
' < result0.txt > data/$1.txt
rm result0.txt
