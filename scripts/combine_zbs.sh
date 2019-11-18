#!/bin/bash

source exec_card.sh

# Parse job name to get command line arguments
args=(${PJM_JOBNAME//_/ })
run=`echo ${args[0]} | cut -c2-`
fnum1=`printf "%04d" ${args[1]}`
#findex=${args[2]}
seed=`./genrand.sh`

#fstart=`expr $findex \* 9`
#ifinal=8
#if [ "$fstart" -eq 4 ]
#then
    #ifinal=9
#fi
for i in `seq 0 45`
do
    # Inputs (vector files) and outputs
    fnum2=`printf "%02d" $i`
    inname=$lowfitzbs_dir/
    infile=$inname/$lowfitzbs_prefix\.$fnum1\_$fnum2\.zbs
    outname=$combinezbs_dir/
    mkdir -p $outname
    outfile=$outname/$lowfitzbs_prefix\.$fnum1\_$fnum2\.root
    rm -f $outfile

    # Check if input file exists
    if [ ! -f $infile ]
    then
        echo $infile" does not exist"
        exit 1
    fi
    i10=`expr $i \* 1000`
    nseed=`expr $seed + $i10`

    # Combine events
    echo "$combine $outfile $infile $run $nseed y"
    time $combine $outfile $infile $run $nseed y

    echo ' Check: '
    ls -haltr $outfile
done

