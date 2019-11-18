#!/bin/bash

ulimit -c 0

source exec_card.sh

# Create log and error dirs
logdir=$leaf_dir/log
errdir=$leaf_dir/err
mkdir -p  $logdir
mkdir -p $errdir

# Parse command line arguments
while getopts ":dh" opt
do
    case $opt in
        h)
            echo '****************     help for qsub_leaf.sh    ****************:'
            echo 'Usage ./qsub_leaf.sh <start run> <end run> <emin> <emax>'
            ;;
    esac
done

repeat=$5

while read run
do
    if [ "$run" -lt $1 ]
    then
        continue
    fi
    if [ "$run" -gt $2 ]
    then
        break
    fi
    for r in `seq 1 $repeat`
    do
        # Job limit
        jobrunning=`qstat -a $queue | grep $USER | wc -l`
        echo $jobrunning" jobs running"
        while [ $jobrunning -gt $maxjobs ]
        do
            jobrunning=`qstat -a $queue | grep $USER | wc -l`
        done
        jobname=$1\_$3\_$4\_$r
        echo "qsub -q $queue -o $logdir/$jobname.out -e $errdir/$jobname.err -r $jobname ./leaf.sh"
        qsub -q $queue -o $logdir/$jobname.out -e $errdir/$jobname.err -r $jobname ./leaf.sh
    done
done < ../runs.txt
