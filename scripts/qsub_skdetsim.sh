#!/bin/bash

source exec_card.sh

# Create log and error dirs
logdir=$lowfit_dir/log
errdir=$lowfit_dir/err
mkdir -p  $logdir
mkdir -p $errdir

# Parse command line arguments
while getopts ":dh" opt
do
    case $opt in
        h)
            echo '****************     help for qsub_skdetsim.sh    ****************:'
            echo 'Usage ./qsub_skdetsim <start run> <end run> <emin> <emax>'
            ;;
    esac
done

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
    jobname=$run\_$3\_$4
    qsub -q $queue -o $logdir/$jobname.out -e $logdir/$jobname.err -r $jobname ./skdetsim.sh
done < ../runs.txt
