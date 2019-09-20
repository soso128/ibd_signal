#!/bin/bash

source exec_card.sh

# Create log and error dirs
logdir=$combine_dir/log
errdir=$combine_dir/err
mkdir -p  $logdir
mkdir -p $errdir

# Parse command line arguments
while getopts ":dh" opt
do
    case $opt in
        h)
            echo '****************     help for qsub_combine.sh    ****************:'
            echo 'Usage ./qsub_combine <start run> <end run> <emin> <emax>'
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
    echo "qsub -q lowe -o $logdir/$jobname.out -e $errdir/$jobname.err -r $jobname ./combine.sh"
    qsub -q $queue -o $logdir/$jobname.out -e $errdir/$jobname.err -r $jobname ./combine.sh
done < ../runs.txt
