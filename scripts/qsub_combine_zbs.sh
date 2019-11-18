#!/bin/bash

ulimit -c 0

source exec_card.sh

# Create log and error dirs
logdir=$combinezbs_dir/log
errdir=$combinezbs_dir/err
mkdir -p  $logdir
mkdir -p $errdir

# Run bins
runbin=(61525 65536 67554 68861 70324 71854 72898 73892 74743 76098)

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

repeat=$5
maxjobs=200
queue=all

for i in `seq 0 499`
do
    irun=`expr $i \/ 50`
    run=${runbin[$irun]}
    # Job limit
    jobrunning=`qstat -a $queue | grep $USER | wc -l`
    echo $jobrunning" jobs running"
    while [ $jobrunning -gt $maxjobs ]
    do
        jobrunning=`qstat -a $queue | grep $USER | wc -l`
    done
    jobname=$run\_$i
    echo "qsub -q $queue -o $logdir/$jobname.out -e $errdir/$jobname.err -r $jobname ./combine_zbs.sh"
    qsub -q $queue -o $logdir/$jobname.out -e $errdir/$jobname.err -r $jobname ./combine_zbs.sh
done
