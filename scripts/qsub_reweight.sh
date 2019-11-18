#!/bin/bash

ulimit -c 0

source exec_card.sh

# Create log and error dirs
logdir=$weight_dir/log
errdir=$weight_dir/err
mkdir -p  $logdir
mkdir -p $errdir

# Parse command line arguments
while getopts ":dh" opt
do
    case $opt in
        h)
            echo '****************     help for qsub_reweight.sh    ****************:'
            echo 'Usage ./qsub_reweight.sh <start run> <end run> <emin> <emax> <repeat> <0 for antinu 1 for positron> <spectrum file>'
            exit 1
            ;;
    esac
done

repeat=$5
export SPECTRUM=$7

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
        jobname=$1\_$3\_$4\_$r\_$6
        echo "qsub -q lowe -o $logdir/$jobname.out -e $errdir/$jobname.err -r $jobname -x ./reweight.sh"
        qsub -q $queue -o $logdir/$jobname.out -e $errdir/$jobname.err -r $jobname -x ./reweight.sh
    done
done < ../runs.txt

