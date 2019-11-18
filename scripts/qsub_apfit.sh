#!/bin/bash

ulimit -c 0

source exec_card.sh

# Create log and error dirs
logdir=$apfit_dir/log
errdir=$apfit_dir/err
mkdir -p  $logdir
mkdir -p $errdir

# Parse command line arguments
while getopts ":dh" opt
do
    case $opt in
        h)
            echo '****************     help for qsub_skdetsim.sh    ****************:'
            echo 'Usage ./qsub_apfit <start run> <end run> <emin> <emax>'
            exit 0
            ;;
    esac
done

repeat=$5

for run in `seq $1 $2`
do
    for r in `seq 1 $repeat`
    do
        inname=$lowfit_dir/$3\_$4
        infile=$inname/$lowfit_prefix\.r$run.$r.mcfit.root
        # Check if zbs file exists
        if [ ! -f $infile ]
        then
            echo $infile" does not exist"
            continue
        fi
        # Job limit
        jobrunning=`qstat -a $queue | grep $USER | wc -l`
        echo $jobrunning" jobs running"
        while [ $jobrunning -gt $maxjobs ]
        do
            jobrunning=`qstat -a $queue | grep $USER | wc -l`
        done
        jobname=$run\_$3\_$4\_$r
        echo "qsub -q $queue -o $logdir/$jobname.out -e $logdir/$jobname.err -r $jobname ./apfit.sh"
        qsub -q $queue -o $logdir/$jobname.out -e $errdir/$jobname.err -r $jobname ./apfit.sh
    done
done 

