#!/bin/bash

ulimit -c 0

source exec_card.sh

# Create log and error dirs
logdir=$vector_dir/log
errdir=$vector_dir/err
mkdir -p  $logdir
mkdir -p $errdir

runfile=runs_darkrate.txt

debug=0
emin=-1
emax=-1
repeat=1

# Parse command line arguments
while getopts ":dhr" opt
do
    case ${opt} in
        h)
            echo '****************     help for qsub_vectgen_run.sh    ****************:'
            echo 'Usage: ./qsub_vectgen_run.sh <run> <emin> <emax> nevents repeat [posx] [posy] [posz]'
            exit 0
            ;;
    esac
done

run=$1
emin=$2
emax=$3
nevents=$4
repeat=$5

setpos=0
if [ "$#" -gt 5 ]
then
   export POSX=$6
   export POSY=$7
   export POSZ=$8
   setpos=1
fi
            

for r in `seq 1 $repeat`
do
    jobname=$run\_$emin\_$emax\_$nevents\_$r\_$setpos
    echo $jobname
    seed=$seed_dir/r$run\_$r\_$setpos
    echo $seed
    if [ -f $seed ]
    then
        echo "Keep former seed"
    else
        ../make_random $seed
    fi
    echo qsub -x -q $queue -o $logdir/$run.$r.vect.out -e $errdir/$run.$r.vect.err -r $jobname ./vectgen_run.sh
    qsub -x -q $queue -o $logdir/$run.$r.vect.out -e $errdir/$run.$r.vect.err -r $jobname ./vectgen_run.sh
    # Job limit
    jobrunning=`qstat -a $queue | grep $USER | wc -l`
    echo $jobrunning" jobs running"
    while [ $jobrunning -gt $maxjobs ]
    do
        jobrunning=`qstat -a $queue | grep $USER | wc -l`
    done
done
