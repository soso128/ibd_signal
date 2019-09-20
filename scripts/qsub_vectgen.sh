#!/bin/bash

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

# Parse command line arguments
while getopts ":dhf:l:e:E:" opt
do
    case ${opt} in
        e)
            emin=$OPTARG
            ;;
        E)
            emax=$OPTARG
            ;;
        f)
            frun=$OPTARG
            ;;
        l)
            lrun=$OPTARG
            ;;
        h)
            echo '****************     help for qsub_vectgen.sh    ****************:'
            echo 'Usage: ./qsub_vectgen.sh <start run> <end run> <emin> <emax>'
            exit 0
            ;;
    esac
done
            

while read run
do
    if [ "$run" -lt $frun ]
    then
        continue
    fi
    if [ "$run" -gt $lrun ]
    then
        break
    fi
    jobname=$run\_$emin\_$emax
    echo $jobname
    mkdir -p $seed_dir/$emin\_$emax
    seed=$seed_dir/$emin\_$emax/r$run
    echo $seed
    if [ -f $seed ]
    then
        rm -f $seed
    fi
    ../make_random $seed
    echo qsub -x -q lowe -o $logdir/$run.vect.out -e $errdir/$run.vect.err -r $jobname ./vectgen.sh
    qsub -x -q $queue -o $logdir/$run.vect.out -e $errdir/$run.vect.err -r $jobname ./vectgen.sh
done < ../runs.txt
