#!/bin/bash

runfile=runs_darkrate.txt

debug=0
spec=0

# Parse command line arguments
while getopts ":ds:hf:l:" opt
do
    case ${opt} in
        s)
            spec=1
            file=$OPTARG
            echo "Spectral file is "$file
            ;;
        f)
            frun=$OPTARG
            ;;
        l)
            lrun=$OPTARG
            ;;
        h)
            echo '****************     help for qsub_vectgen.sh    ****************:'
            echo 'Usage: ./qsub_vectgen.sh <start run> <end run> -s [file for spectral parameters]'
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
    if [ "$spec" -eq 0 ]
    then
        jobname=$run
    else
        jobname=$run\_$file
        echo $jobname
    fi
    seed=../seed/r$run
    echo $seed
    if [ -f $seed ]
    then
        rm -f $seed
    fi
    ../make_random $seed
    echo qsub -x -q lowe -o log/$run.vect.out -e err/$run.vect.err -r $jobname ./vectgen.sh
    qsub -x -q lowe -o log/$run.vect.out -e err/$run.vect.err -r $jobname ./vectgen.sh
done < ../runs.txt
