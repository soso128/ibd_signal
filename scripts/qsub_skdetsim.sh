#!/bin/bash

# Parse command line arguments
while getopts ":dh" opt
do
    case $opt in
        h)
            echo '****************     help for qsub_skdetsim.sh    ****************:'
            echo 'Usage ./qsub_skdetsim <start run> <end run>'
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
    jobname=$run
    qsub -q lowe -o log/$jobname.out -e err/$jobname.err -r $jobname ./skdetsim.sh
done < ../runs.txt
