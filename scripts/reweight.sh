#!/bin/bash

source exec_card.sh

# Parse job name to get command line arguments
args=(${PJM_JOBNAME//_/ })
run=`echo ${args[0]} | cut -c2-`
emin=${args[1]}
emax=${args[2]}
repeat=${args[3]}
isp=${args[4]}
tag=r$run.$repeat

# Inputs (vector files) and outputs
inname=$leaf_dir/$emin\_$emax
infile=$inname/$leaf_prefix\.$tag\.root
outname=$weight_dir/$emin\_$emax
mkdir -p $outname
outfile=$outname/$weight_prefix\.$tag.root
rm -f $outfile

# Check if input file exists
if [ ! -f $infile ]
then
    echo $infile" does not exist"
    exit 1
fi

# leaf events
echo "$reweight $isp $outfile $SPECTRUM $infile"
time $reweight $isp $outfile $SPECTRUM $infile

echo ' Check: '
ls -haltr $outfile

