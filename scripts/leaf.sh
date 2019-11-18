#!/bin/bash

source exec_card.sh

# Parse job name to get command line arguments
args=(${PJM_JOBNAME//_/ })
run=`echo ${args[0]} | cut -c2-`
emin=${args[1]}
emax=${args[2]}
repeat=${args[3]}
tag=r$run.$repeat

# Inputs (vector files) and outputs
inname=$lowfit_dir/$emin\_$emax
infile=$inname/$lowfit_prefix\.$tag\.mcfit.root
outname=$leaf_dir/$emin\_$emax
mkdir -p $outname
outfile=$outname/$leaf_prefix\.$tag.root
rm -f $outfile

# Check if input file exists
if [ ! -f $infile ]
then
    echo $infile" does not exist"
    exit 1
fi

# leaf events
echo "$leaf $infile $outfile"
time $leaf $infile $outfile

echo ' Check: '
ls -haltr $outfile

