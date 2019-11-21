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
inname_zbs=$apfit_dir/$emin\_$emax
infile_zbs=$inname/$apfit_prefix\.$tag\.root
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
echo "$leaf $outfile $infile $infile_zbs"
time $leaf $outfile $infile $infile_zbs

echo ' Check: '
ls -haltr $outfile

