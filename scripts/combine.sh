#!/bin/bash

# Parse job name to get command line arguments
run=`echo $PJM_JOBNAME | cut -c2-`
tag=r$run

# Executables
combine=$PWD/../incorporate

# Inputs (vector files) and outputs
indir=/disk02/usr6/elhedri/SK2p2MeV/signal/srn/skdetsim/
infile=$indir/skdetsim.lowfit.$tag\.mcfit.root
outdir=/disk02/usr6/elhedri/SK2p2MeV/signal/srn/combined/
outfile=$outdir/combined.$tag.root
rm -f $outfile

# Check if input file exists
if [ ! -f $infile ]
then
    echo $infile" does not exist"
    exit 1
fi

# Combine events
time $combine $outfile $infile $run

echo ' Check: '
ls -haltr $outfile

