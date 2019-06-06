#!/bin/bash

# Parse job name to get command line arguments
run=`echo $PJM_JOBNAME | cut -c2-`
tag=r$run

# Executables
SKDETSIM=/home/elhedri/skdetsim/skdetsim
lowfit=$PWD/../lowfit_sk4_gain_corr_mc

# Inputs (vector files) and outputs
dirv=/disk02/usr6/elhedri/SK2p2MeV/signal/srn/vector/
infile=$dirv/vect.$tag.zbs
outdir=/disk02/usr6/elhedri/SK2p2MeV/signal/srn/skdetsim/
outdir2=/disk02/usr6/elhedri/SK2p2MeV/signal/srn/lowfit/
outfile=$outdir/skdetsim.$tag.root
outfile2=$outdir/skdetsim.lowfit.$tag
rm -f $outfile
touch $outfile

# Check if zbs file exists
if [ ! -f $infile ]
then
    echo $infile" does not exist"
    exit 1
fi

# Write card

#set random seed
RAN1=`sh -c 'RANDOM=$0; echo $RANDOM' $1`
RAN2=`sh -c 'RANDOM=$0; echo $RANDOM' $2`

CARD=../card/supersim.card.$tag
if [ -f $CARD  ]
then
    rm $CARD
fi
sed -e 's/RAN1/'$RAN1'/; s/RAN2/'$RAN2'/;' ../supersim.card > $CARD

#card over

# run skdetsim
inputsize=$(wc -c <"$infile")
if [ "$inputsize" -lt 100 ]
then
    echo "Empty input file "$infile
    exit 1
fi
time $SKDETSIM $CARD $outfile $infile
# run lowe reconstruction with gain correction
time $lowfit 0 $outfile $outfile2

echo ' Check: '
ls -haltr $outfile

