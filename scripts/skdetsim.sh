#!/bin/bash

# Parse job name to get command line arguments
args=(${PJM_JOBNAME//_/ })
run=`echo ${args[0]} | cut -c2-`
emin=${args[1]}
emax=${args[2]}
repeat=${args[3]}
tag=r$run.$repeat

source exec_card.sh

# Inputs (vector files) and outputs
inname=$vector_dir/$emin\_$emax
infile=$inname/$vector_prefix\.$tag.zbs
outname=$skdetsim_dir/$emin\_$emax
mkdir -p $outname
outfile=$outname/$skdetsim_prefix\.$tag.root
outname2=$lowfit_dir/$emin\_$emax
mkdir -p $outname2
outfile2=$outname2/$lowfit_prefix.$tag
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
RAN1=`./genrand.sh`
RAN2=`./genrand.sh`
echo $RAN1" "$RAN2

mkdir -p $card_dir
CARD=$card_dir/supersim.card.$PJM_JOBNAME
echo $CARD
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
# run lowe reconstruction without gain correction
time $lowfit 0 $outfile $outfile2

echo ' Check: '
ls -haltr $outfile

