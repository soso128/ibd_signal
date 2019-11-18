#!/bin/bash

# Parse job name to get command line arguments
args=(${PJM_JOBNAME//_/ })
run=`echo ${args[0]} | cut -c2-`
emin=${args[1]}
emax=${args[2]}
repeat=${args[3]}
tag=r$run.$repeat

source exec_card.sh

# Disable core dumps
ulimit -c 0

# Inputs (vector files) and outputs
inname=$lowfit_dir/$emin\_$emax
infile=$inname/$lowfit_prefix\.$tag.mcfit.root
infile_symb=aplinks/$lowfit_prefix\.$tag.root
ln -s $infile $infile_symb
outname=$apfit_dir/$emin\_$emax
mkdir -p $outname
outfile=$outname/$apfit_prefix\.$tag.zbs

# Check if input file exists
if [ ! -f $infile ]
then
    echo $infile" does not exist"
    exit 1
fi

# Convert to ZBS file
zbsfile=$outname/apfit_in.$tag.zbs
zbsfile_symb=aplinks/apfit_in.$tag.zbs
ln -s $zbsfile $zbsfile_symb
echo "$root2zbs $infile $zbsfile"
$root2zbs $infile_symb $zbsfile_symb
rm $infile_symb
echo 'Check ZBS'
ls -lh $zbsfile
echo "$apfit  $zbsfile $outfile"
$apfit  $zbsfile_symb $outfile
rm $zbsfile_symb
rm $zbsfile

echo ' Check: '
ls -haltr $outfile
