#!/bin/bash

args=(${PJM_JOBNAME//_/ })
run=`echo ${args[0]} | cut -c2-`
emin=${args[1]}
emax=${args[2]}

source exec_card.sh

date
outname=$vector_dir/$emin\_$emax
mkdir -p $outname
outfile=$outname/$vector_prefix\.r$run\.zbs
cd ..
echo "./vectgen "$run" "seed/r$run" "$outfile" "$emin" "$emax
$vector $run $seed_dir/$emin\_$emax/r$run $outfile $emin $emax
cd scripts/
date
