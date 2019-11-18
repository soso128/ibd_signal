#!/bin/bash

args=(${PJM_JOBNAME//_/ })
run=`echo ${args[0]} | cut -c2-`
emin=${args[1]}
emax=${args[2]}
repeat=${args[3]}
setpos=${args[4]}

source exec_card.sh

date
outname=$vector_dir/$emin\_$emax
mkdir -p $outname
outfile=$outname/$vector_prefix\.r$run\.$repeat\.zbs
cd ..
if [ "$setpos" -eq 0 ]
then
    echo "./vectgen "$run" "$seed_dir/r$run\_$repeat" "$outfile" "$emin" "$emax
    $vector $run $seed_dir/r$run\_$repeat $outfile $emin $emax
else
    echo "./vectgen $run $seed_dir/r$run\_$repeat $outfile $emin $emax $POSX $POSY $POSZ"
    $vector $run $seed_dir/r$run\_$repeat $outfile $emin $emax $POSX $POSY $POSZ
fi
cd scripts/
date
