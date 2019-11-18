#!/bin/bash

args=(${PJM_JOBNAME//_/ })
run=`echo ${args[0]} | cut -c2-`
emin=${args[1]}
emax=${args[2]}
nevents=${args[3]}
repeat=${args[4]}
setpos=${args[5]}

source exec_card.sh

date
outname=$vector_dir/$emin\_$emax
mkdir -p $outname
outfile=$outname/$vector_prefix\.r$run\.$repeat\.zbs
cd ..
if [ "$setpos" -eq 0 ]
then
    echo "./vectgen_run "$run" "$seed_dir/r$run\_$repeat" "$outfile" "$emin" "$emax" "$nevents
    $vector_run $run $seed_dir/r$run\_$repeat\_$setpos $outfile $emin $emax $nevents
else
    echo "./vectgen_run $run $seed_dir/r$run\_$repeat $outfile $emin $emax $nevents $POSX $POSY $POSZ"
    $vector_run $run $seed_dir/r$run\_$repeat\_$setpos $outfile $emin $emax $nevents $POSX $POSY $POSZ
fi
cd scripts/
date
