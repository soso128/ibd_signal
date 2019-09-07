#!/bin/bash

args=(${PJM_JOBNAME//_/ })
run=`echo ${args[0]} | cut -c2-`
file=""
if [ ${#args[@]} -eq 2 ]
then
    file=${args[1]}
    file=spectrum_input/$file
    echo $file
fi

date
outfile=/disk02/usr6/elhedri/SK2p2MeV/signal/srn/vector/jul17/vect.r$run\.zbs
cd ..
echo "./vectgen "$run" "seed/r$run" "$outfile" "$file
./vectgen $run seed/r$run $outfile $file
cd scripts/
date
