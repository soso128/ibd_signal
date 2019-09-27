#!/bin/bash

input=/disk02/usr6/elhedri/SK2p2MeV/signal/srn/lowfit/
out=/disk02/usr6/elhedri/SK2p2MeV/signal/srn/weighted_nakazato/
spectrum=/disk02/usr6/yashida/relic/mc/relic.flux/nakazato/nakazato_nuebar_ref_nh.txt

mkdir -p $out

emin=$1
emax=$2

input=$input/$emin\_$emax/

for f in $input/*
do
    fname=`basename $f .root`
    outfile=$out/$fname\_$emin\_$emax.root
    echo "Processing $f"
    ./read_reweight $outfile $spectrum $f
    ls -lh $outfile
done
