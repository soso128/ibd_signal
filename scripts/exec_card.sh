#!/bin/bash

# Cluster parameters
queue=all
maxjobs=190

# Executables
vector=./vectgen
vector_run=./vectgen_run
SKDETSIM=/home/elhedri/skdetsim/skdetsim
lowfit=$PWD/../lowfit_sk4_gain_corr_mc
combine=$PWD/../incorporate
leaf=$PWD/../leaf
#leaf=/disk02/usr6/yashida/relic/mc/antinue/lowfit/leaf
apfit=/usr/local/sklib_g77/atmpd_14c/bin/apfit_fc.sh
root2zbs=/home/elhedri/relic/2ring/root2zbs/root2zbs
reweight=$PWD/../reweight

# Folders where to store results
vector_dir=/disk02/usr6/elhedri/SK2p2MeV/signal/srn/vector/
vector_prefix=vect
skdetsim_dir=/disk02/usr6/elhedri/SK2p2MeV/signal/srn/skdetsim/
#lowfit_dir=/disk02/lowe8/relic_sk4/mc/relic/lowfit/
lowfit_dir=/disk02/usr6/elhedri/SK2p2MeV/signal/srn/lowfit/
lowfitzbs_dir=/disk02/usr6/yashida/relic/mc/atmnu/sk4.apr16_lowfit_20191005/zbs/
skdetsim_prefix=skdetsim
lowfit_prefix=skdetsim.lowfit
lowfitzbs_prefix=lowfit.atmnu
combine_dir=/disk02/usr6/elhedri/SK2p2MeV/signal/srn/combined/
combinezbs_dir=/disk02/usr6/elhedri/SK2p2MeV/signal/srn/combinedzbs/
combine_prefix=combined
leaf_dir=/disk02/usr6/elhedri/SK2p2MeV/signal/srn/leaf/
leaf_prefix=leaf
apfit_dir=/disk02/usr6/elhedri/SK2p2MeV/signal/srn/apfit/
apfit_prefix=apfit
#weight_dir=/disk02/usr6/elhedri/SK2p2MeV/signal/srn/weighted_nakazato/
weight_dir=/disk02/usr6/elhedri/SK2p2MeV/background/reactors/
#weight_dir=/disk02/usr6/elhedri/SK2p2MeV/background/li9/
weight_prefix=reweight

# SKDetSim card
card_dir=/disk02/usr6/elhedri/SK2p2MeV/signal/srn/cards/
seed_dir=/disk02/usr6/elhedri/SK2p2MeV/signal/srn/seed/
