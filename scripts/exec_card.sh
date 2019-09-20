#!/bin/bash

# Cluster parameters
queue=all

# Executables
vector=./vectgen
SKDETSIM=/home/elhedri/skdetsim/skdetsim
lowfit=$PWD/../lowfit_sk4_gain_corr_mc
combine=$PWD/../incorporate

# Folders where to store results
vector_dir=/disk02/usr6/elhedri/SK2p2MeV/signal/srn/vector/
vector_prefix=vect
skdetsim_dir=/disk02/usr6/elhedri/SK2p2MeV/signal/srn/skdetsim/
lowfit_dir=/disk02/usr6/elhedri/SK2p2MeV/signal/srn/lowfit/
skdetsim_prefix=skdetsim
lowfit_prefix=skdetsim.lowfit
combine_dir=/disk02/usr6/elhedri/SK2p2MeV/signal/srn/combined/
combine_prefix=combined

# SKDetSim card
card_dir=/disk02/usr6/elhedri/SK2p2MeV/signal/srn/cards/
seed_dir=/disk02/usr6/elhedri/SK2p2MeV/signal/srn/seed/
