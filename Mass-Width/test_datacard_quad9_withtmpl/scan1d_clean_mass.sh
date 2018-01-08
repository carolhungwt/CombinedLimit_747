#!/bin/bash
folder=$1
root=$2
#workdir=$CMSSW_BASE/src/width
workdir=$PWD
#/afs/cern.ch/user/w/wahung/work/public/CMSSW_7_4_7/src/HiggsAnalysis/CombinedLimit/Mass-Width/prepareInputs_with_dbkgkin
#workdir=/afs/cern.ch/user/w/wahung/work/public/CMSSW_7_1_5/src/HiggsAnalysis/CombinedLimit/Mass-Width/prepareInputs
cd $CMSSW_BASE/src
eval `scramv1 runtime -sh`
cd $workdir
mkdir $folder
cd $folder
cp ../$root ../scan1d_clean_mass.lsf ./
#      for i in 29 46
     	for i in {0..1}
        do
          first=$(($i*1))
          last=$(($first+0))
          echo $first $last
					bsub -q 1nd -C 0 scan1d_clean_mass.lsf $first $last $i $folder $root $workdir
					echo $first $last $i
				done

cd $workdir
