#!/bin/env bash
workdir=$1
workdir=$( readlink -e ${workdir} )
cd $CMSSW_BASE/src
eval `scramv1 runtime -sh`
cd ${workdir}

inputdir=${workdir}"/rpdfWS_withCat"

#inputdir=${inputdir}"/rpdfWS_withCat"
quad=9
if [ ! -d "spline_WS_withCat_quad${quad}" ]; then
mkdir -p spline_WS_withCat_quad${quad}
fi


notready=true
echo $PWD
while [ ${notready} == true ]; do
numspline=$(find ./rpdfWS_withCat -maxdepth 1 -type f |wc -l)
if [ ${numspline} -ne 243 ]; then sleep 1m; else notready=false; fi
done

for tag in 4e 4mu 2e2mu
do
cat=0
while [ $cat -lt 9 ]
do
#root -l -q -b loadLib.C make_rpdfWS_withCat.cc\(\"${tag}\",${cat}\)
root -l -q -b -n make_ggH_spline_withCat.c\(\"${tag}\",${cat},${quad},\"${inputdir}\"\)
((cat+=1))
done 
done
outdir=$PWD


