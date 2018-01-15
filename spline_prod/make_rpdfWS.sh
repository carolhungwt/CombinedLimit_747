#!/bin/env bash
outdir=$1
quad=9
workdir=$3
cd $CMSSW_BASE/src
eval `scramv1 runtime -sh`
cd ${workdir}
#cp  ${workdir}/make_rpdfWS_withCat.c .

cd ${workdir}

if [ ! -d "${outdir}/rpdfWS_withCat" ];
then
mkdir -p ${outdir}/rpdfWS_withCat
fi

tag=$2

root -l -q -b -n make_rpdfWS_withCat.c\(\"${workdir}\",\"${tag}\",${quad}\)
