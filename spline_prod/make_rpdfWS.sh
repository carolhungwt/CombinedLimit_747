#!/bin/env bash
outdir=$1
quad=9
workdir=$3
cp  ${workdir}/make_rpdfWS_withCat.c .

cd ${workdir}

if [ ! -d "${outdir}/rpdfWS_withCat" ];
then
mkdir -p ${outdir}/rpdfWS_withCat
fi

tag=$2
cat=$4
root -l -q -b -n make_rpdfWS_withCat.c\(\"${workdir}\",\"${tag}\",${cat},${quad}\)
