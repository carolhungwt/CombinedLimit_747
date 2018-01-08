#tag=$1
#!/bin/env bash
quad=9
if [ ! -d "spline_WS_withCat_quad${quad}" ]; then
mkdir -p spline_WS_withCat_quad${quad}
fi

for tag in 4e 4mu 2e2mu
do
cat=0
while [ $cat -lt 9 ]
do
#root -l -q -b loadLib.C make_rpdfWS_withCat.cc\(\"${tag}\",${cat}\)
root -l -q -b loadLib.C make_ggH_spline_withCat.cc++\(\"${tag}\",${cat},${quad}\)
((cat+=1))
done 
done
outdir=$PWD
#"/afs/cern.ch/user/w/wahung/work/CMSSW_8_0_26_patch1/src/HiggsMassConstraint_new/HiggsMassConstraint/test/rpdfWS_withCat/"
#mv *_rpdfWS_cat*.root $outdir
