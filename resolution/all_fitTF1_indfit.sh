cur=$1
cur=$(readlink -f $cur)
echo "Current dir $cur"
ch=$2
echo "current channel ${ch}"
cat=$3
cd $CMSSW_BASE/src
eval `scramv1 runtime -sh`
cd $cur

rm plotVar_fitInd*cat*.C

#sed -e 's|<workdir>|'$cur'|g' -e 's|<maxCat>|'$quad'|g' < resofit_temp.h > resofit.h

root -l -n -q -b readData.cc\(\"$ch\",${cat}\) indfit.cc\(\"$ch\",${cat}\)



