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

bash readParam_single_${ch}_cat${cat}.sh $ch $cat
root -l -n -b -q "plotVar_fitInd_${ch}_cat${cat}.C"
python checkParams.py 
python rewriteParamstxt.py $ch $cat
bash readParam_ind_pol1.sh $ch $cat
root -l -n -q -b readData.cc\(\"$ch\",${cat}\) simfit_${ch}_cat${cat}.cc\(\"$ch\",${cat}\)



