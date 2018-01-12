workdir=$1
cd $CMSSW_BASE/src
eval `scramv1 runtime -sh`
cd ${workdir}
ch=$2
vbf_cate=$3
onshell=1

quad=9
for folder in workspace125_onshell
do 	
	if [ -d "$folder" ];
	then 
		mkdir -p $folder
	fi
done
#for vbf_cate in 0 1 2; do
#for ch in 4e 4mu 2e2mu; do
	root -l -q -b clean_quad9_mH125.c\(\"${ch}\",${vbf_cate},${onshell},${quad}\)
	root -l -q -b vbf_quad9_mH125.c\(\"${ch}\",${vbf_cate},${onshell},${quad}\)
	root -l -q -b vh_quad9_mH125.c\(\"${ch}\",${vbf_cate},${onshell},${quad}\)
#done;done

