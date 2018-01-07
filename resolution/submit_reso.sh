#!/bin/bash

#channels=(4e 2e2mu 4mu)
channels=(2e2mu 4mu 4e)
icat=0
quad=9
cur=$PWD
echo $cur
paramdir=$PWD
runswitch=$1
#/afs/cern.ch/user/w/wahung/work/public/CombineLimitDbkgkin/resolution_pkg/resolution_Quad9

#echo `awk 'BEGIN{FS="=";OFS=" ";} $0~/'${ch}'_cat'${cat}'=/ {print$2;}' massarray_quad9.txt `
#OR
#echo `awk 'BEGIN{OFS=" ";} $0~/4e_cat1=/ {$1="";print$0;}' massarray_quad9.txt `


dirs=( singleMassFit_Quad9 fixParamFit params )
for subdir in ${dirs[@]};do
if [ ! -d ${subdir} ]; then mkdir $subdir; fi
done

if [[ -z ${runswitch} || ${runswitch} -eq 1 ]];then
production=$2
#if [ ! -d "Quad9_tree" ];then bash rewriteTree.sh $quad ${production};fi
bash rewriteTree.sh $quad ${production}
fi

if [[ -z ${runswitch} || ${runswitch} -eq 2 ]]; then
for ch in ${channels[@]}; do
dmasscat=${icat}
while [ ${dmasscat} -lt $quad ]; do
	bash all_fitTF1_indfit.sh ${cur} ${ch} ${dmasscat}
((dmasscat++))
done;done
	python firstScreening.py
fi


if [[ -z ${runswitch} || ${runswitch} -eq 3 ]]; then
for ch in ${channels[@]}; do
dmasscat=${icat}
while [ ${dmasscat} -lt $quad ]; do
	#obtain the mass array from .txt, possible for fine-tuning 
#	rm params/individual_fit_param_${ch}_cat${dmasscat}.txt
	massarray=$(echo `awk 'BEGIN{FS="=";OFS=" ";} $0~/'${ch}'_cat'${dmasscat}'=/ {print$2;}' massarray_quad9.txt`)
	massstring=$(echo `awk 'BEGIN{OFS=",";} $0~/'${ch}'_cat'${dmasscat}'=/ {$1="";print$0;}' massarray_quad9.txt| cut -d',' -f2-`)
	echo ${massstring}
	# create customized readParam_single.sh for each Dmass Category
	sed -e "s~<massarray>~${massarray}~" -e "s~<paramdir>~${paramdir}~" -e "s~<MassStr>~${massstring}~" < readParam_single_tmpl_pol1.sh > readParam_single_${ch}_cat${dmasscat}.sh

	bash all_fitTF1.sh ${cur} ${ch} ${dmasscat}
	
((dmasscat++))	
done; done
fi
