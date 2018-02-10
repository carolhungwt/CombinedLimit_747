#outdir="/eos/user/w/wahung/Mass_Width_Measurement/171217_refit/width"
outdir=$1
if [ -z "${outdir}" ]; then echo "outdir not valid "; exit 1; fi
if [ ! -d ${outdir} ]; then echo "${outdir} is not a valid path"; fi
workdir=$PWD
quad=9
channels=(4e 4mu 2e2mu)

for ch in ${channels[@]}; do
bsub -q 1nd -J "make_rpdfWS"_${ch} make_rpdfWS.sh ${outdir} ${ch} ${workdir} 
done

donejobs=0
