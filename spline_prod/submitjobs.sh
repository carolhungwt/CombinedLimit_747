#outdir="/eos/user/w/wahung/Mass_Width_Measurement/171217_refit/width"
outdir=$1
if [ -z "${outdir}" ]; then echo "outdir not valid "; exit 1; fi
if [ ! -d ${outdir} ]; then echo "${outdir} is not a valid path"; fi
workdir=$PWD
quad=9

for ch in 4e 4mu 2e2mu; do
cat=0
while [ ${cat} -lt ${quad} ]; do
#bsub -q 1nd make_rpdfWS.sh ${outdir} ${ch} ${workdir} ${cat}
((cat++))
done; done 
