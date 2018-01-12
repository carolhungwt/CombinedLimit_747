workdir=$PWD


for vbf_cate in 0 1 2; do
for ch in 4mu 2e2mu 4e; do

bsub -q 1nd write_clean.sh ${workdir} ${ch} ${vbf_cate}

done; done
cd $workdir
