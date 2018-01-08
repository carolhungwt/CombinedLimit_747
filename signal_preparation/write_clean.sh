onshell=1

quad=9
for folder in workspace125_onshell
do 	
	if [ -d "$folder" ];
	then 
		mkdir -p $folder
	fi
done
for vbf_cate in 0 1 2
do
for ch in 4e
do
	root -l -q -b clean_quad9_mH125.c\(\"${ch}\",${vbf_cate},${onshell},${quad}\)
	root -l -q -b vbf_quad9_mH125.c\(\"${ch}\",${vbf_cate},${onshell},${quad}\)
	root -l -q -b vh_quad9_mH125.c\(\"${ch}\",${vbf_cate},${onshell},${quad}\)
done
done
