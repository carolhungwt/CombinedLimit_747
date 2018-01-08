#replaceReso.py
#!/bin/env/ python2

import argparse,os
channels=['4e','2e2mu','4mu']
files = ['clean_quad9_mH125','vbf_quad9_mH125','vh_quad9_mH125']

if __name__=="__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('-l','--loc',dest='loc',required=True)
	parser.add_argument('-a','--arg')
	args = parser.parse_args()
	tempdict,cmdstr={},''
	for chan in channels:
		resoloc = os.path.join(args.loc,'Reso_'+chan+'_quad9.txt')
		resoloc = os.path.abspath(resoloc)
		try:
			with open(resoloc,'r') as f:
				resostr = f.read()	
				tempdict[chan]=resostr
		except: 
			raise IOError('Cannot not find '+resoloc)
	for ifile in files:
		try:
			with open(ifile+'_tmpl.c','r') as fin:
				cfilestr = fin.read()
				for key,item in tempdict.iteritems():
					replstr = '<dcbPara_{chan}_2nd>'.format(chan=key)
#					print replstr
					cfilestr = cfilestr.replace(replstr,item)
				with open(ifile+'.c','w') as fout:
					fout.write(cfilestr)
		except:
			raise IOError('Cannot not find '+ifile+'_tmpl.c')
		cmdstr+=ifile+'.c '
	if args.arg:
		os.system('mkdir -p '+args.arg+'; cp '+cmdstr+' '+args.arg)
