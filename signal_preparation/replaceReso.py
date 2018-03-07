#!/usr/bin/env python

import argparse,os
channels=['4e','2e2mu','4mu']
files = ['clean_quad9_mH125','vbf_quad9_mH125','vh_quad9_mH125']
withsysfiles = ['clean_quad9_mH125_withsys','vbf_quad9_mH125_withsys','vh_quad9_mH125_withsys']

def readparams(chan,ifile):
  addtag = ''
  if 'vh' in ifile:  addtag='_ZH'
  resoloc = os.path.join(args.loc,'Reso_{chan}_quad9{addtag}.txt'.format(chan=chan,addtag=addtag))
  resoloc = os.path.abspath(resoloc)
  try:
    with open(resoloc,'r') as f:
      resostr = f.read()      
    return resostr
  except:
    raise IOError("cannot locate "+resoloc)
  


if __name__=="__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('-l','--loc',dest='loc',required=True)
	parser.add_argument('-a','--arg')
	parser.add_argument('--withsys', default=True)
	args = parser.parse_args()
	tempdict,cmdstr={},''
	tmplcstr = '_withsys'
	if not args.withsys: tmplcstr=''
	for ifile in files:
	  try:
	    with open(ifile+tmplcstr+'_tmpl.c','r') as fin:	
		cfilestr = fin.read()	
		for chan in channels: 
		  resostr = readparams(chan,ifile)
		  replstr = '<dcbPara_{chan}_2nd>'.format(chan=chan)
		  assert replstr in cfilestr 
		  cfilestr = cfilestr.replace(replstr,resostr)
		with open(ifile+'.c','w') as fout:
		  fout.write(cfilestr)
	  except:
	    raise IOError('Cannot not find '+ifile+'_tmpl.c')
	  cmdstr+=ifile+'.c '
	  if args.arg:
	    os.system('mkdir -p '+args.arg+'; cp '+cmdstr+' '+args.arg)
