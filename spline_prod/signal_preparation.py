#!/usr/bin/env python
#signal_preparation.py

def prepare_rpdfWS_script(var='mass'):
	if var == 'mass':
		ChosenVar = 'mean_pole'
		NotChosenVar = 'sigma_pole'
		varIndex = '0'
		NotChosenVarVal = 0.00407
	elif var == 'width':
		ChosenVar = 'sigma_pole'
		NotChosenVar = 'mean_pole'
		varIndex = '2'
		NotChosenVarVal = 125.
	else:
		raise Exception(var+" must be either 'mass' or 'width'")
	if not os.path.exists('make_rpdfWS_withCat_tmpl.c'):
		raise Exception('make_rpdfWS_withCat_tmpl.c does not exist')
	with open('make_rpdfWS_withCat_tmpl.c','r') as f:
		fstr = f.read()
		fstr = fstr.replace('<ChosenVar>',ChosenVar)
		fstr = fstr.replace('<NotChosenVar>',NotChosenVar)
		fstr = fstr.replace('<NotChosenVarVal>',str(NotChosenVarVal))
		fstr = fstr.replace('<varIndex>',varIndex)
		with open('make_rpdfWS_withCat.c','w') as fout:
			fout.write(fstr)

def submitjobs(outdir):
	if not os.path.isdir(outdir):
		raise Exception(outdir+' does not exist')
	os.system('bash submitjobs.sh {outdir}'.format(outdir=outdir))

def submitmakespline(outdir):
	if not os.path.isdir(outdir):
		raise Exception(outdir+' does not exist')
	if not os.path.exists('make_Sp_WS_withCat_all.sh'):
		raise Exception('make_Sp_WS_withCat_all.sh not found')
	os.system('bash make_Sp_WS_withCat_all.sh {outdir}'.format(outdir))


import os, argparse

if __name__=="__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('-v','--var',required=True,help='Chosen Variable: mass or width')
	parser.add_argument('-p','--process',action='append',help='Process 0: Prepare make_rpdfWS_withCat.c script 	1: submit jobs to make rpdfWS root files	 2: run make_ggH_spline_withCat')
	parser.add_argument('-o','--outdir', help='output dir of rpdfWS root files. Will locate all rpdfWS root files at \{outdir\}/rpdfWS_withCat', default=".")
	args=parser.parse_args()
	print args.process
	if '0' in args.process:
		prepare_rpdfWS_script(args.var)
	if '1' in args.process:
		submitjobs(args.outdir)
	if '2' in args.process:
		submitmakespline(args.outdir)


