massscan = ['m','mass','scan1d_clean_mass.sh']
widthscan = ['w','width','scan1d_clean.sh']

def texttows(txtcard, outroot):
	#text2workspace.py hzz4l_combined_13TeV_onshell.txt -o allch_quad9_withDmass_prod_170623.root -P HiggsAnalysis.CombinedLimit.HighmassModel:HighmassModel --PO=mwAsPOI -v 3
	cmd = 'text2workspace.py {card} -o {root} -P HiggsAnalysis.CombinedLimit.HighmassModel:HighmassModel --PO=mwAsPOI -v 3'.format(card=txtcard,root=outroot)
	os.system(cmd)

def submitjob(var,outroot):
	try:
		rootexists = os.path.exists(outroot)	
		rootstem = re.strip('.root'out,root)
		rootstem = rootstem[0]		
	except:
		print outroot+' production failed'
	cmd = 'bash '
	if var in massscan:
		foldertag, script = massscan[1],massscan[2]
	elif var in widthscan:
		foldertag, script = widthscan[1],widthmass[2]
	else:
		raise Exception('var not valid')
	cmd+='{script} {stem}_{tag} {outroot}'.format(script=script,stem=rootstem,tag=foldertag, outroot=outroot)
	try:
		os.system(cmd)
	except Exception as e:
		print type(e)
		print e.args
		print e
		

import argparse, os

if __name__=='__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-i','--inputcard',required=True)
	parser.add_argument('-o','--outputroot',required=True)
	parser.add_argument('-v','--var',required=True,help='mass(m) or width(w)')
	parser.add_argument('--process',default='all', type=str)
	args = parser.parser_args()

	if args.process == 'all' or args.process=='1':
		texttows(args.inputcard,args.outputroot)
	if args.process =='all' or args.process == '2':
		submitjob(args.var,args.outputroot)
