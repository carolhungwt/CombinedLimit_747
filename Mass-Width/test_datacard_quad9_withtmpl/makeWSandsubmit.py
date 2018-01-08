massscan = ['m','mass','scan1d_clean_mass.sh']
widthscan = ['w','width','scan1d_clean.sh']

def getdate():
	now = datetime.datetime.now()
	year = str(now.year)
	if int(now.month)<10:
		month = '0'+str(now.month)
	else:	month = str(now.month)
	if int(now.day) <10:
		day = '0'+str(now.day)
	else:	day = str(now.day)
	return year+month+day

def texttows(txtcard, outroot):
	#text2workspace.py hzz4l_combined_13TeV_onshell.txt -o allch_quad9_withDmass_prod_170623.root -P HiggsAnalysis.CombinedLimit.HighmassModel:HighmassModel --PO=mwAsPOI -v 3
	cmd = 'text2workspace.py {card} -o {root} -P HiggsAnalysis.CombinedLimit.HighmassModel:HighmassModel --PO=mwAsPOI -v 3'.format(card=txtcard,root=outroot)
	os.system(cmd)

def submitjob(var,outroot):
	rootstem = ''
	if os.path.exists(outroot):
		rootstem = re.split('.root',outroot)
		rootstem = rootstem[0]		
	else:
		raise Exception( outroot+' production failed')
	cmd = 'bash '
	if var in massscan:
		foldertag, script = massscan[1],massscan[2]
	elif var in widthscan:
		foldertag, script = widthscan[1],widthmass[2]
	else:
		raise Exception('var not valid')
	cmd+='{script} {date}_{stem}_{tag} {outroot}'.format(script=script,stem=rootstem,tag=foldertag, outroot=outroot, date=getdate())
	try:
		os.system(cmd)
	except Exception as e:
		print type(e)
		print e.args
		print e
		

import argparse, os, re, datetime

if __name__=='__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-i','--inputcard',required=True,action='append')
	parser.add_argument('-o','--outputroot',required=True,action='append')
	parser.add_argument('-v','--var',required=True,help='mass(m) or width(w)')
	parser.add_argument('--process',default='all', type=str)
	args = parser.parse_args()

	if len(args.inputcard) != len(args.outputroot):
		raise Exception('length of inputcard and outputroot must match')
	for i in range(len(args.inputcard)):
		icard = (args.inputcard)[i]
		iroot = (args.outputroot)[i]
		if args.process == 'all' or args.process=='1':
			texttows(icard,iroot)
		if args.process =='all' or args.process == '2':
			submitjob(args.var,iroot)
