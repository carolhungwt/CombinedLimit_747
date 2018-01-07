#!/bin/env python2

"""individual_fit_param_2e2mu_cat8.txt

2e2mu_cat8_a1_param_0: 1.605651E+00
2e2mu_cat8_a1_param_1: -5.523104E-03

2e2mu_cat8_a2_param_0: 1.408136E+00
2e2mu_cat8_a2_param_1: 5.220488E-03

2e2mu_cat8_mean_param_0: -3.224502E-01
2e2mu_cat8_mean_param_1: -5.690041E-03

2e2mu_cat8_n1_param_0: 1.304504E-01
2e2mu_cat8_n1_param_1: 3.369731E-02

2e2mu_cat8_n2_param_0: 1.157389E+00
2e2mu_cat8_n2_param_1: -5.850516E-06

2e2mu_cat8_sigma_param_0: 3.525836E-01
2e2mu_cat8_sigma_param_1: 2.459571E-02

SingleMassFit_ResoParam_MH125_2e2mu_Category_8.txt

a1_CB =  1.0695 +/- 0.086523 L(0 - 5)
a2_CB =  1.9864 +/- 0.22132 L(0 - 5)
mean_CB = -1.19997 +/- 0.075636 L(-50 - 1)
n1_CB =  3.2981 +/- 0.80768 L(0 - 5)
n2_CB =  2.4719 +/- 1.5765 L(0 - 20)
sigma_CB =  3.5865 +/- 0.073318 L(0 - 500)
"""

chans = ['2e2mu', '4e', '4mu']
maxCat=9
jsonfile='failed_categ'

def backupOldFile(fileloc):
	os.system('mkdir -p ./params_backup')
	os.system('cp '+fileloc+' ./params_backup')
	fileloc = re.split('./params/|\n',fileloc)
	if os.path.exists('./params_backup/'+fileloc[0]):
		return True

class paramTxt(object):
	def __init__(self, chan, cat):
		self.chan = chan
		self.cat = cat
		self.fileloc = self._getfileloc(chan, cat)
		self._mh125paramsdict = self._readmh125param()

	def _getfileloc(self, chan, cat):
		fileloc = './params/individual_fit_param_'+chan+'_cat'+str(cat)+'.txt'
		if not os.path.exists(fileloc):	return 0
		else:
			isbackup = backupOldFile(fileloc)
			if isbackup:		return fileloc
			else:				raise Exception(fileloc+' not backed up yet')

	def _readmh125param(self):
		tempval,tempdict = 0,{}
		with open('massarray_quad9.txt','r') as fread:	
			lines = fread.readlines()
			headerkey = '{chan}_cat{cat}='.format(chan=self.chan,cat=self.cat)
			for line in lines:
				if headerkey in line:
					line  = line.split()
					masses = line[1:]
			for i,mass in enumerate(masses):
				mh125fileloc = 'singleMassFit_Quad{maxCat}/SingleMassFit_ResoParam_MH{mass}_{chan}_Category_{cat}.txt'.format(mass=mass,maxCat=maxCat, chan=self.chan, cat=self.cat)
				if not os.path.exists(mh125fileloc):
					raise Exception(mh125fileloc+' not found')
				with open(mh125fileloc,'r') as f:
					lines = f.readlines()
					for line in lines:
						line = line.split()
						if line and i==0:
							tempdict[line[0]]=float(line[2])
						elif line:
							tempdict[line[0]]+=float(line[2])
						else:
							continue
			for key,item in tempdict.iteritems():
				tempdict[key]=item/len(masses)
			return tempdict

	def readAndreplace(self, param):
		key = param+'_CB'
		tobereplaced={}
		stem='{chan}_cat{cat}_{param}_param'.format(chan=self.chan, cat=self.cat, param=param)
		with open(self.fileloc) as f:
			tempfileloc = self.fileloc+'_temp'
			with open(tempfileloc,'w+') as fout:
				done=False
				first=0
				lines = f.readlines()
				for line in lines:
					line = line.split()
					if line  and not done and stem in line[0]:
						if not first:
							val = self._mh125paramsdict[key]	
							linetowrite = '{stem}_{first}:  {value}'.format(stem=stem, first=first, value=val)
							first=1
						elif first:
							val=0
							done = True
							linetowrite = '{stem}_{first}:  {value}'.format(stem=stem, first=first, value=val)
					else:	
						linetowrite = ' '.join(line)
					linetowrite=linetowrite+'\n'	
					fout.write(linetowrite)
		os.system('mv '+tempfileloc+' '+self.fileloc)
	

def readjson():
	data=0
	try:
		with open(jsonfile+'.json','r') as f:
			data = json.load(f)
	except:
		raise IOError(jsonfile+'.json not found')
	else:
		return data


import re, os, json, sys
if __name__=='__main__':
	data = readjson()
	curchan = str(sys.argv[1])
	curcat = str(sys.argv[2])
#	print data
	for chan, subdict in data.iteritems():	
		for cat, paramslist in subdict.iteritems():	
			if not (chan == curchan and curcat == cat):	continue
			tempparamTxt = paramTxt(chan, cat)
			for param in paramslist:
				tempparamTxt.readAndreplace(param)	
			print ' '.join(paramslist)+' CB params for {chan} category {cat} have been rewritten'.format(chan=chan,cat=cat)	
