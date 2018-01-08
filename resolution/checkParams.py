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
"""

params = ['a1', 'a2', 'n1', 'n2', 'sigma']
chans = ['2e2mu', '4e', '4mu']
maxCat=9

def getParamDict(chan, cat):
	paramdict = {}
	fileloc = './params/individual_fit_param_'+chan+'_cat'+str(cat)+'.txt'
	if not os.path.exists(fileloc):	return 0
	with open(fileloc, 'r') as fin:
		lines = fin.readlines()
	for param in params:
		templist = []
		for line in lines:
			line = re.split('_|: |\n| ', line)	
			firstpar = 1
			if param in line:
				if firstpar!=-1:
					templist.append(float(line[-2]))
					if firstpar==1:		firstpar==0
					elif firstpar==0:	firstpar==-1	
		paramdict[param]=templist
	return paramdict

def checkBound(intercept, slope, bound):
	value = intercept+slope*bound
	if value<0:	
		print value
		return True	
	else:		return False

def checkParams(param, paramlist, ranges=[110,140]):
	intercept, slope = paramlist
	failed = False
	for bound in ranges:
		if not failed:
			failed = checkBound(intercept, slope, bound)
	return failed

def tojson(faileddict, filename):
	with open(filename+'.json','w+') as fout:
		json.dump(faileddict,fout)

import os, re
import json
if __name__=='__main__':
	faileddict = {}
	for chan in chans:
		channeedtorepl = False
		chandict = {}
		for cat in range(maxCat):
			failedlist = []
			needtoreplace = False
			paramdict = getParamDict(chan, cat)
			if paramdict==0:	continue
			for param in params:
				if checkParams(param, paramdict[param]):
					print chan, cat, param, paramdict[param]
					failedlist.append(param)
					needtoreplace = True
			if needtoreplace:	
				channeedtorepl = True
				key = str(cat)
				chandict[key]=failedlist
		if channeedtorepl:
			faileddict[chan] = chandict
	print faileddict
	tojson(faileddict, 'failed_categ')

