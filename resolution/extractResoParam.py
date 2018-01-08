#!/bin/env python
import re
import os
import sys

paramspath = 'params'
npars=6
paramslist = ['a1','a2','mean','n1','n2','sigma']
channels=['4e','4mu','2e2mu']

def chechpositivity(p0,p1):
	positive= True
	if(float(p0)+float(p1)*110<0 or float(p0)+float(p1)*140<0):		positive=False
	return positive

def writeTotxt(chcard):
	with open('Reso_%s_quad%s.txt'%(chcard.channel,chcard.nquad),'w+') as fout:
		strtowrite = 'double dcbPara_%s_2nd[6][%s][2]={\n'%(chcard.channel, chcard.nquad)
		for param in paramslist:
			strtowrite+='{'
			for i,categ in enumerate(chcard):
				templist = categ.paramsdict[param]
				if(param != 'mean'):
					if(not chechpositivity(templist[0],templist[1])):	print chcard.channel,i,param
				if i != int(chcard.nquad)-1:	strtowrite+="{%s,%s},"%(templist[0],templist[1])
				else: 							strtowrite+="{%s,%s}"%(templist[0],templist[1])
			strtowrite+='},\n'
		strtowrite+='};\n'
		fout.write(strtowrite)

class categReso(object):
	def __init__(self,channel,dmasscat):
		self.channel = str(channel)
		self.dmasscat = int(dmasscat)
		self.paramsdict = self.fetchparams()

	def fetchparams(self):
		paramsdict={}
		for param in paramslist:
				paramsdict[param] = []
		paramstxt = '%s/individual_fit_param_%s_cat%s.txt'%(paramspath,self.channel,self.dmasscat)
		try:
			with open(paramstxt) as f:
				for line in f:
					if line.strip():	
						line = re.split('_|:|\n',line)	
						paramsdict[line[2]].append(line[-2])
				return paramsdict
		except IOError as e:
			print paramstxt + ' does not exist' 
			return e


class channelCard(object):
	def __init__(self, channel, nquad):
		self.nquad = nquad
		self.channel = channel 
		self.categlist = self.getReso()

	def getReso(self):
		categlist = []
		for i in range(self.nquad):
			tempcategcard = categReso(self.channel,i)
			categlist.append(tempcategcard)
		return categlist

	def __iter__(self):
		for categ in self.categlist:
			yield categ


if __name__=="__main__":
	for ch in channels:
		tempcard = channelCard(ch, 9)
		tempcard.getReso()
		writeTotxt(tempcard)
		

