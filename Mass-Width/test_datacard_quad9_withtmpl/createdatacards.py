#!/usr/bin/python

import os, sys, argparse

class datacard():
	def __init__(self,tag,cat):
		self.tag = str(tag)
		self.cat = str(cat)

	def get_string(self, fdir):
		tag=self.tag
		cat = self.cat
		ftmpl = open('{1}/hzz4l_{0}_cat_tmpl.txt'.format(tag,fdir),'r')
		tmpl_string = ftmpl.read()
		return tmpl_string

	def read_rate(self):
		frate = open('fyield.txt','r')
		all_lines = frate.readlines()
		j=0		
		dmassstr = '{1}_Dmass_Cat{0}'.format(self.cat, self.tag)
		for line in all_lines:
			if dmassstr in line:
				index = j
				break
			j+=1
		return all_lines[index+1]

	def write_cat(self, tmpl_string):
		tmpl_string = tmpl_string.replace('<cat>',self.cat)
		return tmpl_string

	def write_rate(self, tmpl_string):
		rate_str = self.read_rate()
		tmpl_string = tmpl_string.replace('<rate>', rate_str)
		return tmpl_string

	def write_datacards(self,fdir,splineloc,bkgloc):
		tmpldir = '{0}/hzz4l_{1}_cat{2}_13TeV_onshell.txt'.format(fdir,self.tag,self.cat)
		#assert os.path.exists(tmpldir)
		fnew_card = open(tmpldir,'w')
		string_to_write = self.get_string(fdir)
		string_to_write = self.write_cat(string_to_write)
		string_to_write = self.write_rate(string_to_write)
		string_to_write = string_to_write.replace('<splinelocation>',splineloc)
		string_to_write = string_to_write.replace('<bkglocation>',bkgloc)
		fnew_card.write(string_to_write)


channels=["4e","4mu","2e2mu"]
#curdir = os.system("echo $PWD")
curdir = os.getcwd()
all_cards=[]
parser = argparse.ArgumentParser()
parser.add_argument('-s','--splineloc',required=True)
parser.add_argument('-b','--bkgloc',required=True)
args = parser.parse_args()
splineloc = os.path.abspath(args.splineloc)
bkgloc = os.path.abspath(args.bkgloc)
for chan in channels:
	for dcat in range(0,9):
		new_member = datacard(chan,dcat)
		new_member.write_datacards(curdir,splineloc,bkgloc)
		#all_cards.append(new_member)


