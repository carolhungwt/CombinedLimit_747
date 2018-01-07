#firstScreening.py
#!/bin/env/ python2

channels = ['4e','4mu','2e2mu']
maxCat = 9
masseslist=[115,120,124,125,126,130,135,140]

class firstScreenObj(object):

	def __init__(self, chan, cat, delmasses):
		self.delmasses = delmasses
		self.chan = chan
		self.cat = cat

	def MassesTxtStr(self):
		headerkey = '{chan}_cat{cat}='.format(chan=self.chan,cat=self.cat)
		masses = masseslist
		masses = map(str, masses)
		for delmass in self.delmasses:
			if delmass in masses:
				masses.remove(delmass)
		masses.insert(0,headerkey)
		returnstr = ' '.join(masses)
		returnstr+='\n'
		return returnstr


import itertools

if __name__=="__main__":
	deletedict = {}
	try:
		with open('checkchi2.txt') as f:
			lines = f.readlines()
			for line in lines:
				line = line.split()
				mass,chan,cat = str(line[0]), str(line[1]), str(line[2])
				tempdict_sub = {cat:[mass]}
				tempdict = {chan:tempdict_sub}
				if chan in deletedict.keys() and cat in deletedict[chan].keys():
					deletedict[chan][cat].append(mass)
				elif chan in deletedict.keys():
					deletedict[chan].update(tempdict_sub)
				else:
					deletedict.update(tempdict)
	except:
		raise IOError("checkchi2.txt does not exist")
	else:
		print deletedict
	with open('massarray_quad9.txt','w') as fout:
		for ichan, icat in itertools.product(channels,range(maxCat)):	
			try:
				fsolist = deletedict[ichan][str(icat)]
				#print fsolist
			except:
				fsolist = []
			tempobj = firstScreenObj(ichan,icat,fsolist)
			tempstr = tempobj.MassesTxtStr()
#			print tempstr
			fout.write(tempstr)

