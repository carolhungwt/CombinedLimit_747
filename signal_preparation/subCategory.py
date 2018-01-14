#subCategory.py
import ROOT as root

class subCategory(object):
	def __init__(self,channel,dmasscat,prodcate,vbfcate):
		self.channel = channel
		self.dmasscat = dmasscat
		self.prodcate = prodcate
		self.vbfcate = vbfcate

	def printYield(self):
		with open('yield.txt', 'a+') as f:
			printstr = '{ch}_cat{cat}{prodtag}_{vbf}jet:	{syield:>10}	{oldsyield:>10}		{oldmasssyield:>10}\n'.format(ch=self.channel, cat=self.dmasscat, prodtag=self.prodtag, vbf=self.vbfcate, syield=round(float(self.sigYield),5),oldsyield=round(self.getSigYield(0),5), oldmasssyield=round(self.getSigYield(1),5))
			f.write(printstr)

	def getSigYield(self, old=1):
		if old:		rootAddress = self.getOldAddress
		else:		rootAddress = self.getNewAddress
		f = root.TFile.Open(self.rootAddress)
		try:
			ws = f.Get("w")
			meanpole = ws.var('mean_pole')
			sigmapole = ws.var('sigma_pole')
			meanpole.setVal(125)
			sigmapole.setVal(0.004)
			sigintegral = ws.obj(self.getWidthIntegral)
			return float(sigintegral.getVal())
		except Exception as e:		
			print self.getNewAddress
			print self.getWidthIntegral
			print e

	def getSigYieldRootObj(self):
		f = root.TFile.Open(self.getOldAddress)
		ws = f.Get("w")
		sigintegral = ws.obj(self.getWidthIntegral)
		return sigintegral

	@property
	def getWidthIntegral(self):
	#vhsig_integral4e1_1jet
	#vbfsig_integral4e5_0jet
	#sig_integral4e2_2jetNominal	
		bodystem = "sig_integral"
		bodystem += self.channel+str(self.dmasscat)+"_"+str(self.vbfcate)+'jet'
		if self.prodcate == 0: 		bodystem += 'Nominal'
		elif self.prodcate == 1:	bodystem = 'vbf' + bodystem 
		elif self.prodcate == 2:	bodystem = 'vh' + bodystem
		return bodystem

	@property
	def prodtag(self):
		if self.prodcate==0:		return ""
		elif self.prodcate==1:		return "_vbf"
		elif self.prodcate==2:		return "_vh"

	@property 
	def getOldAddress(self):
		returnstr = "/eos/user/w/wahung/Mass_Width_Measurement/spline_WS_withCat_quad9_mass_only_171205_qqh/{chan}_spline_WS_cat{dmasscat}.root".format(chan=self.channel,dmasscat=self.dmasscat)
		return returnstr

	@property
	def getNewAddress(self):
		returnstr = "/afs/cern.ch/user/w/wahung/work/public/CombineLimitDbkgkin/CMSSW_7_4_7/src/HiggsAnalysis/CombinedLimit/signal_preparation/workspace125_onshell/hzz{chan}_13TeV.input_func{prodtag}_{vbfcate}jet_cat{dmasscat}.root".format(chan=self.channel, prodtag=self.prodtag, vbfcate=self.vbfcate, dmasscat=self.dmasscat)
		return returnstr

	@classmethod
	def getAllChannels(cls):
		for ch in ["4e", "2e2mu", "4mu"]:
			for dmasscat in range(9):
				for prodcate in range(3):
					for vbfcate in range(3):
						yield cls(ch, dmasscat, prodcate, vbfcate)


