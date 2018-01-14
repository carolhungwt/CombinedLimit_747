#groupSignalYields.py

import ROOT as root
from subCategory import subCategory

if __name__=="__main__":
	fout = root.TFile("all_sig_integral.root","recreate")
	w = root.RooWorkspace("w")
	for subcatobj in subCategory.getAllChannels():
		sigyieldobj = subcatobj.getSigYieldRootObj()
		getattr(w,'import')(sigyieldobj,root.RooFit.RecycleConflictNodes())
	fout.WriteTObject(w)
	fout.Close()
