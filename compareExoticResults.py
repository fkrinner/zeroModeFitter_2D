# compareExoticResults.py
import os, sys
import ROOT
from rootfabi import root_open
from cmath import phase

def parseDataFile(inFileName):
	data = []
	with open(inFileName, 'r') as inFile:
		for line in inFile.readlines():
			data.append([float(v) for v in line.split()])
	return data

def loadTxtFile(inFileName, mode = "intens"):
	data = parseDataFile(inFileName)
	intensVals = []
	for dp in data:
		m       = dp[ 0 ]
		ii      = dp[ 2 ]
		re      = dp[ 3 ]
		im      = dp[ 4 ]
		Er      = dp[ 5 ]
		Ei      = dp[ 6 ]
		cr      = dp[ 7 ]

		if mode == "intens":
			value  = re**2 + im**2
			uncer  = 2*(Er*re**2 + 2*cr*re*im + Ei*im**2)**.5
			value *= ii
			uncer *= ii		
		elif mode == "real":
			value = re
			uncer = Er**.5
			value *= ii**.5
			uncer *= ii**.5
		elif mode == "imag":
			value = im
			uncer = Et**.5
			value *= ii**.5
			uncer *= ii**.5
		elif mode == "phase":
			value = phase(re + 1.j*im)
			jac = [0.,0.]
			if re == 0.:
				if im > 0.:
					jac[0] = -1./im
				else:
					jac[0] =  1./im
			else:
				common = 1. + im**2/re**2
				jac[0] = -im/re**2/common
				jac[1] = 1./re/common
			uncer = (jac[0]*jac[0]*Er + jac[1]*jac[1]*Ei + 2*jac[1]*jac[0]*cr)**.5

		intensVals.append((m,value, uncer))
	return intensVals

def makeStdHistograms(data, name = "hist"):
	hist = ROOT.TH1D(name, name, 50, .5, 2.5)
	for val in data:
		bin = hist.GetXaxis().FindBin(val[0])
		hist.SetBinContent(bin, val[1])
		if len(val) > 2:
			hist.SetBinError(bin, val[2])
	return hist

def loadDataFromROOT(rootFileName, histName):
	data = []
	with root_open(rootFileName, "READ") as inFile:
		hist = inFile.Get(histName)
		if not hist:
			raise IOError("Could not get '" + histName + "' from '" + rootFileName + "'")

		for i in range(hist.GetNbinsX()):
			m = hist.GetXaxis().GetBinCenter(i+1)
			v = hist.GetBinContent(i+1)
			e = hist.GetBinError(i+1)
			data.append((m,v,e))
	return data

def main():
	inFileNameSameBase = "./1mp_rho_cpls_sameParams_<tBin>.dat"
	inFileNameAdjuBase = "./1mp_rho_cpls_adjustedParams_<tBin>.dat"

	freedFileName = "/nfs/mds/user/fkrinner/extensiveFreedIsobarStudies/results_exotic.root"
	fixedFileName = "/nfs/mds/user/fkrinner/extensiveFreedIsobarStudies/results_allIsob.root"
	elvenFileName = "/nfs/mds/user/fkrinner/extensiveFreedIsobarStudies/std11_withIntegral.root"

	freedHistName = "total_1-+1+[pi,pi]1--PiP_<tBin>"
	fixedHistName = "1-+1+rhoPiP_<tBin>_intens"

	mode = "phase"

	for tBin in range(4):
		intensesSame = loadTxtFile(inFileNameSameBase.replace("<tBin>", str(tBin)), mode = mode)
		intensesAdju = loadTxtFile(inFileNameAdjuBase.replace("<tBin>", str(tBin)), mode = mode)

		histSame = makeStdHistograms(intensesSame, "same"+str(tBin))
		histAdju = makeStdHistograms(intensesAdju, "adju"+str(tBin))

		histAdju.SetLineColor(2)

		if mode == "intens":
			intensesFixd = loadDataFromROOT(fixedFileName, fixedHistName.replace("<tBin>", str(tBin)))
			intensesFred = loadDataFromROOT(freedFileName, freedHistName.replace("<tBin>", str(tBin)))
			intensesSt11 = loadDataFromROOT(elvenFileName, fixedHistName.replace("<tBin>", str(tBin)))

			histFixd = makeStdHistograms(intensesFixd, "fixd"+str(tBin))
			histFred = makeStdHistograms(intensesFred, "fred"+str(tBin))
			histSt11 = makeStdHistograms(intensesSt11, "st11"+str(tBin))

			histFixd.SetLineColor(3)
			histFred.SetLineColor(4)
			histSt11.SetLineColor(5)

		histSame.Draw()
		histAdju.Draw("SAME")
		if mode == "intens":
			histFixd.Draw("SAME")
			histFred.Draw("SAME")
			histSt11.Draw("SAME")

		raw_input()

if __name__ == "__main__":
	main()
