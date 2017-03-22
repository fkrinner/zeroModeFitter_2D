import pyRootPwa
import numpy as np

def parseTH1D(fileName):
	binning = []
	vals    = []
	errs    = []
	with open(fileName, 'r') as inin:
		for line in inin.readlines():
			parsed = [float(v) for v in line.split()]
			vals.append(parsed[2])
			errs.append(parsed[3])
			if len(binning) == 0:
				binning.append(parsed[0])
			else:
				if not parsed[0] == binning[-1]:
					raise ValueError("Binning is not continuous")
			binning.append(parsed[1])
	binning = np.asarray(binning, dtype = np.float64)
	hist    = pyRootPwa.ROOT.TH1D('hist', 'hist', len(binning)-1, binning)
	for i in range(len(vals)):
		hist.SetBinContent(i+1, vals[i])
		hist.SetBinError(  i+1, errs[i])
	return hist

def parseArgand(fileName, skipZero = False):
	X  = []
	EX = []
	Y  = []
	EY = []
	with open(fileName, 'r') as inin:
		for line in inin.readlines():
			parsed = [float(v) for v in line.split()]
			if skipZero and parsed[0] == 0. and parsed[1] == 0. and parsed[2] == 0. and parsed[3] == 0.:
				continue
			X.append(parsed[0])
			EX.append(parsed[1])
			Y.append(parsed[2])
			EY.append(parsed[3])
	X  = np.asarray(X , dtype = np.float64)
	EX = np.asarray(EX, dtype = np.float64)
	Y  = np.asarray(Y , dtype = np.float64)
	EY = np.asarray(EY, dtype = np.float64)
	return X,EX, Y,EY

def parseTGraph(fileName):
	X,EX,Y,EX = parseArgand(fileName, False)
	graph = pyRootPwa.ROOT.TGraphErrors(len(X), X,Y,EX,EY)
	return graph

def main():
	fileName = "./compare.argand"
	hist     = parseTGraph(fileName)
	hist.Draw()
	raw_input()

if __name__ == "__main__":
	main()
	
