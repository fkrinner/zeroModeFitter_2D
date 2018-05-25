import pyRootPwa
import numpy as np

def parseTH1D(fileName, fakk = 1., addX = None):
	binning = []
	vals    = []
	errs    = []
	with open(fileName, 'r') as inin:
		for line in inin.readlines():
			parsed = [float(v) for v in line.split()]
			vals.append(parsed[2]*fakk)
			errs.append(parsed[3]*fakk)
			if len(binning) == 0:
				binning.append(parsed[0])
			else:
				if not parsed[0] == binning[-1]:
					raise ValueError("Binning is not contiguous")
			binning.append(parsed[1])
	if addX is not None:
		for i in range(len(binning)):
			binning[i] += addX

	binning = np.asarray(binning, dtype = np.float64)
	hist    = pyRootPwa.ROOT.TH1D('hist', 'hist', len(binning)-1, binning)
	for i in range(len(vals)):
		hist.SetBinContent(i+1, vals[i])
		hist.SetBinError(  i+1, errs[i])
	return hist

def parseArgand(fileName, skipZero = False, fakk = 1.):
	X  = []
	EX = []
	Y  = []
	EY = []

	with open(fileName, 'r') as inin:
		for line in inin.readlines():
			parsed = [float(v) for v in line.split()]
			if skipZero and parsed[0] == 0. and parsed[1] == 0. and parsed[2] == 0. and parsed[3] == 0.:
				continue
			X.append(parsed[0]*fakk)
			EX.append(parsed[1]*fakk)
			Y.append(parsed[2]*fakk)
			EY.append(parsed[3]*fakk)
	X  = np.asarray(X , dtype = np.float64)
	EX = np.asarray(EX, dtype = np.float64)
	Y  = np.asarray(Y , dtype = np.float64)
	EY = np.asarray(EY, dtype = np.float64)
	return X,EX, Y,EY

def parseTGraph(fileName, fakk = 1.):
	X,EX,Y,EX = parseArgand(fileName, False, fakk)
	graph = pyRootPwa.ROOT.TGraphErrors(len(X), X,Y,EX,EY)
	return graph

def main():
	fileName = "./compare.argand"
	hist     = parseTGraph(fileName)
	hist.Draw()
	raw_input()

if __name__ == "__main__":
	main()
	
