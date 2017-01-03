import numpy as np
import numpy.linalg as la

numLim = 1.E-10

nF0   = 62
nRho  = 56
nF2   = 56
nRho3 = 56

def zeroForSectors(sectors, title):
	for sector in sectors:
		if sector in title:
			return True
	return False

def countSectors(hist):
	return len(getZeroHistSectors(hist))

def getZeroHistSectors(hist):
	return  hist.GetTitle().split("<|sec++>")[0].split("<|>")

def getNforSector(sector):
	if "[pi,pi]0++" in sector:
		return nF0
	elif "[pi,pi]1--" in sector:
		return nRho
	elif "[pi,pi]2++" in sector:
		return nF2
	elif "[pi,pi]3--" in sector:
		return nRho3
	else:
		raise ValueError("Could nto detect nBins for '" + sector + "'")

def getZeroHistBorders(zeroHist):
	sectors = getZeroHistSectors(zeroHist)
	borders = [0]
	for sector in sectors:
		borders.append(borders[-1] + getNforSector(sector))
	return borders

def renormToBinWidth(hist, exponent = 1.):
	for i in range(hist.GetNbinsY()):
		binWidth = hist.GetYaxis().GetBinWidth(i+1)
		for j in range(hist.GetNbinsX()):
			hist.SetBinContent(j+1,i+1, hist.GetBinContent(j+1,i+1)/binWidth**exponent)
			hist.SetBinError(j+1, i+1, hist.GetBinError(j+1, i+1)/binWidth**exponent)

def pinv(matrix, numLim = 1.e-13):
	dim = len(matrix)
	val, vec = la.eig(matrix)
	for i in range(dim):
		if abs(val[i]) < numLim:
			val[i] = 0.
		elif val[i].real < 0.:
			raise ValueError("Negative eingenvalue: " + str(val[i]))
		else:
			val[i] = 1./val[i]
	return np.dot(vec, np.dot(np.diag(val), np.transpose(vec)))

def isValidPhaseSpace(m3pi, m2pi, mPi = 0.139):
	if m3pi >= m2pi + mPi:
		return True
	return False
