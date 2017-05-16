import numpy as np
import numpy.linalg as la

numLim = 1.E-10
INF = float("inf")

def getNDp(sector):
	nDpfile = "/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/ppppppppp/build/nDp.dat"
	ns = []
	with open(nDpfile, 'r') as inin:
		for line in inin.readlines():
			ns += [int(chunk) for chunk in line.split()]
	if not len(ns) == 3:
		raise IOError("Dit not get exactly 3 values")
	if "0++" in sector:
		return ns[0]
	elif "1--" in sector:
		return ns[1]
	elif "2++" in sector:
		return ns[2]
	else:
		raise ValueError("Unknown sector '" + sector + "'")	

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
	if "Dp" in sector:
		return getNDp(sector)
	elif "[pi,pi]0++" in sector:
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

def LtoInt(L):
	if L == 'S':
		return 0
	if L == 'P':
		return 1
	if L == 'D':
		return 2
	if L == 'F':
		return 3
	if L == 'G':
		return 4

def getOnePhaseDirection(ampl):
	return (-ampl.imag, ampl.real)

def normVector(vector):
	dim  = len(vector)
	norm = 0.
	for v in vector:
		norm += v**2
	norm**=.5
	for i in range(dim):
		vector[i] /= norm
	return vector

def getPhaseDirection(ampls):
	dim = len(ampls)
	if dim%2 != 0:
		raise ValueError("Odd number of values... Error, cannot determine phase-direction")
	retVal = np.zeros(dim)
	for i in range(dim/2):
		pd     = getOnePhaseDirection(ampls[2*i] + 1.j * ampls[2*i+1])
		retVal[2*i  ] = pd[0]
		retVal[2*i+1] = pd[1]
	return normVector(retVal)

def getZeroModeNumber(hist):
	name = hist.GetName()
	if not name.startswith('zero'):
		raise NameError("'" + name + "' is not the name of a zeroMode histogram")
	return int(name.split('_')[0][4:])

def sumUp(ll):
	if hasattr(ll, '__len__'):
		summ = 0.
		for l in ll:
			summ += sumUp(l)
		return summ
	return ll

def cloneZeros(lst):
	retVal = []
	for val in lst:
		if hasattr(val, '__len__'):
			retVal.append(cloneZeros(val))
		else:
			retVal.append(0.)
	return retVal

def fillAwithB(A,B, weight = 1.):
	if not len(A) == len(B):
		raise ValueError("fillAwithB(...): Size mismatch ('" + str(len(A)) + "' != '" + str(len(B)) + "')" )
	for i in range(len(A)):
		if hasattr(A[i], '__len__'):
			fillAwithB(A[i],B[i], weight)
		else:
			A[i] += weight*B[i]

def weightedSum(weights, params):
	if len(params) == 0:
		return []
	for m in weights:
		retVal = cloneZeros(params[m])
		break
	for m in weights:
		fillAwithB(retVal, params[m], weights[m])
	return retVal

def main():
	lst = [[[1,2,3],[1,2]],[1,2,3]]
	zers =  cloneZeros(lst)
	weight = 2
	fillAwithB(zers, lst, weight)
	fillAwithB(zers, lst)
	print zers

if __name__ == "__main__":
	main()

	

