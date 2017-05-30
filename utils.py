import numpy as np
import numpy.linalg as la

numLim = 1.E-10
INF = float("inf")

def getNDp(sector):
	"""
	Gets the number of freed isobar bins for D decays
	"""
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
	"""
	Determines, whether a zero mode cintributes to the sector at hand
	"""
	for sector in sectors:
		if sector in title:
			return True
	return False

def countSectors(hist):
	"""
	Counts the sectors, a zero mode is contributung to
	"""
	return len(getZeroHistSectors(hist))

def getZeroHistSectors(hist):
	"""
	Gets the list of sectors for a zero mode
	"""
	return  hist.GetTitle().split("<|sec++>")[0].split("<|>")

def getNforSector(sector):
	"""
	Gets the number of bind for a sctor
	"""
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
	"""
	Gets the sector borders for a zero mode histogram
	"""
	sectors = getZeroHistSectors(zeroHist)
	borders = [0]
	for sector in sectors:
		borders.append(borders[-1] + getNforSector(sector))
	return borders

def renormToBinWidth(hist, exponent = 1.):
	"""
	Renorms a histogram to the bin width (with a certain exponent)
	"""
	for i in range(hist.GetNbinsY()):
		binWidth = hist.GetYaxis().GetBinWidth(i+1)
		for j in range(hist.GetNbinsX()):
			hist.SetBinContent(j+1,i+1, hist.GetBinContent(j+1,i+1)/binWidth**exponent)
			hist.SetBinError(j+1, i+1, hist.GetBinError(j+1, i+1)/binWidth**exponent)

def pinv(matrix, numLim = 1.e-13):
	"""
	Own method fot the pseudo inverse of a matrix
	"""	
	dim = len(matrix)
	val, vec = la.eig(matrix)
	for i in range(dim):
		if abs(val[i]) < numLim:
			val[i] = 0.
		elif val[i].real < -2.e-8:
			raise ValueError("Negative eingenvalue: " + str(val[i]))
		else:
			val[i] = 1./val[i]
	return np.dot(vec, np.dot(np.diag(val), np.transpose(vec)))

def isValidPhaseSpace(m3pi, m2pi, mPi = 0.139):
	"""
	Ckecks, if the process is kinametically valid
	"""
	if m3pi >= m2pi + mPi:
		return True
	return False

def LtoInt(L):
	"""
	Gets the spin quantum number for the corresponding letter
	"""
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
	if L == 'H':
		return 5


def getOnePhaseDirection(ampl):
	"""
	Returns the direction of a phase roatation of the complex inpit number
	"""
	return (-ampl.imag, ampl.real)

def normVector(vector):
	"""
	Normalizes a the vector to 1
	"""
	dim  = len(vector)
	norm = 0.
	for v in vector:
		norm += v**2
	norm**=.5
	for i in range(dim):
		vector[i] /= norm
	return vector

def getPhaseDirection(ampls):
	"""
	Gets the normalized direction of a global phase rotation in all complex numbers
	"""
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
	"""
	Returns the ID number of a zero mode histograms
	"""
	name = hist.GetName()
	if not name.startswith('zero'):
		raise NameError("'" + name + "' is not the name of a zeroMode histogram")
	return int(name.split('_')[0][4:])

def sumUp(ll):
	"""
	Sums up the scalar content of an n-dimensional array
	"""
	if hasattr(ll, '__len__'):
		summ = 0.
		for l in ll:
			summ += sumUp(l)
		return summ
	return ll

def cloneZeros(lst):
	"""
	Clones the shape of an array, filled with zeros
	"""
	retVal = []
	for val in lst:
		if hasattr(val, '__len__'):
			retVal.append(cloneZeros(val))
		else:
			retVal.append(0.)
	return retVal

def addBtoA(A,B, weight = 1.):
	"""
	Adds the entries od array B at the corresponding palces of the equal shaped array A, unsing the weight as factor
	"""
	if not len(A) == len(B):
		raise ValueError("addBtoA(...): Size mismatch ('" + str(len(A)) + "' != '" + str(len(B)) + "')" )
	for i in range(len(A)):
		if hasattr(A[i], '__len__'):
			addBtoA(A[i],B[i], weight)
		else:
			A[i] += weight*B[i]

def divideAbyB(A,B):
	"""
	Divides the entries of A by the corresponding entris of B
	"""
	if not len(A) == len(B):
		raise ValueError("divideAbyB(...): Size mismatch ('" + str(len(A)) + "' != '" + str(len(B)) + "')" )
	for i in range(len(A)):
		if hasattr(A[i], "__len__"):
			divideAbyB(A[i], B[i])
		else:
			A[i] /= B[i]

def weightedSum(weights, params):
	"""
	Calculates the weightd sum of the arrays
	"""
	if len(params) == 0:
		return []
	for m in weights:
		retVal = cloneZeros(params[m])
		break
	for m in weights:
		addBtoA(retVal, params[m], weights[m])
	return retVal

def invertAllEntries(A):
	for i in range(len(A)):
		if hasattr(A[i], "__len__"):
			invertAllEntries(A[i])
		else:
			A[i] = 1./A[i]

def main():
	lst = [[[1,2,3],[1,2]],[1,2,3]]
	for a in multiYield(lst):
		print a
	return

	zers =  cloneZeros(lst)
	weight = 2
	addBtoA(zers, lst, weight)
	addBtoA(zers, lst)
	print zers

if __name__ == "__main__":
	main()

	

