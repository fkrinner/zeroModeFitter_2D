import numpy as np
import physUtils
from cmath import phase

class parameter:
	def __init__(self, value, name, error = 0.):
		self.name  = name
		self.value = value
		self.lock  = False
		self.error = error

	def __str__(self, nameWidth = 10):
		retVal = self.name
		while len(retVal) < nameWidth:
			retVal += ' '
		if len(retVal) > nameWidth:
			retVal = retVal[:nameWidth]
		retVal += ' = '
		retVal += "{:10.6f}".format(self.value) + " +- " + "{:10.6f}".format(self.error)
		return retVal

class PRparameter:
	def __init__(self, value, name, error = 0.):
		self.name  = name
		self.value = value
		self.lock  = False
		self.error = error

	def __str__(self, nameWidth = 10):
		retVal = self.name
		while len(retVal) < nameWidth:
			retVal += ' '
		if len(retVal) > nameWidth:
			retVal = retVal[:nameWidth]
		retVal += ' = '
		retVal += "{:10.6f}".format(self.value) + " +- " + "{:10.6f}".format(self.error)
		retVal += " (Interaction radius of {:10.6f} fm)".format(0.1973/self.value)

		return retVal

def loadAmplitudeFile(ampFileName):
	masses = []
	reals  = []
	imags  = []
	with open(ampFileName, 'r') as inin:
		for line in inin.readlines():
			chunks = line.split()
			masses.append(float(chunks[0]))
			reals.append(float(chunks[1]))
			imags.append(float(chunks[2]))
	return np.asarray(masses), np.asarray(reals), np.asarray(imags)

OMNES_LOADED = False
def loadOmnes():
	global OMNES_LOADED
	global S_OMNES
	global R_OMNES
	global I_OMNES
	omnesFileName = "/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/zeroModeFitter_2D/Omnes11_new.dat"
	if OMNES_LOADED:
		print "Omnes function already loaded"
	else:
		S_OMNES, R_OMNES, I_OMNES = loadAmplitudeFile(omnesFileName)
	OMNES_LOADED = True
	
def lookupOmnes(s):
	if not OMNES_LOADED:
		loadOmnes()
	for i in range(len(S_OMNES) - 1):
		if S_OMNES[i] > s and S_OMNES[i+1] >= s:
			x = (s - S_OMNES[i])/(S_OMNES[i+1] - S_OMNES[i])
			return (R_OMNES[i] + 1.j * I_OMNES[i])*(1-x) + (R_OMNES[i+1] + 1.j * I_OMNES[i+1])*x
	print "No lookup value found for s =",s
	return 0.+0.j

class parametersTrackingParameterization:
	def setParameters(self, params):
		if not len(params) == self.nPar:
			raise IndexError("Number of parameters does not match "+str(len(params))+" vs "+str(self.nPar))
		count = 0
		for i,v in enumerate(self.loadMap):
			if v:
				self.parameters[i].value = params[count]
				count += 1

	def makeLoadMap(self, parameters):
		self.nPar       = 0
		self.loadMap    = []
		self.parameters = parameters
		if not len(parameters) == self.nParAll:
			raise IndexError("Number of parameters does not match (" + str(len(parameters))+"/"+str(self.nParAll)+")")
		for p in parameters:
			if not p.lock:
				self.loadMap.append(True)
				self.nPar += 1
				p.lock = True
			else:
				self.loadMap.append(False)
#		print "Made loadMap:",str(self.loadMap)

	def setParametersAndErrors(self, params, errors):
		if not len(errors) == self.nPar:
			print "Number of errors does not match"
			return False
		if not len(params) == self.nPar:
			print "Number of parameters does not match"
			return False
		count = 0
		for i,v in enumerate(self.loadMap):
			if v:
				self.parameters[i].value = params[count]
				self.parameters[i].error = errors[count]
				count += 1
		return True

	def getParameters(self):
		retVal = []
		for i,v in enumerate(self.loadMap):
			if v:
				retVal.append(self.parameters[i].value)
		return retVal

	def returnParameters(self):
		return self.parameters

	def writeToFile(self, outFileName, masses, externalKinematicVariables = []):
		vals = self(masses, externalKinematicVariables)
		with open(outFileName, 'w') as outFile:
			for i in range(len(masses)):
				outFile.write(str(masses[i]) + ' ' + str(abs(vals[i])**2) + ' ' + str(vals[i].real) + ' ' + str(vals[i].imag) + ' ' + str(phase(vals[i])) + '\n')

class simpleOneChannelKmatrix(parametersTrackingParameterization):
	def __init__(self, parameters, nPole, polyDegPlusOne, sThresh, nComplex = 0):
		self.sThresh = sThresh
		self.nPole   = nPole
		self.nComx   = nComplex
		self.pdpo    = polyDegPlusOne
		self.nParAll = 2*nPole + 3*nComplex + polyDegPlusOne
		self.makeLoadMap(parameters)

		self.use_CM       = False
		self.secondSheet  = False
		self.iEpsilon     = 1.e-5j

	def __call__(self, ms, externalKinematicVariables = []):
		par = [p.value for p in self.parameters]
		retVals = np.zeros((len(ms)), dtype = complex)
		for i,m in enumerate(ms):
			s     = m**2 + self.iEpsilon

			if self.use_CM:
				iRhoPS = physUtils.threshChewMandelstam(s,self.sThresh)
			else:
				iRhoPS = -(self.sThresh/s-1.+0.j)**.5 # phase space
			K     = 0.
			count = 0
			for p in range(self.nPole):
				K     += par[count+1]/(par[count] - s)
				count += 2
			for p in range(self.nComx):
				K += par[count+2]/(par[count] + 1.j*par[count+1] - s)
				count += 3
			for p in range(self.pdpo):
				K     += par[count]*s**p
				count += 1
			retVals[i] = K/(1. - iRhoPS*K)
		return retVals

	def complexCall(self, s):
		par    = [p.value for p in self.parameters]
		if self.use_CM:
			iRhoPS = physUtils.threshChewMandelstam(s,self.sThresh)
		else:
			iRhoPS = -(self.sThresh/s-1.+0.j)**.5 # phase space
		K      = 0.
		count  = 0
		for p in range(self.nPole):
			K     += par[count+1]/(par[count] - s)
			count += 2
		for p in range(self.nComx):
				K += par[count+2]/(par[count] + 1.j*par[count+1] - s)
				count += 3
		for p in range(self.pdpo):
			K     += par[count]*s**p
			count += 1
		if not self.secondSheet:
			return K/(1.-iRhoPS*K)
		else:
#			iRhoPS = - (self.sThresh/s-1.+0.j)**.5
			return K/(1.- K*(iRhoPS+2*(self.sThresh/s-1.+0.j)**.5))

	def absInverse(self, reim):
		val = self.complexCall(reim[0] + 1.j*reim[1])
		return 1./abs(val)**2

	def makeMathematicaString(self, functionName = "K", var = 's'):
		par = [p.value for p in self.parameters]
		retVal = functionName + '[' + var + '_] := '
		count  = 0
		first  = True
		for p in range(self.nPole):
			if not first:
				retVal += " + "
			else:
				first = False
			retVal += '('+str(par[count+1].real)+")/("+str(par[count].real)+' - '+var+')'
			count  += 2
		for p in range(self.nComx):
			if not first:
				retVal += " + "
			else:
				first = False
			retVal += '('+str(par[count+2].real)+")/("+str(par[count].real) + ' +i ' + str(par[count+1]) + ' - ' +var + ')'
			count += 3
			
		for p in range(self.pdpo):
			if not first:
				retVal += " + "
			else:
				first = False
			retVal += '(' + str(par[count].real) + ')'
			count  += 1
			if p == 1:
				retVal += var
			if p > 1:
				retVal += var+'^'+str(p)
		rhoString = r"\[Rho]["+var+"_] := Sqrt[ 1 - " + str(self.sThresh)+'/'+var + ']'
		return retVal, rhoString

class binnedPolynomial(parametersTrackingParameterization):
	"""
	Independent polynomial for all bins in m3pi
	"""
	def __init__(self, parameters, mMin, mMax, nBins, degree, baseExponent, real = True):
		self.mMin  = mMin
		self.mMax  = mMax
		self.nBins = nBins	
		self.step  = (mMax-mMin)/nBins

		self.deg   = degree # 1 + c1 x + c2 x*2 is degree 2 (zertoh order is take care of by the normalization)
		self.base  = baseExponent
		self.real  = real

		self.parPerBin = self.deg
		if not self.real: # is complex valued
			self.parPerBin *= 2
		self.nParAll = self.nBins*self.parPerBin
		self.makeLoadMap(parameters)

	def __call__(self, ms, externalKinematicVariables = []):
		if len(externalKinematicVariables) == 0:
			raise ValueError("No 3pi mass bin given")
		nBin    = self.findBin(externalKinematicVariables[0])
		par     = [p.value for p in self.parameters]
		retVals = np.full((len(ms)),1., dtype = complex)
		cpls = [0.]*self.deg
		if self.real:
			for i in range(self.deg):
				cpls[i] = par[nBin*self.parPerBin + i]
		else:
			for i in range(self.deg):
				cpls[i] = par[nBin*self.parPerBin + 2*i] + 1.j*par[nBin*self.parPerBin + 2*i+1]
		for i,m in enumerate(ms):
			val = 1.
			for i in range(self.deg):
				val += cpls[i]*m**(self.base*(i+1))
			retVals[i] = val
		return retVals

	def findBin(self, m):
		if m < self.mMin or m > self.mMax:
			raise ValueError("Mass m = "+str(m)+" out of range: "+str(self.mMin) +" < m < "+str(self.mMax))
		return int((m-self.mMin)/self.step)

class monomial(parametersTrackingParameterization):
	def __init__(self, degree, parameters = []): # no base exponent, make the degre explicit here
		self.degree  = degree
		self.nParAll = 0
		self.makeLoadMap(parameters)

	def __call__(self, ms, externalKinematicVariables = []):
		retVals = np.full((len(ms)),1., dtype = complex)
		for i,m in enumerate(ms):
			retVals[i] = m**self.degree
		return retVals

class multiply(parametersTrackingParameterization):
	"""
	Functions to multiply may not be used elsewhere
	To do this, initialize the same function twice and use the same 
	instances of the 'parameters' class
	"""
	def __init__(self, functions):
		"""
		Do not need the paramters here!!!
		"""
		self.nPar       = 0
		self.nParAll    = 0
		self.functions  = functions
		self.borders    = [0]
		self.bordersAll = [0]
		for f in self.functions:
			self.nPar    += f.nPar
			self.borders.append(self.nPar)
			self.nParAll += f.nParAll
			self.bordersAll.append(self.nParAll)

	def setParameters(self, params):
		if not len(params) == self.nPar:
			raise IndexError("Number of parameters does not match")
		for i,f in enumerate(self.functions):
			f.setParameters(params[self.borders[i]:self.borders[i+1]])

	def setParametersAndErrors(self, params, errors):
		if not len(errors) == self.nPar:
			print "Number of errors does not match"
			return False
		if not len(params) == self.nPar:
			print "Number of parameters does not match"
			return False
		for i,f in enumerate(self.functions):
			f.setParametersAndErrors(params[self.borders[i]:self.borders[i+1]], errors[self.borders[i]:self.borders[i+1]])
		return True	

	def getParameters(self):
		retVal = []
		for f in self.functions:
			retVal += f.getParameters()
		return retVal

	def returnParameters(self):
		retVal = []
		for f in self.functions:
			retVal += f.parameters
		return retVal

	def __call__(self, ms, externalKinematicVariables = []):
		retVals = np.full((len(ms)),1., dtype = complex)
		for f in self.functions:
			retVals *= f(ms, externalKinematicVariables = externalKinematicVariables) # numpy arrays are cool
		return retVals

class productOfSimpleFunctions(parametersTrackingParameterization):
	def __init__(self, functions, parameters):
		self.nParAll = len(functions)
		self.fs      = functions
		self.makeLoadMap(parameters)

	def __call__(self,ms, externalKinematicVariables = []):
		par = [p.value for p in self.parameters]
		retVals = np.full((len(ms)),1., dtype = complex)
		for i,m in enumerate(ms):
			for nF,f in enumerate(self.fs):
				retVals[i] *= f(par[nF]*m)
		return retVals

class omnesFunctionPolynomial(parametersTrackingParameterization):
	"""
	Omnes polynomial allows to set parameters to real. For complex coefficients it it more efficient to use monomials
	"""
	def __init__(self, parameters, nDimPol = 0, shift = True, stretch = True, complexPolynomial = False):
		# Since a complex factor will be mutliplied in the fit, set c0 = 1., then:
		# c0 + c1 m + c2 m**2 + ... + cn m**n -> 1. + c1 m + ... + cn m**n
		self.complexPolynomial = complexPolynomial
		self.polDeg = nDimPol
		if complexPolynomial:
			self.nParAll = 2 * nDimPol
		else:
			self.nParAll = nDimPol
		self.shift = shift
		if shift:
			self.nParAll += 1
		self.stretch = stretch
		if stretch:
			self.nParAll += 1
		self.makeLoadMap(parameters)

	def __call__(self, ms, externalKinematicVariables = []):
		par = [p.value for p in self.parameters]
		if not len(par) == self.nParAll:
			raise IndexError("omnesFunctionPolynomial: Wrong number of parameters")
		retVals = np.zeros((len(ms)), dtype = complex)
		for i,m in enumerate(ms):
			skip = 0

			s = m**2
			x = s
#			print s,'+',
			if self.shift:
#				print par[skip],'=',
				s += par[skip]
#				print s,
				skip += 1
			if self.stretch:
#				print 
				s *= par[skip]
				skip += 1
			poly = 1.
			for d in range(self.polDeg):
				if self.complexPolynomial:
					poly += x * (par[2*d+skip] + 1.j*par[2*d+1+skip])
				else:
					poly += x * par[d+skip]
				x *= s
			retVals[i] = poly * lookupOmnes(s)
		return retVals

class fixedParameterization(parametersTrackingParameterization):
	def __init__(self, path):
		self.path       = path
		self.parameters = []
		self.nParAll    = 0
		self.makeLoadMap(parameters)
	
	def __call__(self, ms, externalKinematicVariables = []):
		retVals = np.zeros((len(ms)), dtype = complex)
		for i in range(len(retVals)):
			retVals[i] = 1.
		return retVals

class breitWigner(parametersTrackingParameterization):
	def __init__(self, parameters):
		self.nParAll    = 2
		self.makeLoadMap(parameters)

	def __call__(self, ms, externalKinematicVariables = []):
		par = [p.value for p in self.parameters]
		retVals = np.zeros((len(ms)), dtype = complex)
		num = par[0]*par[1]
		den = par[0]**2 - 1.j * par[0] * par[1]
		for i, m in enumerate(ms):
			retVals[i] = num/(den - m**2)
		return retVals

class pole(parametersTrackingParameterization):
	def __init__(self, parameters):
		self.nParAll = 2
		self.makeLoadMap(parameters)
		
	def __call__(self, ms, externalKinematicVariables = []):
		par = [p.value for p in self.parameters]
		retVals = np.zeros((len(ms)), dtype = complex)
		polPos = par[0] + 1.j*par[1]
		for i,m in enumerate(ms):
			s = m**2
			retVals[i] = 1./(polPos-s)
		return retVals

class pietarinenExpansion(parametersTrackingParameterization):
	def __init__(self, parameters, nPole, threshsolds = [], polDeg3Pi = 0):
		"""
		With an empty list of thresholds, the function is a linear combination of poles+1
		Since the fitter will multiply everything with a complex coefficient, set the 0th order in the expansions to 1.+0.j

		Paparameters are : [pole1ResuduumReal, pole1ResiduumImag, pole1PositionReal, pole1PositionImag, pole2..., thresh1Alpha, thres1C1, thresh1C2, ..., thres2...]
		"""
		self.nPole       = nPole
		self.threshsolds = threshsolds
		self.nThresh     = len(threshsolds)
		self.nParAll     = 4*self.nPole
		self.polDeg3Pi   = polDeg3Pi # Constant is polynomial of 0th degree... index handling more tedious, but what can you do
		for t in range(self.nThresh):
			if self.threshsolds[t][1] > 0:
				self.nParAll += 1+self.threshsolds[t][1]*(polDeg3Pi+1) # 1 for the alpha and self.threshsolds[t][1] for the polynomial coefficients
		self.makeLoadMap(parameters)

	def __call__(self, ms, externalKinematicVariables):
		retVals = np.zeros((len(ms)), dtype = complex)
		par = [p.value for p in self.parameters]
#		print ">>>",par,"<<<"
		if not len(par) == self.nParAll:
			raise IndexError("pietarinenExpansion.__call__(...): Parameter number mismatch")
		s3pi = externalKinematicVariables[0]**2
		for b,m in enumerate(ms):
			val      = 1.+0.j
			s        = m**2
			parCount = 0
			for _ in range(self.nPole):
				val      += (par[parCount] + 1.j*par[parCount+1])/(par[parCount+2] + 1.j*par[parCount+3] - s)
				parCount += 4
			for t in range(self.nThresh):
				nn = self.threshsolds[t][1]
				if nn == 0:
					continue
				sqrt = (self.threshsolds[t][0] - s+ 0.j)**.5 # + 0.j, that python handels complexity right. (First Riemann sheet will be used: (-1.+0.j)**.5 = 0.+1.j
				Z    = (par[parCount]-sqrt)/(par[parCount]+sqrt)
				parCount += 1
				for i in range(nn):
					C = 0.
					for d in range(self.polDeg3Pi+1):
						C += par[parCount+d]*s3pi**d
					val      += C*Z**(i+1)
					parCount += self.polDeg3Pi + 1
			retVals[b] = val
		if not parCount == self.nParAll:
			raise IndexError("pietarinenExpansion.__call__(...): Some parameter number mismatch:" + str(parCount) + " vs. " + str(self.nParAll))
#		print ">>>",retVals,"<<<"
		return retVals

globalLfactor = 1 # probably not needed anymore

class complexPolynomial(parametersTrackingParameterization):
	def __init__(self, degree, parameters, baseExponent = 1):
		self.degree = degree
		if not len(parameters) == 2*self.degree:
			raise IndexError("Number of parameters does not match")
		self.nParAll      = 2*degree

		self.baseExponent = baseExponent # Can be adjusted (2 e.g. gives a polynomial in s = m**2)
		self.makeLoadMap(parameters)

	def __call__(self,ms, externalKinematicVariables):
		retVals = np.zeros((len(ms)), dtype = complex)
		par = [p.value for p in self.parameters]
		for i, m in enumerate(ms):
			val     = 1.+0.j
			massExp = 1.
			for j in range(self.degree):
				massExp *= m**self.baseExponent
				val     += (par[2*j] + 1.j*par[2*j+1]) * massExp
			retVals[i] = val
		return retVals

class realPolynomial(parametersTrackingParameterization):
	def __init__(self, degree, parameters, baseExponent = 1):
		self.degree = degree
		if not len(parameters) == self.degree:
			raise IndexError("Number of parameters does not match")
		self.nParAll      = degree
		self.baseExponent = baseExponent # Can be adjusted (2 e.g. gives a polynomial in s = m**2)
		self.makeLoadMap(parameters)

	def __call__(self, ms, externalKinematicVariables = []):
		retVals = np.zeros((len(ms)), dtype = complex)
		par = [p.value for p in self.parameters]
		for i, m in enumerate(ms):
			val     = 1.
			massExp = 1.
			for j in range(self.degree):
				massExp *= m**self.baseExponent
				val     += par[j]*massExp
			retVals[i] = val
		return retVals

class twoDimensionalRealPolynomial(parametersTrackingParameterization):
	def __init__(self, degree2Pi, degree3Pi, parameters, baseExponent = 1):
		self.degree2Pi    = degree2Pi
		self.degree3Pi    = degree3Pi
		self.baseExponent = baseExponent
		self.nParAll      = degree2Pi*(degree3Pi+1)
		self.makeLoadMap(parameters)

	def __call__(self, ms, externalKinematicVariables = []):
		retVals = np.zeros((len(ms)), dtype = complex)
		par     = [p.value for p in self.parameters]
		m3      = externalKinematicVariables[0]
		for i, m in enumerate(ms):
			val   = 1.
			count = 0
			for i2 in range(self.degree2Pi):
				c = 0.
				for i3 in range(self.degree3Pi+1):
					c     += par[count]*m3**(i3*self.baseExponent)
					count += 1
				val += c*m**((i2+1)*self.baseExponent)			
			retVals[i] = val
		return retVals

class relativisticBreitWigner(parametersTrackingParameterization):
	def __init__(self, parameters, m1, m2, m3, J, L, fitPr = False, barrierFactors = True):
		self.L            = L * globalLfactor
		self.J            = J * globalLfactor
		self.m1           = m1
		self.m2           = m2
		self.m3           = m3
		self.Pr           = 0.1973
		self.nParAll      = 2 # Number of parameters, this function uses
		self.fitPr        = fitPr
		self.barrier      = barrierFactors
		self.compensatePr = True # If True and self.fitPr, the resultung
					 # amplitude will be compensated for 
					 # self.pr. Meaning, that it will be 
					 # assumed, that in the normalization, 
					 # barrier factors with Pr = 0.1973 have
					 # been used, which will be divided out 
					 # again.
		if self.fitPr:
			self.nParAll += 2
		self.makeLoadMap(parameters)

	def __call__(self, ms, externalKinematicVariables = []):
		par = [p.value for p in self.parameters]

#		print par,r"|||\\\||"

		if len(externalKinematicVariables) < 1:
			raise IndexError("No external kinematic variable given")

		M0 = par[0].real
		G0 = par[1].real

		q0 = physUtils.breakupMomentum(M0, self.m1, self.m2)

		if q0 == 0.:
			return np.asarray([0.+0.j]*len(ms), dtype = complex)
		retVals = np.zeros((len(ms)), dtype = complex)
		if self.fitPr:
			PrMother = par[2]
			PrIsob = par[3]
		else:
			PrMother = self.Pr
			PrIsob = self.Pr

		for i,M in enumerate(ms):
			q = physUtils.breakupMomentum(M, self.m1, self.m2)
			retVals[i] = physUtils.breitWigner(M,M0,G0, self.J, q, q0, PrIsob)
			if self.fitPr and self.compensatePr:
				compensationFactor = physUtils.barrierFactor(self.J, q, PrIsob)/physUtils.barrierFactor(self.J, q, self.Pr)
				Q = physUtils.breakupMomentum(externalKinematicVariables[0],  M, self.m3)
				if Q == 0.:
					compensationFactor = 0.
				else:
					compensationFactor *= physUtils.barrierFactor(self.L, Q, PrMother)/physUtils.barrierFactor(self.L, Q, self.Pr)
				retVals[i] *= compensationFactor
		return retVals

class flatte(parametersTrackingParameterization):
	def __init__(self, parameters, m1, m2, m3, m2nd1,m2nd2, J, L, fitPr = False):
		self.L            = L * globalLfactor
		self.J            = J * globalLfactor
		self.m1           = m1
		self.m2           = m2
		self.m3           = m3
		self.m2nd1        = m2nd1 # Masses for the second channel
		self.m2nd2        = m2nd2 # Masses for the second channel
		self.Pr           = 0.1973
		self.nParAll      = 3 # Number of parameters, this function uses
		self.fitPr        = fitPr
		self.compensatePr = True # If True and self.fitPr, the resultung
					 # amplitude will be compensated for 
					 # self.pr. Meaning, that it will be 
					 # assumed, that in the normalization, 
					 # barrier factors with Pr = 0.1973 have
					 # been used, which will be divided out 
					 # again.
		if self.fitPr:
			self.nParAll += 2
		self.makeLoadMap(parameters)

	def __call__(self, ms, externalKinematicVariables = []):
		par = [p.value for p in self.parameters]

		if len(externalKinematicVariables) < 1:
			raise IndexError("No external kinematic variable given")

		M0 = par[0].real
		g1 = par[1].real
		g2 = par[2].real

		q0 = physUtils.breakupMomentum(M0, self.m1, self.m2)

		if q0 == 0.:
			return np.asarray([0.+0.j]*len(ms), dtype = complex)
		retVals = np.zeros((len(ms)), dtype = complex)
		if self.fitPr:
			PrMother = par[3]
			PrIsob = par[4]
		else:
			PrMother = self.Pr
			PrIsob = self.Pr

		for i,M in enumerate(ms):
			q = physUtils.breakupMomentum(M, self.m1, self.m2)
			retVals[i] = physUtils.flatte(M,M0,g1,g2,self.m1,self.m2,self.m2nd1,self.m2nd2)
			if self.fitPr and self.compensatePr:
				compensationFactor = physUtils.barrierFactor(self.J, q, PrIsob)/physUtils.barrierFactor(self.J, q, self.Pr)
				Q = physUtils.breakupMomentum(externalKinematicVariables[0],  M, self.m3)
				if Q == 0.:
					compensationFactor = 0.
				else:
					compensationFactor *= physUtils.barrierFactor(self.L, Q, PrMother)/physUtils.barrierFactor(self.L, Q, self.Pr)
				retVals[i] *= compensationFactor
		return retVals

class relativisticBreitWignerRatio(parametersTrackingParameterization):
	def __init__(self, parameters, m1, m2, m3, J, L, fitPr = False):
		self.L            = L * globalLfactor
		self.J            = J * globalLfactor
		self.m1           = m1
		self.m2           = m2
		self.m3           = m3
		self.Pr           = 0.1973
		self.nParAll      = 3 # Number of parameters, this function uses
		self.fitPr        = fitPr
		self.compensatePr = True # If True and self.fitPr, the resultung
					 # amplitude will be compensated for 
					 # self.pr. Meaning, that it will be 
					 # assumed, that in the normalization, 
					 # barrier factors with Pr = 0.1973 have
					 # been used, which will be divided out 
					 # again.
		if self.fitPr:
			self.nParAll += 2
		self.makeLoadMap(parameters)

	def __call__(self, ms, externalKinematicVariables = []):
		par = [p.value for p in self.parameters]

		if len(externalKinematicVariables) < 1:
			raise IndexError("No external kinematic variable given")

		M0  = par[0].real
		G0  = par[1].real
		rat = abs(par[2].real)**.5 # Scaling ratio of real/imag

		q0 = physUtils.breakupMomentum(M0, self.m1, self.m2)

		if q0 == 0.:
			return np.asarray([0.+0.j]*len(ms), dtype = complex)
		retVals = np.zeros((len(ms)), dtype = complex)
		if self.fitPr:
			PrMother = par[3]
			PrIsob   = par[4]
		else:
			PrMother = self.Pr
			PrIsob   = self.Pr

		for i,M in enumerate(ms):
			q = physUtils.breakupMomentum(M, self.m1, self.m2)
			retVals[i] = physUtils.breitWigner(M,M0,G0, self.J, q, q0, PrIsob)
			if self.fitPr and self.compensatePr:
				compensationFactor = physUtils.barrierFactor(self.J, q, PrIsob)/physUtils.barrierFactor(self.J, q, self.Pr)
				Q = physUtils.breakupMomentum(externalKinematicVariables[0],  M, self.m3)
				if Q == 0.:
					compensationFactor = 0.
				else:
					compensationFactor *= physUtils.barrierFactor(self.L, Q, PrMother)/physUtils.barrierFactor(self.L, Q, self.Pr)
				retVals[i] *= compensationFactor
			retVals[i] = retVals[i].real*rat + 1.j*retVals[i].imag/rat
		return retVals

class integratedRelativisticBreitWigner(parametersTrackingParameterization):
	def __init__(self, parameters, m1, m2, m3, J, L, binning, intensWeight = False, nPoints = 10, fitPr = False, reweightInverseBW = False):
		self.L            = L * globalLfactor
		self.J            = J * globalLfactor
		self.m1           = m1
		self.m2           = m2
		self.m3           = m3
		self.Pr           = 0.1973
		self.nParAll      = 2 # Number of parameters, this function uses
		self.fitPr        = fitPr
		self.compensatePr = True # If True and self.fitPr, the resultung
					 # amplitude will be compensated for 
					 # self.pr. Meaning, that it will be 
					 # assumed, that in the normalization, 
					 # barrier factors with Pr = 0.1973 have
					 # been used, which will be divided out 
					 # again.
		self.binning      = binning
		self.intensWeight = intensWeight
		self.rwibw        = reweightInverseBW
		self.nPoints      = nPoints
		if self.fitPr:
			self.nParAll += 2
		self.makeLoadMap(parameters)

	def __call__(self, ms, externalKinematicVariables = []):
		par = [p.value for p in self.parameters]

		if len(externalKinematicVariables) < 1:
			raise IndexError("No external kinematic variable given")

		M0 = par[0].real
		G0 = par[1].real

		q0 = physUtils.breakupMomentum(M0, self.m1, self.m2)

		if q0 == 0.:
			return np.asarray([0.+0.j]*len(ms), dtype = complex)
		retVals = np.zeros((len(ms)), dtype = complex)
		if self.fitPr:
			PrMother = par[2]
			PrIsob = par[3]
		else:
			PrMother = self.Pr
			PrIsob = self.Pr

		for i,Mcntr in enumerate(ms):
			found = False
			for i in range(len(self.binning)-1):
				if self.binning[i] < Mcntr and self.binning[i+1] >= Mcntr:
					mMin = self.binning[i  ]
					mMax = self.binning[i+1]
					found = True
					break
			if not found:
				raise ValueError("Mass " + str(Mcntr) + " invalid")
			retVal = 0.
			weight = 0.
			step   = (mMax - mMin)/self.nPoints
			for b in range(self.nPoints):
				M = mMin + (b + 0.5)*step
				q = physUtils.breakupMomentum(M, self.m1, self.m2)
				ampl = physUtils.breitWigner(M,M0,G0, self.J, q, q0, PrIsob)
				if self.fitPr and self.compensatePr:
					compensationFactor = physUtils.barrierFactor(self.J, q, PrIsob)/physUtils.barrierFactor(self.J, q, self.Pr)
					Q = physUtils.breakupMomentum(externalKinematicVariables[0],  M, self.m3)
					if Q == 0.:
						compensationFactor = 0.
					else:
						compensationFactor *= physUtils.barrierFactor(self.L, Q, PrMother)/physUtils.barrierFactor(self.L, Q, self.Pr)
					ampl *= compensationFactor
				ww = 1.
				if self.intensWeight:
					ww *= abs(ampl)**2
				if self.rwibw:
					Q    = physUtils.breakupMomentum(externalKinematicVariables[0],  M, self.m3)
					fakk = physUtils.barrierFactor(self.J, q, PrIsob) * physUtils.barrierFactor(self.L, Q, PrMother)
					ww  *= fakk
				retVal += ww*ampl
				weight += ww
			retVals[i] = retVal/weight
		return retVals

class timesBF:
	def __init__(self, BW):
		self.BW   = BW
		self.m1   = BW.m1
		self.m2   = BW.m2
		self.m3   = BW.m3
		self.L    = BW.L
		self.J    = BW.J
		self.Pr   = BW.Pr
		self.nPar = BW.nPar

	def __call__(self, ms, externalKinematicVariables = []):
		retVals = self.BW(ms, externalKinematicVariables)
		for i,m in enumerate(ms):
			q   = physUtils.breakupMomentum(m,self.m1, self.m2)
			Q   = physUtils.breakupMomentum(externalKinematicVariables[0], m, self.m3)
			BFs = physUtils.barrierFactor(self.L, Q, self.Pr)*physUtils.barrierFactor(self.J, q, self.Pr)
			retVals[i] *= BFs
		return retVals

	def setParameters(self, params):
		self.BW.setParameters(params)

	def setParametersAndErrors(self, params, errors):
		return self.BW.setParametersAndErrors(params, errors)

	def getParameters(self):
		return self.BW.getParameters()
	
	def returnParameters(self):
		return self.BW.returnParameters()

def main():
	commonMass = parameter(1., 'commonMass')
	width1     = parameter(.1, 'width1'    )
	width2     = parameter(.2, 'width2'    )

	BW1 = 	relativisticBreitWigner([commonMass, width1], 0.1, 0.1, 0.1, 1, 1, False)
	BW2 = 	relativisticBreitWigner([commonMass, width2], 0.1, 0.1, 0.1, 1, 1, False)

	print BW1.nPar
	print BW2.nPar

if __name__ == "__main__":
	main()
			
