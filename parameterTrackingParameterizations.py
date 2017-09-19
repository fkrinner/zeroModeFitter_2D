import numpy as np
import physUtils

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
			raise IndexError("Number of parameters does not match")
		count = 0
		for i,v in enumerate(self.loadMap):
			if v:
				self.parameters[i].value = params[count]
				count += 1

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
		if not len(parameters) == self.nParAll:
			raise IndexError("Number of parameters does not match")
		self.loadMap = []
		self.nPar    = 0
		self.parameters = parameters
		for p in self.parameters:
			if not p.lock:
				self.loadMap.append(True)
				self.nPar += 1
				p.lock = True
			else:
				self.loadMap.append(False)

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
		self.nPar       = 0
	
	def __call__(self, ms, externalKinematicVariables = []):
		retVals = np.zeros((len(ms)), dtype = complex)
		for i in range(len(retVals)):
			retVals[i] = 1.
		return retVals


class breitWigner(parametersTrackingParameterization):
	def __init__(self, parameters):
		self.nPar       = 0
		self.nParAll    = 2
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
	def __call__(self, ms, externalKinematicVariables = []):
		par = [p.value for p in self.parameters]
		retVals = np.zeros((len(ms)), dtype = complex)
		num = par[0]*par[1]
		den = par[0]**2 - 1.j * par[0] * par[1]
		for i, m in enumerate(ms):
			retVals[i] = num/(den - m**2)
		return retVals

globalLfactor = 1 # probobly not needed anymore

class relativisticBreitWigner(parametersTrackingParameterization):
	def __init__(self, parameters, m1, m2, m3, J, L, fitPr = False):
		self.L            = L * globalLfactor
		self.J            = J * globalLfactor
		self.m1           = m1
		self.m2           = m2
		self.m3           = m3
		self.Pr           = 0.1973
		self.nPar         = 0 # Number of parameters, this function handels
		self.nParAll      = 2 # Number of parameters, this function uses
		self.loadMap      = []
		self.fitPr        = fitPr
		self.parameters   = parameters
		self.compensatePr = True # If True and self.fitPr, the resultung
					 # amplitude will be compensated for 
					 # self.pr. Meaning, that it will be 
					 # assumed, that in the normalization, 
					 # barrier factors with Pr = 0.1973 have
					 # been used, which will be divided out 
					 # again.
		if self.fitPr:
			self.nParAll += 2
		if not len(parameters) == self.nParAll:
			raise IndexError("Number of parameters does not match (" + str(len(parameters))+"/"+str(self.nParAll)+")")
		for p in parameters:
			if not p.lock:
				self.loadMap.append(True)
				self.nPar += 1
				p.lock = True
			else:
				self.loadMap.append(False)

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
		self.nPar         = 0 # Number of parameters, this function handels
		self.nParAll      = 3 # Number of parameters, this function uses
		self.loadMap      = []
		self.fitPr        = fitPr
		self.parameters   = parameters
		self.compensatePr = True # If True and self.fitPr, the resultung
					 # amplitude will be compensated for 
					 # self.pr. Meaning, that it will be 
					 # assumed, that in the normalization, 
					 # barrier factors with Pr = 0.1973 have
					 # been used, which will be divided out 
					 # again.
		if self.fitPr:
			self.nParAll += 2
		if not len(parameters) == self.nParAll:
			raise IndexError("Number of parameters does not match (" + str(len(parameters))+"/"+str(self.nParAll)+")")
		for p in parameters:
			if not p.lock:
				self.loadMap.append(True)
				self.nPar += 1
				p.lock = True
			else:
				self.loadMap.append(False)

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
		self.nPar         = 0          # Number of parameters, this function handels
		self.nParAll      = 3 # Number of parameters, this function uses
		self.loadMap      = []
		self.fitPr        = fitPr
		self.parameters   = parameters
		self.compensatePr = True # If True and self.fitPr, the resultung
					 # amplitude will be compensated for 
					 # self.pr. Meaning, that it will be 
					 # assumed, that in the normalization, 
					 # barrier factors with Pr = 0.1973 have
					 # been used, which will be divided out 
					 # again.
		if self.fitPr:
			self.nParAll += 2
		if not len(parameters) == self.nParAll:
			raise IndexError("Number of parameters does not match (" + str(len(parameters))+"/"+str(self.nParAll)+")")
		for p in parameters:
			if not p.lock:
				self.loadMap.append(True)
				self.nPar += 1
				p.lock = True
			else:
				self.loadMap.append(False)

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
		self.nPar         = 0 # Number of parameters, this function handels
		self.nParAll      = 2 # Number of parameters, this function uses
		self.loadMap      = []
		self.fitPr        = fitPr
		self.parameters   = parameters
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
		if not len(parameters) == self.nParAll:
			raise IndexError("Number of parameters does not match (" + str(len(parameters))+"/"+str(self.nParAll)+")")
		for p in parameters:
			if not p.lock:
				self.loadMap.append(True)
				self.nPar += 1
				p.lock = True
			else:
				self.loadMap.append(False)

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
			
