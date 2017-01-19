import numpy as np
import pyRootPwa

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

class breitWigner:
	def __init__(self):
		self.nPar       = 2
		self.parNames   = ["mass", "width"]
		self.parameters = [0.] * self.nPar

	def __call__(self, ms, par = None, externalKinematicVariables = []):
		if par is None:
			par = self.parameters
		if not len(par) == self.nPar:
			raise IndexError("breitWigner: Wrong number of parameters")
		retVals = np.zeros((len(ms)), dtype = complex)
		num = par[0]*par[1]
		den = par[0]**2 - 1.j * par[0] * par[1]
		for i, m in enumerate(ms):
			retVals[i] = num/(den - m**2)
		return retVals

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

class omnesFunctionPolynomial:
	"""
	Omnes polynomial allows to set parameters to real. For complex coefficients it it more efficient to use monomials
	"""
	def __init__(self, nDimPol = 0, complexPolynomial = False):
		# Since a complex factor will be mutliplied in the fit, set c0 = 1., then:
		# c0 + c1 m + c2 m**2 + ... + cn m**n -> 1. + c1 m + ... + cn m**n
		self.complexPolynomial = complexPolynomial
		self.polDeg = nDimPol
		if complexPolynomial:
			self.nPar = 2 * nDimPol
		else:
			self.nPar = nDimPol
		self.parNames = []
		for i in range(nDimPol):
			if complexPolynomial:
				self.parNames.append("c"+str(i+1)+"_real")
				self.parNames.append("c"+str(i+1)+"_imag")
			else:
				self.parNames.append("c"+str(i+1))
		self.parameters = [0.]*self.nPar

	def __call__(self, ms, par = None, externalKinematicVariables = []):
		if par is None:
			par = self.parameters
		if not len(par) == self.nPar:
			raise IndexError("omnesFunctionPolynomial: Wrong number of parameters")
		retVals = np.zeros((len(ms)), dtype = complex)
		for i,m in enumerate(ms):
			s = m**2
			x = s
			poly = 1.
			for d in range(self.polDeg):
				if self.complexPolynomial:
					poly += x * (par[2*d] + 1.j*par[2*d+1])
				else:
					poly += x * par[d]
				x *= s
			retVals[i] = poly * lookupOmnes(s)
		return retVals

class omnesFunctionMonomial:
	"""
	Several omnes monomiala can use the intrinsic linearity in the parameters and avoid fitting (Real coefficients, whoever are not possible)
	"""
	def __init__(self, degree):
		self.degree     = degree
		self.nPar       =  0
		self.parameters = [ ]
		self.parNames   = [ ]

	def __call__(self, ms, par = None, externalKinematicVariables = []):
		if par is None:
			par = self.parameters
		if not len(par) == self.nPar:
			raise IndexError("omnesFunctionMonomial: Wrong number of parameters")
		retVals = np.zeros((len(ms)), dtype = complex)
		for i,m in enumerate(ms):
			retVals[i] =  m**(2*self.degree) * lookupOmnes(m**2)
		return retVals

class rpwaBreitWigner:
	"""
	Breit-Wigner as defined in rootPwa
	"""
	def __init__(self, m1, m2, m3, J, L):
		self.L          = 2*L
		self.J          = 2*J
		self.m1         = m1
		self.m2         = m2
		self.m3         = m3
		self.nPar       = 2
		self.parameters = [1.,.1]
		self.parNames   = ["mass", "width"]

	def __call__(self, ms, par = None, externalKinematicVariables = []):
		if par is None:
			par = self.parameters
		if not len(par) == self.nPar:
			raise IndexError("rpwaBreitWigner: Wring number of parameters")
		if len(externalKinematicVariables) < 1:
			raise IndexError("No external kinematic variable given")

		M0 = par[0].real
		G0 = par[1].real

		q0 = abs(pyRootPwa.core.breakupMomentumSquared(M0, self.m1, self.m2, True))**.5
		if q0 == 0.:
			return [0.]*len(ms)
		retVals = np.zeros((len(ms)), dtype = complex)
		for i,M in enumerate(ms):
			q = pyRootPwa.core.breakupMomentum(M,  self.m1, self.m2)
			if M + self.m3 >= externalKinematicVariables[0]:
				Q = 0.
			else:
				Q = pyRootPwa.core.breakupMomentum(externalKinematicVariables[0],  M, self.m3)
			retVals[i] = pyRootPwa.core.breitWigner(M, M0, G0, self.J, q, q0) # * pyRootPwa.core.barrierFactorSquared(self.J, q)**.5*pyRootPwa.core.barrierFactorSquared(self.L, Q)**.5
		return retVals

class rpwaBreitWignerInt:
	"""
	Breit-Wigner as defined in rootPwa integrated over the bin
	"""
	def __init__(self, m1, m2, m3, J, L, nPoints = 10):
		self.L            = 2*L
		self.J            = 2*J
		self.m1           = m1
		self.m2           = m2
		self.m3           = m3
		self.nPar         = 2
		self.nPoints      = nPoints
		self.parameters   = [1.,.1]
		self.parNames     = ["mass", "width"]
		self.weightBF     = False
		self.weightIntens = False

	def __call__(self, ms, par = None, externalKinematicVariables = []):
		if par is None:
			par = self.parameters
		if not len(par) == self.nPar:
			raise IndexError("rpwaBreitWigner: Wrong number of parameters")
		if len(externalKinematicVariables) < 1:
			raise IndexError("No external kinematic variable given")

		M0 = par[0].real
		G0 = par[1].real
		q0 = abs(pyRootPwa.core.breakupMomentumSquared(M0, self.m1, self.m2, True))**.5
		if q0 == 0.:
			return [0.]*len(ms)
		retVals = np.zeros((len(ms)), dtype = complex)
		for i,MM in enumerate(ms):

			if i > 0: 
				binWidth = 2*(MM - ms[i-1]) - binWidth
			else:
				binWidth = 0.042 # First bin width is different... use fixed known value here HARDCODE!!! <BAD> <<!! HARDCORE !!>>
			ampl = 0.+0.j
			step   = binWidth/(self.nPoints)
			totalWeight = 0.
			for p in range(self.nPoints):
				M     = MM - binWidth/2 + step*(p + .5)
				if M + self.m3 >= externalKinematicVariables[0]:
					Q = 0.
				else:
					Q = pyRootPwa.core.breakupMomentum(externalKinematicVariables[0],  M, self.m3)

				if M <= self.m1 + self.m2:
					continue
				q = pyRootPwa.core.breakupMomentum(M,  self.m1, self.m2)
				if self.weightBF:
					try:
						weight = (pyRootPwa.core.barrierFactorSquared(self.J, q)**.5*pyRootPwa.core.barrierFactorSquared(self.L, Q)**.5)
					except ZeroDivisionError:
						weight = 0.
				elif self.weightIntens:
					weight = abs(pyRootPwa.core.breitWigner(M, M0, G0, self.J, q, q0))**2
				else:
					weight = 1.

				ampl += pyRootPwa.core.breitWigner(M, M0, G0, self.J, q, q0) * weight
				totalWeight += weight
			retVals[i] = ampl/totalWeight
		return retVals

class m3PiHalf:
	"""
	Function to test the external kinematic variables, setting everything below m3Pi/2 to 1., above to 0.
	"""
	def __init__(self):
		self.nPar       = 0
		self.parameters = []
		self.parNames   = []

	def __call__(self, ms, par = None, externalKinematicVariables = []):
		if par is None:
			par = self.parameters
		if not len(par) == self.nPar:
			raise IndexError("rpwaBreitWigner: Wrong number of parameters")
		if len(externalKinematicVariables) < 1:
			raise IndexError("No external kinematic variable given")
		retVals = np.zeros((len(ms)), dtype = complex)
		for i,m in enumerate(ms):
			if m < externalKinematicVariables[0]/2.:
				retVals[i] = 1.
		return retVals

		
class fixedParameterization:
	"""
	Fixed parameterization parsed from a text file, with mass real imag in each line
	"""
	def __init__(self, ampFileName, polynomialDegree = 0, complexPolynomial = False, baseExponent = 1):
		self.ampFileName       = ampFileName
		self.polynomialDegree  = polynomialDegree
		self.complexPolynomial = complexPolynomial
		self.baseExponent      = baseExponent
		if self.complexPolynomial:
			self.nPar = 2*polynomialDegree
		else:
			self.nPar = polynomialDegree
		self.parNames = []
		for i in range(self.polynomialDegree):
			if complexPolynomial:
				self.parNames.append("c"+str(i+1)+"_real")
				self.parNames.append("c"+str(i+1)+"_imag")
			else:
				self.parNames.append("c"+str(i+1))
		self.parameters = [0.]*self.nPar
		self.m,self.r,self.i = loadAmplitudeFile(self.ampFileName)

	def lookupValue(self, m):
		found = False	
		for i in range(len(self.m)-1):
			if m > self.m[i] and m <= self.m[i+1]:
				x = (m - self.m[i])/(self.m[i+1] - self.m[i])
				found = True
				break
		if not found:
			return 0.+0.j
		return x * (self.r[i+1] + 1.j*self.i[i+1]) + (1-x)*(self.r[i] + 1.j * self.i[i])

	def __call__(self, ms, par = None, externalKinematicVariables = []):
		if par is None:
			par = self.parameters
		if not len(par) == self.nPar:
			raise IndexError("fixedParameterization: Wrong number of parameters")
		retVals = np.zeros((len(ms)), dtype = complex)
		for i, m in enumerate(ms):
			polynom = 1.
			for j in range(self.polynomialDegree):
				monom = m**((j+1)*self.baseExponent)
				if self.complexPolynomial:
					monom *= par[2*j] + 1.j * par[2*j+1]
				else:
					monom *= par[j]
				polynom += monom
			retVals[i] =  polynom * self.lookupValue(m)
		return retVals

def main():
	ms = []
	of = omnesFunction()
	for i in range(2000):
		m = 0.001 * i
		ms.append(m)
	amps = of(ms)
	with open("omnes", 'w') as out:
		for i in range(len(ms)):
			out.write(str(ms[i]) + ' ' + str(abs(amps[i])**2) + ' ' + str(amps[i].real) + ' ' + str(amps[i].imag) + '\n')

if __name__ == "__main__":
	main()
