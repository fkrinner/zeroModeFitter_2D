import numpy as np

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

	def __call__(self, ms, par = None):
		if not par:
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

	def __call__(self, ms, par = None):
		if not par:
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

	def __call__(self, ms, par = None):
		if not par:
			par = self.parameters
		if not len(par) == self.nPar:
			raise IndexError("omnesFunctionMonomial: Wrong number of parameters")
		retVals = np.zeros((len(ms)), dtype = complex)
		for i,m in enumerate(ms):
			retVals[i] =  m**(2*self.degree) * lookupOmnes(m**2)
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

	def __call__(self, ms, par = None):
		if not par:
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
