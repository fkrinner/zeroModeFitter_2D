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


class relativisticBreitWigner:
	def __init__(self, parameters, m1, m2, m3, J, L, fitPr = False):
		self.L            = L
		self.J            = J
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
			
