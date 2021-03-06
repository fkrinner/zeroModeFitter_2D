from modes import INTENS, PHASE
import numpy as np
class tBinHolder:
	"""
	Class to hold an globally manipulate several instances of 'allBins' corresponting to t' bins
	"""
	def __init__(self):
		"""
		Simplemost initializer
		"""
		self.bestChi2UpToNow = float('inf')
		self.bestPar         = None
		self.bins            = []
		
	def addBin(self, bin):
		"""
		Simplemost append method
		"""
		self.bins.append(bin)

	def chi2(self, params):
		"""
		Evaluation of the chi2 function, summed over t'
		"""
		chi2 = 0
#		print "calllolo"
		for bin in self.bins:
			chi2 += bin.chi2(params)
#		print "chi2 call with", params, 'giving', chi2
#		print ">>>",chi2,"<<<"
		if chi2 < self.bestChi2UpToNow:
			self.bestChi2UpToNow = chi2
			self.bestPar = params[:]
		return chi2

	def addComaValueForZeroMode(self, val, unitsOf = 'smallestComaValue'):
		"""
		Simple t' bin loop of the 'addComaValueForZeroMode(...)' method
		"""
		for tb in self.bins:
			tb.addComaValueForZeroMode(val, unitsOf = unitsOf)

	def getTheoryTotalMatrices(self,binRange = None):
		"""
		Simple t' bin loop
		"""
		retVals = []
		for tb in self.bins:
			retVals.append(tb.getTheoryTotalMatrices(binRange = binRange))
		return retVals

	def compareTwoZeroModeCorrections(self, params1, params2):
		"""
		Compares two different corrections of the zero-modes
		"""
		if not len(params1) == len(self.bins):
			raise ValueError("Number of t' bins does not match (1)" )
		if not len(params2) == len(self.bins):
			raise ValueError("Number of t' bins does not match (2)" )
		retVal = []
		print "----in tbin holder"
		print '----',params1
		print '----',params2
		print "----out tbin holder"
		for t, tb in enumerate(self.bins):
			retVal.append(tb.compareTwoZeroModeCorrections(params1[t], params2[t]))
		return retVal

	def nParAll(self):
		"""
		Gives the number of total parameters in the model, assuming an equal number of parameters in all t' 
		"""
		nPar = 0
		return len(self.bins) * (self.bins[0].nParAll() - self.bins[0].nPar()) + self.bins[0].nPar()

	def nPar(self):
		"""
		Gives the number of shape-parameters in the model
		"""
		return self.bins[0].nPar()

	def getNDF(self):
		"""
		Return NDF as list of (NDF, 2*nZero, 2*nFunc, nPar) (Simple t' bin loop)
		"""
		retVal = []
		for tb in self.bins:
			retVal.append(tb.getNDF())
		return retVal

	def removeAllCorrelations(self):
		"""
		Removes ALL correlations from the covariance matrix (Simple t' bin loop)
		"""
		for tb in self.bins:
			tb.removeAllCorrelations()

	def modeChi2(self, pars, mode = PHASE, tBinResolved = False, mBinResolved = False):
		"""
		Evaluates a chi2 for a specific mode, other than the amplitude (PHASE or INTENS), summed over all t' bins
		"""
		if not len(pars) == self.nParAll():
			print len(pars), self.nParAll()
			raise IndexError("Number of parameters does not match")
		if mBinResolved and not tBinResolved:
			raise ValueError("tBinResolved must be True, if mBinResolved is true")

		nPar  = self.bins[0].nPar()
		nNon  = self.bins[0].nParAll() - nPar
		count = 0
		if tBinResolved:
			chi2 = []
		else:
			chi2  = 0.
		for tb in self.bins:
			par   = pars[count:count+nNon].tolist()
			if nPar > 0:
				par += pars[-nPar:].tolist()
			c2    = tb.modeChi2(par, mode = mode, mBinResolved = mBinResolved)
			if tBinResolved:
				chi2.append(c2)
			else:
				chi2 += c2
			count += nNon
		return chi2

	def getFixedZMndf(self):
		"""
		Returns the NDF for the chi2 function of fixedZMPchi2
		"""
		if len(self.binsToEvaluate) == 0:
			raise RuntimeError("No bins to evaluate set")
		NDF    = 0
		nShape = None
		for i in self.binsToEvaluate:
			ndf, nSh = self.bins[i].getFixedZMndf()
			NDF += ndf
			if not nShape:
				nShape = nSh
			elif not nShape == nSh:
				raise IndexError("nShape mismatch in t' sum NDF method")
		return NDF - nSh

	def fixedZMPchi2(self, pars):
		"""
		Returns a chi2 for the shape parameters and self.zeroModeParameters. The couplings are calculated. Sums over all bins in self.binsToEvaluate
		"""
		if len(self.binsToEvaluate) == 0:
			raise RuntimeError("No bins to evaluate set")
		chi2 = 0.
		for i in self.binsToEvaluate:
			chi2 += self.bins[i].fixedZMPchi2(pars)
		return chi2

	def fixedZMPchi2_realCouplings(self, pars, phases):
		"""
		Returns a chi2 for the shape parameters and self.zeroModeParameters. The couplings are calculated. Sums over all bins in self.binsToEvaluate
		"""
		if len(self.binsToEvaluate) == 0:
			raise RuntimeError("No bins to evaluate set")
		chi2 = 0.
		phasesPerTbin = len(phases)/len(self.binsToEvaluate)

		for i,t in enumerate(self.binsToEvaluate):
			chi2 += self.bins[t].fixedZMPchi2_realCouplings(pars, phases[i*phasesPerTbin:(i+1)*phasesPerTbin])
		return chi2

	def setBinsToEvalueate(self, tBinsToEvaluate,mBinsToEvaluate):
		"""
		Sets the t' and m bins to evaluate in "fitRho.fitShapeParametersForBinRange()" in the analysisclass
		"""
		self.binsToEvaluate = tBinsToEvaluate
		for bin in self.bins:
			bin.setBinsToEvalueate(mBinsToEvaluate)

	def setZeroModeParameters(self, zmp):
		"""
		Sets zero mode parameters
		"""
		if not len(zmp) == len(self.bins):
			raise IndexError("Mismatch in number of t' bins")
		for i,pp in enumerate(zmp):
			self.bins[i].setZeroModeParameters(pp)
	
	def __getitem__(self, index):
		"""
		Enables looping through this class
		"""
		if index < 0 or index >= len(self.bins):
			raise IndexError("Index out of range: "+str(0) +" <= i < " + str(len(self.bins))) # This works with for ... in ...
		return self.bins[index]
	
	def __len__(self):
		"""
		Enables the use of len()
		"""
		return len(self.bins)
