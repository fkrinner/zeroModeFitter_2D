import opp_data_zeroModeFixingComparison as opp
import os, sys
from studyFileNames import fileNameMap

from vectorFormFactor2000 import vectorFormFactor
import parameterTrackingParameterizations as ptc

from globalDefinitions import referenceWave, mPi, mK,mRho,Grho,mRhoPrime,GrhoPrime,mF0,g1,g2,mF2,GF2,Pr0,m1500,G1500
from utils import checkLaTeX

from iminuit import Minuit
from globalDefinitions import referenceWave

def explicit_MINUIT_function():
	pass

def function():
	pass

def MIGRAD(argString, lstString, function2):
	global function
	function = function2
	command  = "def eMf("+lstString+"):\n\treturn function(["+lstString+"])"
	names    = lstString.split(',')
	global explicit_MINUIT_function
	exec command
	explicit_MINUIT_function = eMf
	command ="m = Minuit(explicit_MINUIT_function,"+argString+")"
	exec command
	print "Starting migrad"
	m.migrad()
	vals = [m.values[name.strip()] for name in names]
	return vals

def minuit_function(*args):
	return explicit_function(args)

def migrad(function, parameters):
	names = [p.name for p in parameters]
	global explicit_function
	explicit_function = function
	kwargs = {}
	for p in parameters:
		kwargs[p.name] = p.value
	
	itovaller = [p.value for p in parameters]
	print "Chi2 at start values:", function(itovaller)

	m = Minuit(minuit_function, forced_parameters = names, pedantic = False, **kwargs)
	fmin, param = m.migrad()
	vals = [m.values[name.strip()] for name in names]
	return vals

def main():
	checkLaTeX()
	inFileName = fileNameMap['std11']

	mMin  = None
	mMax  = None
	tBins = []
	for a in sys.argv:
		if a.startswith("mMin"):
			if mMin is not None:
				raise RuntimeError("mMin set twice")
			else:
				mMin = int(a[4:])

		if a.startswith("mMax"):
			if mMax is not None:
				raise RuntimeError("mMax set twice")
			else:
				mMax = int(a[4:])

		if a.startswith('t'):
			t = int(a[1:])
			if not t in range(4):
				raise RuntimeError("invelide tBin: " + a)
			if t in tBins:
				raise RuntimeError("tBin speciofitd twice: " + a)
			tBins.append(t)
	if mMin is None:
		raise RuntimeError("mMin not set")
	if mMax is None:
		raise RuntimeError("mMax not set")
	if not mMin < mMax:
		raise RuntimeError("m borders not ordered")
	if len(tBins) == 0:
		raise RuntimeError("No t' bin set")		

	tBins.sort()
	acv = None
	sectorRangeMap = {}
	lockSecondRho  = True
	
	parameters = [	ptc.parameter( .36674e-01, "lambP" ),
			ptc.parameter( .31230e-02, "lambPP"),
			ptc.parameter( .83386    , "M"     ),
			ptc.parameter( .19771    , "G"     ),
			ptc.parameter(1.4974     , "MP"    ),
			ptc.parameter( .78518    , "GP"    ),
			ptc.parameter(1.6855     , "MPP"   ),
			ptc.parameter( .80109    , "GPP"   )]

	lock = []
	if lockSecondRho:
		# fix m and G of the second rho, due to (1)
		lock.append(4)
		lock.append(5)

	fitters  = []
	binCount = 0
	fitters  = []
	for t in tBins:
		for m in range(mMin, mMax):
			parameters.append(ptc.parameter( .17254, "alP_"+str(t)+"_"+str(m)))
			parameters.append(ptc.parameter(-.97591, "phP_"+str(t)+"_"+str(m)))
			parameters.append(ptc.parameter( .23374, "alPP_"+str(t)+"_"+str(m)))
			parameters.append(ptc.parameter(2.1982,  "phPP_"+str(t)+"_"+str(m)))

			if lockSecondRho:
				# remove the second rho, since it just might modulate the intermediate region (1)
				lock.append(len(parameters)-2)
				lock.append(len(parameters)-3)

				parameters[-2].value = 0.
				parameters[-3].value = 0.

	for p in parameters:
		p.lock = True

	allowSubThr = True
	NDF = 0
	for t in tBins:
		for m in range(mMin, mMax):
			vff = vectorFormFactor(parameters[:8] + parameters[8+binCount*4:8+(binCount+1)*4])
			vff.allowSubThr = allowSubThr
			fitter = opp.produceSingleMassTbinFitter(inFileName, [vff], m, t, sectorRangeMap = sectorRangeMap, referenceWave = referenceWave, acv = acv, removeZM = False)
			NDF += fitter.getNDFforMode()
			fitters.append(fitter)
			binCount += 1

	def chi2(par):
		if not len(par) == len(parameters) - len(lock):
			raise ValueError("parameter size mismatch")
		count = 0
		for i,p in enumerate(parameters):
			if not i in lock:
				p.value = par[count]
				count  += 1
		c2 = 0.
		for f in fitters:
			c2 += f.model[0][0].chi2()
		return c2

	startParameters = []
	for i,p in enumerate(parameters):
		if not i in lock:
			startParameters.append(p)
	NDF -= len(startParameters)
	vals = migrad(chi2, startParameters)

	outFileName = "./global_vff_fits_lockRhoPrime/m_"+str(mMin)+'-'+str(mMax)+"_t"
	for t in tBins:
		outFileName += str(t)
	outFileName += ".dat"

	with open(outFileName, 'w') as outFile:
		for i,p in enumerate(startParameters):
			outFile.write(p.name + ' ' + str(vals[i]) + '\n')
		outFile.write(str(chi2(vals))+'/'+str(NDF))

	doPlots = True
	if doPlots:
		fitterCount = 0
		for t in tBins:
			for m in range(mMin, mMax):
				fitters[fitterCount].model.chi2([])
				fitters[fitterCount].SET('hasFitResult')
				fitters[fitterCount].fitParameters = []

				zeroModeParameters = fitters[fitterCount].calculateNonShapeParameters()
				zeroModeParameters = fitters[fitterCount].getZeroModeParametersForMode()
				RV = fitters[fitterCount].produceResultViewer(zeroModeParameters,"1++0+[pi,pi]1--PiS", noRun = True, plotTheory = True)
				RV.plotData = True
				RV.connectTheoPoints = range(100)
				plotNameBase = "./vectorFormFactor_lockRhoPrime_plots/1pp0p1mmPiS_<mode>_"+str(m)+"_"+str(t)+".pdf"
				RV.writeBinToPdf(m, stdCmd = ["", plotNameBase.replace("<mode>","intens"), [],  plotNameBase.replace("<mode>","argand"), []])
				fitterCount += 1
	return

if __name__ == "__main__":


	main()
