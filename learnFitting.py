from waveNameClass import waveName
from rootfabi import root_open, GetKeyNames
import numpy as np
import numpy.linalg as la
from allBinsClass import allBins
from massBinClass import massBin
from modes import REAL, IMAG, PHASE, INTENS, INTENSNORM, INTENSTHEO, REALTHEO, IMAGTHEO, PHASETHEO
from resultViewerClass import resultViewer
import os, sys
from utils import zeroForSectors, getZeroHistBorders, renormToBinWidth, LtoInt
import parameterizationClasses as pc
import parameterTrackingParameterizations as ptc
import scipy.optimize
from removeZeroModes import removeCertainZeroModes
from fixedparameterizationPaths import getFileNameForSector
from tBinHolder import tBinHolder

def main(rhoFileName = ""):
	inFileName = "/nfs/mds/user/fkrinner/extensiveFreedIsobarStudies/results_4pp.root"
#	sectors = ["4++1+[pi,pi]1--PiG", "4++1+[pi,pi]2++PiF"]
#	sectors = ["4++1+[pi,pi]1--PiG"]
#	sectors = ["4++1+[pi,pi]2++PiF"]
	sectors = ["1++0+[pi,pi]1--PiS"]
#	sectors = ["2++1+[pi,pi]1--PiD"]
#	sectors = ["2-+0+[pi,pi]2++PiS"]
#	sectors = ["2-+0+[pi,pi]2++PiS","4++1+[pi,pi]2++PiF" ]
#	sectors = ["0-+0+[pi,pi]1--PiP"]
#	sectors = ["1++0+[pi,pi]1--PiS", "2++1+[pi,pi]1--PiD", "4++1+[pi,pi]1--PiG", "4++1+[pi,pi]2++PiF"]
	sectorRangeMap = {}
	for sector in sectors:
		if '1--' in sector:
			pass
#			sectorRangeMap[sector] = (0.55,1.00)
		if '2++' in sector:
			pass

	tBins            = [0]#,1,2,3]
	startBin         = 20
	stopBin          = 50

	fitShape         = True
	fitPr            = True

	### Masses and widths

	mPi      = 0.13957018

	rhoMass  = ptc.parameter( .77549, "rhoMass" )
	rhoWidth = ptc.parameter( .1491 , "rhoWidth")

	f2Mass   = ptc.parameter( 1.2   , "f2Mass"  )
	f2Width  = ptc.parameter(  .2   , "f2Width" )

	### Interatction radii

	Pr0   = 0.1973

	PrA1  = ptc.PRparameter( Pr0, "PrA1" )
	PrA2  = ptc.PRparameter( Pr0, "PrA2" )
	PrPi2 = ptc.PRparameter( Pr0, "PrPi2")
	PrA4  = ptc.PRparameter( Pr0, "PrA4" )

	PrRho = ptc.PRparameter( Pr0, "PrRho")
	PrF2  = ptc.PRparameter( Pr0, "PrF2" )

	omnesParameters = []

	shift   = True
	stretch = True

	compl            = True
	polyNomialDegree = 1	

	omnesShift       = ptc.parameter( 0., "shift"  )
	omnesStretch     = ptc.parameter( 1., "stretch")



	if shift:
		omnesParameters.append(omnesShift)
	if stretch:
		omnesParameters.append(omnesStretch)
	for i in range(polyNomialDegree):
		if compl:
			omnesParameters.append(ptc.parameter(0., "cRe_"+str(i+1)))
			omnesParameters.append(ptc.parameter(0., "cIm_"+str(i+1)))
		else:
			omnesParameters.append(ptc.parameter(0., "c_"+str(i+1)))

	omnes = ptc.omnesFunctionPolynomial(omnesParameters, nDimPol = polyNomialDegree, shift = shift, stretch = stretch, complexPolynomial = compl)

	waveModel = {}
	PrMap     = { "4++1+" : PrA4, "2++1+" : PrA2, "1++0+": PrA1, "2-+0+": PrPi2 }

	allParameters = omnesParameters

	for sector in sectors:
#		L     = LtoInt(sector[-1])
#		JPC   = sector[:5]
#		Pr3Pi = PrMap[JPC]
#		if ']1--' in sector:
#			rhoParameters = [rhoMass, rhoWidth]
#			if fitPr:
#				rhoParameters.append(Pr3Pi)
#				rhoParameters.append(PrRho)
#			rho = ptc.relativisticBreitWigner(rhoParameters, mPi, mPi, mPi, 1, L, fitPr)
#			allParameters += rhoParameters
#			waveModel[sector] = [rho]
#		elif "]2++" in sector:
#			f2Parameters = [f2Mass, f2Width]
#			if fitPr:
#				f2Parameters.append(Pr3Pi)
#				f2Parameters.append(PrF2)
#			f2 = ptc.relativisticBreitWigner(f2Parameters, mPi, mPi, mPi, 2, L, fitPr)
#			allParameters += f2Parameters
#			waveModel[sector] = [f2]
		waveModel[sector] = [omnes]

	totalModel = tBinHolder()

	with root_open(inFileName, "READ") as inFile:
		for tBin in tBins:
			histNames = GetKeyNames(inFile)
			histListReal = []
			histListImag = []
			histListNorm = []
			histListIndx = []
			for sector in sectors:
				realHistName = sector+"_"+str(tBin)+"_real"
				histReal = inFile.Get(realHistName)
				if not histReal:
					raise IOError("Could not get '" + realHistName + "' from '" + inFileName + "'")
				histListReal.append(histReal)

				imagHistName = sector+"_"+str(tBin)+"_imag"
				histImag = inFile.Get(imagHistName)
				if not histImag:
					raise IOError("Could not get '" + imagHistName + "' from '" + inFileName + "'")
				histListImag.append(histImag)

				normHistName = sector+"_"+str(tBin)+"_norm"
				histNorm = inFile.Get(normHistName)
				if not histNorm:
					raise IOError("Could not get '" + normHistName + "' from '" + inFileName + "'")
				histListNorm.append(histNorm)

				indexHistName = sector+"_"+str(tBin)+"_index"
				histIndx = inFile.Get(indexHistName)
				if not histIndx:
					raise IOError("Could not get '" + indexHistName + "' from '" + inFileName + "'")
				histListIndx.append(histIndx)
			histsToFill = []
			comaHists = []
			for mBin in range(startBin, stopBin):
				comaHistName = "COMA_"+str(tBin)+"_"+str(mBin)
				comaHist = inFile.Get(comaHistName)
				if not comaHist:
					raise IOError("Could not get '" + comaHistName + "' from '" + inFileName + "'")
				comaHists.append(comaHist)
			ab = allBins(startBin, stopBin, histListReal, histListImag, histListNorm, histListIndx, comaHists)

			zeroCount = 0
			zeroHistList  = []
			eigenHistList = []
			while True:
				zeroName  = "zero"+str(zeroCount)+"_"+str(tBin)
				eigenName = "eigen"+str(zeroCount)+"_"+str(tBin)
				zeroHist = inFile.Get(zeroName)
				if not zeroHist:
					break
				print "Adding zero-mode" 
				zeroCount += 1
				if not zeroForSectors(sectors, zeroHist.GetTitle()):
					continue
				zeroHistList.append(zeroHist)

				eigenHist = inFile.Get(eigenName)
				if eigenHist:
					eigenHistList.append(eigenHist)
			if (not len(eigenHistList) == 0) and (not len(eigenHistList) == len(zeroHistList)):
				raise ValueError("Number of eigenvalue histograms does not match, but is also nonzero")
			removeCertainZeroModes(zeroHistList, eigenHistList)
			for zeroHist in zeroHistList:
				borders = getZeroHistBorders(zeroHist) 
				ab.addZeroMode(borders, zeroHist)
#			ab.rotateToPhaseOfBin(10)
			ab.removeZeroModeFromComa()

#			ab.unifyComa()

			ab.initChi2(waveModel)
			ab.setMassRanges(sectorRangeMap)

			totalModel.addBin(ab)

		totalPars = []
		for k in waveModel:
			for f in waveModel[k]:
				totalPars += f.getParameters()
		if not len(totalPars) == 0 and fitShape:
			res    = scipy.optimize.minimize(totalModel.chi2, totalPars)
			hi     = res.hess_inv
			errs   = []
			for i in range(len(res.x)):
				errs.append(hi[i,i]**.5)
			if not ab.setParametersAndErrors(res.x, errs):
				raise RuntimeError("Could not set parameter errors")
			print "The final Chi2 is:",res.fun

			for parameter in allParameters:
				print parameter



#		for k in waveModel:
#			for f in waveModel[k]:
#				print "function parameters",f.getParameters()
		for tBin in tBins:
			ab = totalModel[tBin]
			chi2, params = ab.chi2(returnParameters = True)
			paramsZ      = ab.linearizeZeroModeParameters(params)
			ab.setTheoryFromOwnFunctions(params, True)
			print "The final chi2 =",chi2

			intenses = []
			reals    = []
			imags    = []
			phases   = []
			intensD  = []
			realsD   = []
			imagsD   = []
			phasesD  = []
			intensT  = []
			realsT   = []
			imagsT   = []
			phasesT  = []

			for rh in histListReal:
				Ih = rh.Clone()
				Ih.Reset()
				intenses.append(Ih)

				realH = rh.Clone()
				realH.Reset()
				reals.append(realH)
	
				imagH = rh.Clone()
				imagH.Reset()
				imags.append(imagH)

				phaseH = rh.Clone()
				phaseH.Reset()
				phases.append(phaseH)

				ID = rh.Clone()
				ID.Reset()
				intensD.append(ID)
	
				rD = rh.Clone()
				rD.Reset()
				realsD.append(rD)

				iD = rh.Clone()
				iD.Reset()
				imagsD.append(iD)

				pD = rh.Clone()
				pD.Reset()
				phasesD.append(pD)

				IT = rh.Clone()
				IT.Reset()
				intensT.append(IT)

				rT = rh.Clone()
				rT.Reset()
				realsT.append(rT)

				iT = rh.Clone()
				iT.Reset()
				imagsT.append(iT)

				pT = rh.Clone()
				pT.Reset()
				phasesT.append(pT)

			ab.fillHistograms(paramsZ, intenses            )
			ab.fillHistograms(paramsZ, reals,  mode = REAL )
			ab.fillHistograms(paramsZ, imags,  mode = IMAG )
			ab.fillHistograms(paramsZ, phases, mode = PHASE)
			zeroP = [0.]*len(paramsZ)
			ab.fillHistograms(zeroP, intensD              )
			ab.fillHistograms(zeroP, realsD,  mode = REAL )
		        ab.fillHistograms(zeroP, imagsD,  mode = IMAG )
		        ab.fillHistograms(zeroP, phasesD, mode = PHASE)

			ab.fillHistograms(zeroP, intensT,  mode = INTENSTHEO)
			ab.fillHistograms(zeroP, realsT,   mode = REALTHEO  )
	       	        ab.fillHistograms(zeroP, imagsT,   mode = IMAGTHEO  )
	       	        ab.fillHistograms(zeroP, phasesT,  mode = PHASETHEO )

			for i in range(len(histListReal)):
				renormToBinWidth(intenses[i]  )
				renormToBinWidth(intensD[i]   )
				renormToBinWidth(intensT[i]   )
				renormToBinWidth(reals[i] , .5)
				renormToBinWidth(realsD[i], .5)
				renormToBinWidth(realsT[i], .5)
				renormToBinWidth(imags[i] , .5)
				renormToBinWidth(imagsD[i], .5)
				renormToBinWidth(imagsT[i], .5)

				allIsZero = True
				for binX in range(intensT[i].GetNbinsX()):
					for binY in range(intensT[i].GetNbinsY()):
						if not intensT[i].GetBinContent(binX+1, binY+1) == 0.:
							allIsZero = False
							break
					if not allIsZero:
						break # # # # # # # # # # # # # # # # # 
				if not allIsZero:
					rv = resultViewer([intenses[i], intensD[i], intensT[i]],[reals[i], realsD[i], realsT[i]],[imags[i], imagsD[i], imagsT[i]], [phases[i], phasesD[i], phasesT[i]])
				else:
					rv = resultViewer([intenses[i], intensD[i]],[reals[i], realsD[i]],[imags[i], imagsD[i]],[phases[i], phasesD[i]])
				rv.run()

if __name__ == "__main__":
	main()
