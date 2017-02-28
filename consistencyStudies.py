from waveNameClass import waveName
from rootfabi import root_open, GetKeyNames
import numpy as np
import numpy.linalg as la
from allBinsClass import allBins
from massBinClass import massBin
from modes import REAL, IMAG, PHASE, INTENS, INTENSNORM, INTENSTHEO, REALTHEO, IMAGTHEO, PHASETHEO
from resultViewerClass import resultViewer
import os, sys
from utils import zeroForSectors, getZeroHistBorders, renormToBinWidth
import parameterizationClasses as pc
import scipy.optimize
from removeZeroModes import removeCertainZeroModes
from fixedparameterizationPaths import getFileNameForSector

def main():
	inFileName = "/nfs/mds/user/fkrinner/extensiveFreedIsobarStudies/results_std11.root"
#	sectors    = ["1++0+[pi,pi]0++PiP","1++0+[pi,pi]1--PiS"]
#	sectors    = ["2-+0+[pi,pi]0++PiD","2-+0+[pi,pi]1--PiP","2-+0+[pi,pi]1--PiF","2-+0+[pi,pi]2++PiS"]
#	sectors    = ["2-+0+[pi,pi]1--PiP"]
#	sectors    = ["1-+1+[pi,pi]1--PiP"]
#	sectors    = ["1++0+[pi,pi]0++PiP"]
#	sectors    = ["2++1+[pi,pi]1--PiD"]
	sectors    = ["0-+0+[pi,pi]0++PiS", "0-+0+[pi,pi]1--PiP"]
#	sectors    = ["0-+0+[pi,pi]1--PiP"]
#	sectors    = ["0-+0+[pi,pi]0++PiS"]
#	sectors    = ["2-+0+[pi,pi]2++PiS"]
#	sectors    = ["2-+1+[pi,pi]1--PiP"]
#	sectors    = ["1++0+[pi,pi]1--PiS","1++0+[pi,pi]0++PiP"]
#	sectors    = ["1++1+[pi,pi]1--PiS"]
#	sectors    = ["1++0+[pi,pi]1--PiS"]
#	sectors    = ["3++0+[pi,pi]3--PiS"]
#	sectors    = ["3++0+[pi,pi]2++PiP"]
#	sectors    = ["3++0+[pi,pi]1--PiD"]
#	sectors    = ["3++0+[pi,pi]3--PiS","3++0+[pi,pi]2++PiP","3++0+[pi,pi]1--PiD"]
#	sectors    = ["4++1+[pi,pi]1--PiG", "4++1+[pi,pi]2++PiF"]
#	sectors    = ["2++1+[pi,pi]1--PiD","2++1+[pi,pi]2++PiP"]

#	sectors    = [sectors[1]]

	doSpecialOneBinFit = -15 # negative values turn it off
	sectorUseMap       = { # Defines, if for the given sector a theory curve will be used
		"0-+0+[pi,pi]0++PiS" : False,
		"0-+0+[pi,pi]1--PiP" : True,
		"1++0+[pi,pi]0++PiP" : True, 
		"1++0+[pi,pi]1--PiS" : True, 
		"2-+0+[pi,pi]0++PiD" : True, 
		"2-+0+[pi,pi]1--PiP" : True, 
		"2-+0+[pi,pi]1--PiF" : True, 
		"2-+0+[pi,pi]2++PiS" : True
	}

	sectorRangeMap = {}
	outFileName    = "resultFilesConsistency/out"
	modeDef        = False
	special        = False
	useSmooth      = False
	fitShape       = True
	for i in range(1,len(sys.argv)):
		arg = sys.argv[i]
		if arg.startswith('t:'):
			tBin = int(arg[2:])
			if tBin < 0 or tBin > 3:
				raise ValueError("Invalid t' bin: '" +str(tBin)+"'")
			outFileName += "_t"+str(tBin)
		elif arg.startswith('b:'):
			chunks = arg.split(':')
			startBin = int(chunks[1])
			if len(chunks) > 2:
				stopBin = int(chunks[2])
			else:
				stopBin = startBin +1
			if startBin < 0 or stopBin < 1 or startBin > 49 or stopBin > 50:
				raise ValueError("Invalid m bins: '" + str(startBin) + "' and '" + str(stopBin) + "'")
			outFileName += "_m"+str(startBin)+'-'+str(stopBin)
		elif arg.startswith("u:"):
			use = [int(c) for c in arg.split(':')[1:]]
			outFileName += "_u"
			for i, s in enumerate(sectors):
				if i in use:
					sectorUseMap[s] = True
					outFileName += str(i)
				else:
					sectorUseMap[s] = False
		elif arg.startswith("r:"):
			chunks = arg.split(":")[1:]
			if not len(chunks)%3 == 0:
				raise ValueError("Range definitions invalid")
			for i in range(len(chunks)/3):
				i    = int(chunks[3*i])
				mMin = float(chunks[3*i+1])
				mMax = float(chunks[3*i+2])
				outFileName += "_r"+str(i)+"-"+str(mMin)+"-"+str(mMax)
				sectorRangeMap[sectors[i]] = (mMin, mMax)
		elif arg == '-f':
			if modeDef:
				raise RuntimeError("Mode already defined")
			outFileName += "_fixedShapes"
			modelMode    = "fixedShapes"
			modeDef      = True
		elif arg == '-bw':
			if modeDef:
				raise RuntimeError("Mode already defined")
			outFileName += "_BW"
			modelMode    = "BW"			
			modeDef      = True
		elif arg =='-S':
			if modeDef:
				raise RuntimeError("Mode already defined")
			outFileName += "_smooth"
			useSmooth    = True
			modeDef      = True
		elif arg == "-P":
			special = True
			if not stopBin - startBin == 1:
				raise ValueError("Phase fit only possible for a single bin")
			outFileName += "_phase"
		else:
			raise ValueError("Unknown argument: '" + arg + "'")

	automatic = len(sys.argv) > 1
	if automatic:
		print outFileName
	else:
		modelMode = ""
		tBin      = 0
		startBin  = 31
		stopBin   = 32

	if modelMode == "BW" or modelMode == "2BW":
#		rho = pc.breitWigner()
		rho = pc.rpwaBreitWignerInt(0.13957018,0.13957018,0.13957018,1,1)
#		rhoParameters = [.77549, .1491]
		rhoParameters = [1.2, .2]
		rho.setParameters(rhoParameters)
#		rho.weightIntens = True
#		rhoParameters = [1.27549, .100]
		rhoPrime = pc.breitWigner()
		rhoPrimeParameters = [1.570, 0.144]
		rhoPrime.setParameters(rhoPrimeParameters)
		rhoModel = [rho]
		if modelMode == "2BW":
			rhoModel.append(rhoPrime)
		waveModel = {}
		for sector in sectors:
			if '1--' in sector:
				waveModel[sector] = rhoModel
			if "3--" in sector:
				waveModel[sector] = rhoModel

	elif modelMode == "fixedShapes":
		useBF       = False
		merge0pp    = False
		polyDegree  = 0
		polyComplex = True
		waveModel   = {}
		for sector in sectors:
#			if "3--" in sector:
#				rho = pc.rpwaBreitWignerInt(0.13957018,0.13957018,0.13957018,1,1)
#				rho.setParameters([1.2, .2])
#				model = [rho]
#				waveModel[sector] = model
#				continue
			if sector in sectorUseMap:
				if not sectorUseMap[sector]:
					continue
			model = []
			fileNames = getFileNameForSector(sector, useBF, merge0pp)
			print fileNames
			for fn in fileNames:
				param = pc.fixedParameterization(fn, polynomialDegree  = polyDegree, complexPolynomial = polyComplex)
				model.append(param)
			waveModel[sector] = model

	with root_open(inFileName, "READ") as inFile:
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
#		ab.rotateToPhaseOfBin(10)
		ab.removeZeroModeFromComa()

#		ab.unifyComa()

		if special:
			from random import random
			ab.setMassRanges(sectorRangeMap)
#			with open(outFileName, 'w') as out:
			if True:
				for mb in ab.massBins:
					mb.setZeroTheory()
				ab.initChi2(waveModel)
				for i,nBin in enumerate(range(startBin, stopBin)):
					print "at bin", nBin
					leBin = ab.massBins[i]
					pars  = [random() for _ in range(leBin.nParAll())]
					print leBin.phaseChi2(pars)
					res = scipy.optimize.minimize(leBin.phaseChi2, pars)
					print "phaseFit for bin",nBin 
					print "Giving a Chi2 of:",res.fun
				paramsZ = [res.x[0], res.x[1]] * 100 # Very, very bad hack...
#			sys.exit(0)	

		if doSpecialOneBinFit >= 0:
			specialBin = ab.massBins[doSpecialOneBinFit]
			specialBin.initChi2(waveModel)
			print "Doing special fit for: ",+ specialBin.bin3pi
			pars = [.77549, .1491]
			if modelMode == "explicitRhoFile":
				return specialBin.chi2([])
			res = scipy.optimize.minimize(specialBin.chi2, pars)
			hi = res.hess_inv
			print "m0 = ",res.x[0],"+-",hi[0,0]**.5
			print "G0 = ",res.x[1],"+-",hi[1,1]**.5
			print "Giving a Chi2 of:",res.fun

			sys.exit(0)

		if not useSmooth and not special and automatic:
			ab.initChi2(waveModel)
			ab.setMassRanges(sectorRangeMap)
			totalPars = []
			for k in waveModel:
				for f in waveModel[k]:
					totalPars += f.getParameters()
			if not len(totalPars) == 0 and fitShape:
				res = scipy.optimize.minimize(ab.chi2, totalPars)
				hi = res.hess_inv
				print "m0 = ",res.x[0],"+-",hi[0,0]**.5
				print "G0 = ",res.x[1],"+-",hi[1,1]**.5
				print "Giving a Chi2 of:",res.fun
			for k in waveModel:
				for f in waveModel[k]:
					print "function parameters",f.getParameters()
			chi2, params = ab.chi2(returnParameters = True)
			paramsZ      = ab.linearizeZeroModeParameters(params)
			print params
#			with open(outFileName, 'w') as out:
#				for i, b in enumerate(range(startBin, stopBin)):
#					out.write(str(b) + ' ' + str(paramsZ[2*i]) + ' ' + str(paramsZ[2*i+1]) + '\n')
#			print ab.massBins[0].phaseChi2(params[0])
# # ### ### ### ### ### ### ### 						#
# # It seems, that the phase fit is "not so good", since the COMA does something#
# # weird. The \deltas for the normal fit result fed into phaseChi2() are a 	#
# # factor 4 smaller, but after the COMA, the Chi2 is bigger than the one of the#
# # Phase fit									#
# # ### ### ### ### ### ### ### 						#
#			sys.exit(0)
		
			ab.setTheoryFromOwnFunctions(params, True)
			print "The final chi2 =",chi2
		elif not special and automatic:
			A,B,C   =  ab.getSmoothnessABC()
			paramsZ = -np.dot(la.inv(A + np.transpose(A)), B)
			with open(outFileName, 'w') as out:
				for i, b in enumerate(range(startBin, stopBin)):
					out.write(str(b) + ' ' + str(paramsZ[2*i]) + ' ' + str(paramsZ[2*i+1]) + '\n')
#			sys.exit(0)
		if not automatic:
			print "Automatic override"
			paramsZ = [0.] * 200
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

		if not useSmooth and automatic:
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
			if automatic:
				startCommand = "wq:"+outFileName+"_sector"+str(i)
			else:
				startCommand = ""
			if not allIsZero:
				rv = resultViewer([intenses[i], intensD[i], intensT[i]],[reals[i], realsD[i], realsT[i]],[imags[i], imagsD[i], imagsT[i]], [phases[i], phasesD[i], phasesT[i]], startBin = startBin, startCommand = startCommand)
			else:
				rv = resultViewer([intenses[i], intensD[i]],[reals[i], realsD[i]],[imags[i], imagsD[i]],[phases[i], phasesD[i]], startBin = startBin, startCommand = startCommand)
			rv.run()

if __name__ == "__main__":
	main()
