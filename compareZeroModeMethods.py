from utils import sumUp
from analysisClass import amplitudeAnalysis
import parameterizationClasses as pc
import parameterTrackingParameterizations as ptc
from modes import PHASE, AMPL, SMOOTH, NONE

def doSeriesOfComparisons(models):
	"""
	use a map of name and model as input, to be forced to give names and thus know what you did
	"""
	for key in models:
		print "Fitting:",key
		models[key].fitParametersForMode()
	evals = {}
	for key1 in models:
		params = models[key1].getZeroModeParametersForMode()
		for key2 in models:
			evals[key2, key1] = models[key2].evaluateZeroModeParametersForMode(params)
	return evals
		
def main():
	inFileName       = "/nfs/mds/user/fkrinner/extensiveFreedIsobarStudies/results_std11.root"
	tBins            = [0]
	startBin         = 20
	stopBin          = 30

	mPi      = 0.13957018

	rhoMass  = ptc.parameter( .77549, "rhoMass" )
	rhoWidth = ptc.parameter( .1491 , "rhoWidth")
	rho      = ptc.relativisticBreitWigner([rhoMass,rhoWidth], mPi, mPi, mPi, 1, 0, False)
	rhoModel = amplitudeAnalysis(inFileName, ["1++0+[pi,pi]1--PiS"], {"1++0+[pi,pi]1--PiS":[rho]}, startBin, stopBin, tBins)
	rhoModel.loadData()
	rhoModel.finishModelSetup()
	rhoModel.mode = AMPL

	pipiSW    = pc.fixedParameterization("/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/anything/zeroModes/bwAmplitudes_noBF/amp_0mp0pSigmaPiS")
	pipiSphaseModel = amplitudeAnalysis(inFileName, ["1++0+[pi,pi]0++PiP"], {"1++0+[pi,pi]0++PiP":[pipiSW]}, startBin, stopBin, tBins)
	pipiSphaseModel.loadData()
	pipiSphaseModel.finishModelSetup()
	pipiSphaseModel.mode = PHASE

	smoothModel = amplitudeAnalysis(inFileName, ["1++0+[pi,pi]0++PiP", "1++0+[pi,pi]1--PiS"], {}, startBin, stopBin, tBins)
	smoothModel.loadData()
	smoothModel.finishModelSetup()
	smoothModel.mode = SMOOTH

	smoothModelP = amplitudeAnalysis(inFileName, ["1++0+[pi,pi]0++PiP"], {}, startBin, stopBin, tBins)
	smoothModelP.loadData()
	smoothModelP.finishModelSetup()
	smoothModelP.mode = SMOOTH

	smoothModelS = amplitudeAnalysis(inFileName, ["1++0+[pi,pi]1--PiS"], {}, startBin, stopBin, tBins)
	smoothModelS.loadData()
	smoothModelS.finishModelSetup()
	smoothModelS.mode = SMOOTH


	modellist = {
	        "rhoModel"        : rhoModel,
		"pipiSphaseModel" : pipiSphaseModel,
	        "smoothModel"     : smoothModel,
	        "smoothModelS"    : smoothModelS,
	        "smoothModelP"    : smoothModelP
	}
	evld = doSeriesOfComparisons(modellist)
	for key in modellist:
		norm = sumUp(evld[key,key])
		print
		for key2 in modellist:
			val = sumUp(evld[key, key2])
			print key2,'in',key,':',val/norm

if __name__ == "__main__":
	main()
