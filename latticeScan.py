import analyze as a
import os

def main():
	folderName = "/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/anything/zeroModes/rhoShapes" 
	results = []
	bestChi2 = float('inf')
	bestMass  = 999.999
	bestWidth = 999.999
	for fn in os.listdir(folderName):
		if not os.path.isdir(folderName + os.sep + fn):
			raise IOError("folder does not exist")
		mass = float(fn)
		for pn in os.listdir(folderName + os.sep + fn):
			width = float(pn)
			fileName = folderName + os.sep + fn + os.sep + pn
			if not os.path.isfile(fileName):
				raise IOError("file does not exist '" + fileName + "'")
			chi2 = a.main(fileName)
			if chi2 < bestChi2:
				bestChi2  = chi2
				bestMass  = mass
				bestWidth = width
			print mass, width, chi2
			results.append((chi2, mass, width))
			print "best so far", bestChi2, bestMass, bestWidth
	results.sort()
	print results
	print results[0]

if __name__ == "__main__":
	main()


