import os

mapp = {
	"0-+0+[pi,pi]0++PiS" : ["amp_0mp0pf0980PiS","amp_0mp0pSigmaPiS","amp_0mp0pf01500PiS"],
	"0-+0+[pi,pi]1--PiP" : ["amp_0mp0pRhoPiP"],
	"1++0+[pi,pi]0++PiP" : ["amp_1pp0pf0980PiP","amp_1pp0pSigmaPiP"],
	"1++0+[pi,pi]1--PiS" : ["amp_1pp0pRhoPiS"],
	"1++1+[pi,pi]1--PiS" : ["amp_1pp1pRhoPiS"],
	"2-+0+[pi,pi]0++PiD" : ["amp_2mp0pf0980PiD","amp_2mp0pSigmaPiD"],
	"2-+0+[pi,pi]1--PiP" : ["amp_2mp0pRhoPiP"],
	"2-+0+[pi,pi]1--PiF" : ["amp_2mp0pRhoPiF"],
	"2-+0+[pi,pi]2++PiS" : ["amp_2mp0pf2PiS"],
	"2-+1+[pi,pi]1--PiP" : ["amp_2mp1pRhoPiP"],
	"2++1+[pi,pi]1--PiD" : ["amp_2pp1pRhoPiD"]}

def getFileNameForSector(sector, withBarrierFactors = False):
	if withBarrierFactors:
		folder = "/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/anything/zeroModes/bwAmplitudes"
	else:
		folder = "/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/anything/zeroModes/bwAmplitudes_noBF"
	return [folder + os.sep + fn for fn in mapp[sector]]

def main():
	for key in mapp:
		fns1 = getFileNameForSector(key, False)
		for fn1 in fns1:
			if not os.path.isfile(fn1):
				print fn1, "does not exist"
			else:
				print 't'
		fns2 = getFileNameForSector(key, True)
		for fn2 in fns2:
			if not os.path.isfile(fn2):
				print fn2, "does not exist"
			else:
				print 'q'

if __name__ == "__main__":
	main()
