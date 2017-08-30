import os

mapp = {
	"0-+0+[pi,pi]0++PiS" : ["amp_0mp0pf0980PiS","amp_0mp0pSigmaPiS","amp_0mp0pf01500PiS"],
	"0-+0+[pi,pi]1--PiP" : ["amp_0mp0pRhoPiP" ],
	"1++0+[pi,pi]0++PiP" : ["amp_1pp0pf0980PiP","amp_1pp0pSigmaPiP"],
	"1++0+[pi,pi]1--PiS" : ["amp_1pp0pRhoPiS" ],
	"1++0+[pi,pi]1--PiD" : ["amp_1pp0pRhoPiD" ],
	"1++0+[pi,pi]2++PiP" : ["amp_1pp0pf2PiP"  ],
	"1++1+[pi,pi]1--PiS" : ["amp_1pp1pRhoPiS" ],
	"1-+1+[pi,pi]1--PiP" : ["amp_1mp1pRhoPiP" ],
	"2-+0+[pi,pi]0++PiD" : ["amp_2mp0pf0980PiD","amp_2mp0pSigmaPiD"],
	"2-+0+[pi,pi]1--PiP" : ["amp_2mp0pRhoPiP" ],
	"2-+0+[pi,pi]1--PiF" : ["amp_2mp0pRhoPiF" ],
	"2-+0+[pi,pi]2++PiS" : ["amp_2mp0pf2PiS"  ],
	"2-+0+[pi,pi]2++PiD" : ["amp_2mp0pf2PiD"  ],
	"2-+1+[pi,pi]1--PiP" : ["amp_2mp1pRhoPiP" ],
	"2++1+[pi,pi]1--PiD" : ["amp_2pp1pRhoPiD" ],
	"2++1+[pi,pi]2++PiP" : ["amp_2pp1pf2PiP"  ],
	"3++0+[pi,pi]1--PiD" : ["amp_3pp0pRhoPiD" ],
	"3++0+[pi,pi]2++PiP" : ["amp_3pp0pf2PiP"  ],
	"3++0+[pi,pi]3--PiS" : ["amp_3pp0pRho3PiS"],
	"4++1+[pi,pi]1--PiG" : ["amp_4pp1pRhoPiG" ],
	"4++1+[pi,pi]2++PiF" : ["amp_4pp1pf2PiF"  ], 
	"4-+0+[pi,pi]1--PiF" : ["amp_4mp0pRhoPiF" ],
	"6-+0+[pi,pi]1--PiH" : ["amp_6pm0pRhoPiH" ]}

mappMerged = {
	"0-+0+[pi,pi]0++PiS" : ["amp_0mp0p0ppPiS" ],
	"0-+0+[pi,pi]1--PiP" : ["amp_0mp0pRhoPiP" ],
	"1++0+[pi,pi]0++PiP" : ["amp_1pp0p0ppPiP" ],
	"1++0+[pi,pi]1--PiS" : ["amp_1pp0pRhoPiS" ],
	"1++0+[pi,pi]1--PiD" : ["amp_1pp0pRhoPiD" ],
	"1++0+[pi,pi]2++PiP" : ["amp_1pp0pf2PiP"  ],
	"1++1+[pi,pi]1--PiS" : ["amp_1pp1pRhoPiS" ],
	"1-+1+[pi,pi]1--PiP" : ["amp_1mp1pRhoPiP" ],
	"2-+0+[pi,pi]0++PiD" : ["amp_2mp0p0ppPiD" ],
	"2-+0+[pi,pi]1--PiP" : ["amp_2mp0pRhoPiP" ],
	"2-+0+[pi,pi]1--PiF" : ["amp_2mp0pRhoPiF" ],
	"2-+0+[pi,pi]2++PiS" : ["amp_2mp0pf2PiS"  ],
	"2-+0+[pi,pi]2++PiD" : ["amp_2mp0pf2PiD"  ],
	"2-+1+[pi,pi]1--PiP" : ["amp_2mp1pRhoPiP" ],
	"2++1+[pi,pi]1--PiD" : ["amp_2pp1pRhoPiD" ],
	"2++1+[pi,pi]2++PiP" : ["amp_2pp1pf2PiP"  ],
	"3++0+[pi,pi]1--PiD" : ["amp_3pp0pRhoPiD" ],
	"3++0+[pi,pi]2++PiP" : ["amp_3pp0pf2PiP"  ],
	"3++0+[pi,pi]3--PiS" : ["amp_3pp0pRho3PiS"],
	"4++1+[pi,pi]1--PiG" : ["amp_4pp1pRhoPiG" ],
	"4++1+[pi,pi]2++PiF" : ["amp_4pp1pf2PiF"  ],
	"4-+0+[pi,pi]1--PiF" : ["amp_4mp0pRhoPiF" ],
	"6-+0+[pi,pi]1--PiH" : ["amp_6pm0pRhoPiH" ]}


def getFileNameForSector(sector, withBarrierFactors = False, merged = False):
	if withBarrierFactors:
		folder = "/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/anything/zeroModes/bwAmplitudes"
	else:
		folder = "/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/anything/zeroModes/bwAmplitudes_noBF"
	if merged:
		return [folder + os.sep + fn for fn in mappMerged[sector]]	
	else:
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
