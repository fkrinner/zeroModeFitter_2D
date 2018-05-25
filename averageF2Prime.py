#!/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/Python_ultra/Python-2.7.10/bin/python
# averageF2Prime.py
# Created: 2018-05-08 15:15:16.757157
# Author: Fabian Krinner
import os, sys
import ROOT
def parseFile(inFileName):
	values = []
	with open(inFileName, 'r') as inFile:
		for line in inFile.readlines():
			vals = [float(v) for v in line.split()]
			values.append(vals)
	return values

def weightFunction(val):
	if val[2] < 1.8:
		return 0.
	if val[2] > 2.2:
		return 0.
	if val[4] > .5:
		return 0.
	if val[4] < .2:
		return 0.
	return 1.

def main():
	values = []
	for t in range(4):
		inFileName = "./2mp_f2primeMassesAndWidths_"+str(t)+".dat"
		values    += parseFile(inFileName)
#	mMax = max(v[2] for v in values)
#	mMin = min(v[2] for v in values)
#	Gmax = max(v[4] for v in values)
#	Gmin = min(v[4] for v in values)
	mMin = 1.
	mMax = 2.5
	Gmin = 0.
	Gmax = 1.

	print "h","h", 100, mMin, mMax, 100, Gmin, Gmax
	hist = ROOT.TH2D("h","h", 100, mMin, mMax, 100, Gmin, Gmax)
	for v in values:
		hist.Fill(v[2], v[4])
	hist.Draw("COLZ")
	raw_input()
#	return

	avgMass     = 0.
	avgWidth    = 0.
	cumulWeight = 0.
	for val in values:
		weight       = weightFunction(val)
		avgMass     += val[2]*weight
		avgWidth    += val[4]*weight
		cumulWeight += weight
	avgMass  /= cumulWeight
	avgWidth /= cumulWeight
	massVar   = 0.
	widthVar  = 0.
	for val in values:
		weight       = weightFunction(val)
		massVar     += (val[2]-avgMass)**2*weight
		widthVar    += (val[4]-avgWidth)**2*weight
	print "m =",avgMass,"+-",massVar**.5
	print "G =",avgWidth,"+-+", widthVar**.5

if __name__ == "__main__":
	main()
