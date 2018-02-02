import os


def checkFitResult(vals, massLimits = None, widthLimits = None):
	if all(v == 0. for v in vals):
		return False
	if massLimits:
		if massLimits[0] > vals[0]:
			return False
		if massLimits[1] < vals[0]:
			return False
	if widthLimits:
		if widthLimits[0] > vals[2]:
			return False
		if widthLimits[1] < vals[2]:
			return False
	return True

def readFile(inFileName, massLimits = None, widthLimits = None):
	retVal = [None]*50
	with open(inFileName, 'r') as inFile:
		for line in inFile.readlines():
			chunks   = line.split()
			binIndex = int(chunks[0])
			vals     = [float(v) for v in chunks[2:]]
			if not checkFitResult(vals, massLimits, widthLimits):
				continue
			retVal[binIndex] = vals
	return retVal

def main():
	folder  = "./massesAndWidths"
	outFile = open("./allRhoMasses.dat", 'w')
	count   = 0
	vals    = []
	massLimits = (.7,0.95)
	widthLimits = (0.05, 0.3)

#	massLimits  = None
#	widthLimits = None

	for fn in os.listdir(folder):
		if not "rhoMassesAndWidths" in fn:
			continue
		inFileName = folder + os.sep + fn
		vals  += readFile(inFileName, massLimits = massLimits, widthLimits = widthLimits)
		
	mMean = 0.
	Gmean = 0.
	mSys  = 0.
	Gsys  = 0.
	for i,v in enumerate(vals):
		if v:
			count += 1
			outFile.write(str(i) + ' ' + str(v[0]) + ' ' + str(v[1]) + ' ' + str(v[2]) + ' ' + str(v[3]) + '\n')
			mMean += v[0]
			Gmean += v[2]
			mSys  += v[1]**2
			Gsys  += v[3]**2
#	return
	mMean /= count
	Gmean /= count
	mStat  = 0.
	Gstat  = 0.
	corr   = 0.

	for v in vals:
		if v:
			mStat += (v[0] - mMean)**2
			Gstat += (v[2] - Gmean)**2
			corr  += (v[0] - mMean)*(v[2] - Gmean)
	mStat /= count-1
	Gstat /= count-1
	corr  /= count-1 * mStat**.5 *Gstat**.5

	print "m = "+str(mMean)+" +- "+str(mSys**.5)+ "sys +- "+str(mStat**.5)+"stat"
	print "G = "+str(Gmean)+" +- "+str(Gsys**.5)+ "sys +- "+str(Gstat**.5)+"stat"
	print "corr =",corr


		
	print count,"masses and width written"

if __name__ == "__main__":
	main()

