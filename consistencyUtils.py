from utils import cloneZeros, addBtoA, divideAbyB, invertAllEntries

def makeStructureDiff(A, B, retVal):
	if not len(A) == len(B) or not len(A) == len(retVal):
		raise ValueError("Size mismatch A/B/retVal")
	for i in range(len(A)):
		if hasattr(A[i], "__len__"):
			makeStructureDiff(A[i], B[i], retVal[i])
		else:
			retVal[i] = ((A[i] - B[i])/B[i]).real

def getmBinResolvedDiffs(methodMap):
	params    = {} 
	selfEvals = {}
	for m in methodMap:
		params[m]    = methodMap[m].getZeroModeParametersForMode()
		selfEvals[m] = methodMap[m].evaluateResolvedZeroModeParametersForMode(params[m])
	diff = {}
#	for m in params:
#		print m, params[m]
#	raise Exception
	for m in methodMap:
		for n in methodMap:
			evals = methodMap[m].evaluateResolvedZeroModeParametersForMode(params[n])
			diff[m,n] = cloneZeros(evals)
			makeStructureDiff(evals, selfEvals[m], diff[m,n])
	return diff

def getCompositions(diffs):
	lst = []
	for pair in diffs:
		if not pair[0] in lst:
			lst.append(pair[0])
	diffSums = {}
	totals   = None
	for m in lst:
		diffSum = cloneZeros(diffs[m,m])
		if not totals:
			totals = cloneZeros(diffs[m,m])
		for n in lst:
			addBtoA(diffSum, diffs[n,m], noBelowZero = True)
		diffSums[m] = diffSum
	for m in lst:
		invertAllEntries(diffSums[m], nanToZero = True)

	for m in lst:
		addBtoA(totals, diffSums[m])
	for m in lst:
		divideAbyB(diffSums[m], totals, ignoreZero = True)
	return diffSums
	
def doAllComparisons(methodMap, startBin, methodBinRanges = {}, noBelowZero = False):
	"""
	Does everything theo toher methods do, bit with the possibility of restricted ranged for several methods
	Start bin has to be given to match the bin ranges, parameters are just counted in the methods
	"""
	params = {}
	evals  = {}
	for m in methodMap:
		params[m] = methodMap[m].getZeroModeParametersForMode()

#		for mb in params[m][0]:
#			print len(mb),
#		print "<-----",m

		for n in methodMap:
			evals[n,m] = methodMap[n].evaluateResolvedZeroModeParametersForMode(params[m])
	diffs         = {}
	resolvedDiffs = {}
	for m in methodMap:
		for n in methodMap:
			sumM = 0.
			sumN = 0.
			resolvedDiffs[m,n] = []
			for tBin in range(len(evals[m,n])):
				resolvedDiffs[m,n].append([])
				for mBin in range(len(evals[m,n][tBin])):
					nBin = mBin + startBin
					if m in methodBinRanges:
						if nBin < methodBinRanges[m][0] or nBin >= methodBinRanges[m][1]:
							resolvedDiffs[m,n][tBin].append(float('nan'))
							continue
					if n in methodBinRanges:
						if nBin < methodBinRanges[n][0] or nBin >= methodBinRanges[n][1]:
							resolvedDiffs[m,n][tBin].append(float('nan'))
							continue
					resolvedDiffs[m,n][tBin].append((evals[m,n][tBin][mBin]-evals[m,m][tBin][mBin])/evals[m,m][tBin][mBin])
					sumM += evals[m,m][tBin][mBin]
					sumN += evals[m,n][tBin][mBin]
			diffs[m,n] = (sumN - sumM)/sumM
	for m in methodMap:
		resolvedWA    = cloneZeros(params[m])
		nonResolvedWA = cloneZeros(params[m])
		break
	compositions = {}
	for m in methodMap:
		compositions[m] = [[0.] * len(resolvedWA[i]) for i in range(len(resolvedWA))]
	for tBin in range(len(resolvedWA)):
		for mBin in range(len(resolvedWA[tBin])):
			totalWeight    = 0.
			totalWeightNon = 0.
			for m in methodMap:
				nBin = mBin + startBin
				if m in methodBinRanges:
					if nBin < methodBinRanges[m][0] or nBin >= methodBinRanges[m][1]:
						continue
				totalDiff    = 0.
				totalDiffNon = 0.
				for n in methodMap:
					if n in methodBinRanges:
						if nBin < methodBinRanges[n][0] or nBin >= methodBinRanges[n][1]:
							continue
					toadd = (evals[n,m][tBin][mBin] - evals[n,n][tBin][mBin])/evals[n,n][tBin][mBin]
					if not toadd < 0. or not noBelowZero:
						totalDiff    += toadd
					totalDiffNon += diffs[n,m]
				for i in range(len(resolvedWA[tBin][mBin])):

					resolvedWA[tBin][mBin][i] += params[m][tBin][mBin][i]/totalDiff
					nonResolvedWA[tBin][mBin][i] += params[m][tBin][mBin][i]/totalDiffNon
				totalWeight    += 1./totalDiff
				totalWeightNon += 1./totalDiffNon
				compositions[m][tBin][mBin] = 1./totalDiff
			for i in range(len(resolvedWA[tBin][mBin])):
				resolvedWA[tBin][mBin][i]/=totalWeight
				nonResolvedWA[tBin][mBin][i]/=totalWeightNon
			for m in methodMap:
				compositions[m][tBin][mBin]/=totalWeight
				
	resolvedWAdiffs    = {}
	nonResovledWAdiffs = {}
	noCorrDiffs        = {}
	for m in methodMap:
		resEvals = methodMap[m].evaluateResolvedZeroModeParametersForMode(resolvedWA)
		nonEvals = methodMap[m].evaluateResolvedZeroModeParametersForMode(nonResolvedWA)
		zerEvals = methodMap[m].evaluateResolvedZeroModeParametersForMode(cloneZeros(nonResolvedWA))
		aSumRes = 0.
		aSumNon = 0.
		aSumZer = 0.
		oSum    = 0.
		resolvedDiffs[m,"WAres"] = []
		resolvedDiffs[m,"WAnon"] = []
		for tBin in range(len(resolvedWA)):
			resolvedDiffs[m,"WAres"].append([])
			resolvedDiffs[m,"WAnon"].append([])
			for mBin in range(len(resolvedWA[tBin])):
				nBin = mBin + startBin
				if m in methodBinRanges:
					if nBin < methodBinRanges[m][0] or nBin >= methodBinRanges[m][1]:
						resolvedDiffs[m,"WAres"][tBin].append(float('nan'))
						resolvedDiffs[m,"WAnon"][tBin].append(float('nan'))
						continue
				resolvedDiffs[m,"WAres"][tBin].append((resEvals[tBin][mBin]-evals[m,m][tBin][mBin])/evals[m,m][tBin][mBin])
				resolvedDiffs[m,"WAnon"][tBin].append((nonEvals[tBin][mBin]-evals[m,m][tBin][mBin])/evals[m,m][tBin][mBin])
				aSumRes += resEvals[tBin][mBin]
				aSumNon += nonEvals[tBin][mBin]
				aSumZer += zerEvals[tBin][mBin]
				oSum    += evals[m,m][tBin][mBin]
		resolvedWAdiffs[m]    = (aSumRes-oSum)/oSum
		nonResovledWAdiffs[m] = (aSumNon-oSum)/oSum
		noCorrDiffs[m]        = (aSumZer-oSum)/oSum
	return diffs, resolvedWA, nonResolvedWA, compositions, resolvedWAdiffs, nonResovledWAdiffs, resolvedDiffs,noCorrDiffs

