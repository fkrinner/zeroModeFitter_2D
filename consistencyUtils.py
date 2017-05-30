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
	totals = None
	for m in lst:
		diffSum = cloneZeros(diffs[m,m])
		if not totals:
			totals = cloneZeros(diffs[m,m])
		for n in lst:
			addBtoA(diffSum, diffs[n,m])
		diffSums[m] = diffSum
	for m in lst:
		invertAllEntries(diffSums[m])
	for m in lst:
		addBtoA(totals, diffSums[m])
	for m in lst:
		divideAbyB(diffSums[m], totals)
	return diffSums
	


