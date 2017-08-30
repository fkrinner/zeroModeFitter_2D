import scipy.optimize

def estimateErrors(func, minPoint, startErrors):
	if not len(startErrors) == 2*len(minPoint):
		raise ValueError("Dimension mismatch")
	stepSizes = [v for v in startErrors]
	x = [v for v in minPoint]
	c0 = func(minPoint)
	retVal = []
	for i in range(len(minPoint)):
		print "at parameter plus ",i
		count = 0
		minn  = 1.
		minI  = 0 
		while True:
			x[i]  += stepSizes[2*i]
			count += 1
			c1 = func(x)
			if abs(c1-c0-1.) < minn:
				minn = abs(c1-c0-1.)
				minI = count
#			out.write(str(i)+' '+ str(stepSizes[2*i]*count) + ' ' + str(abs(c1-c0-1.)) + '\n')
			if abs(c1-c0-1.) > 2.:
				break
			if count > 50 and count%1000 == 0:
				print count, abs(c1-c0-1.), minn, minI

		x[i] -= stepSizes[2*i]*count
		retVal.append(stepSizes[2*i]*minI)
		count = 0
		minn  = 1.
		maxI  = 0
		print "at parameter minus ",i
		while True:
			x[i]  -= stepSizes[2*i+1]
			count += 1
			c1 = func(x)
			if abs(c1-c0-1.) < minn:
				minn = abs(c1-c0-1.)
				maxI = count
#			out.write(str(i)+' '+ str(stepSizes[2*i]*count) + ' ' + str(abs(c1-c0-1.)) + '\n')
			if abs(c1-c0-1.) > 2.:
				break
			if count > 50 and count%1000 == 0:
				print count, abs(c1-c0-1.), minn, maxI

		x[i] += stepSizes[2*i+1]*count
		retVal.append(stepSizes[2*i+1]*maxI)
	return retVal


def estimateErrors2(func, pars, errsFromFitter, recursionLevel = 0, tolerance = 0.1, maxRecursionLevel = 10):
	return errsFromFitter
	centerEval = func(pars)
	if recursionLevel > maxRecursionLevel:
		print "Maximum recursion of",maxRecursionLevel,"exceeded, returning previous value"
		return errsFromFitter
	errs = []
	allEstimationInRange = True
	for p in range(len(pars)):
		pars[p]  += errsFromFitter[p]
		plusEval  = func(pars)
		pars[p]  -= 2*errsFromFitter[p]
		minusEval = func(pars)
		pars[p]  += errsFromFitter[p]

		errs.append(abs(errsFromFitter[p]*(-16 * centerEval + minusEval *(8+minusEval)-2*(minusEval-4) * plusEval + plusEval**2)**.5/(4*centerEval - 2*(minusEval + plusEval))))
		
		pars[p]  += errs[p]
		plusEval  = func(pars)
		pars[p]  -= 2*errs[p]
		minusEval = func(pars)
		pars[p]  += errs[p]
		estimator =  abs(plusEval+minusEval-2*centerEval-2)
		if estimator > 2*tolerance:
#			print 
#			print p," estimated ",abs(plusEval+minusEval-2*centerEval-2),"this is too large"
			allEstimationInRange = False
			if estimator > 10.:
				print "Too far off, scaling down, by:",estimator
				errs[-1] /= estimator

	if allEstimationInRange:
		return errs
	else:
		return estimateErrors2(func, pars, errs, recursionLevel = recursionLevel+1, tolerance = tolerance, maxRecursionLevel = maxRecursionLevel)


def main():
	def f(x):
		retVal = 0.
		for i,v in enumerate(x):
			retVal += v**2/(i+1)
		return retVal
	print estimateErrors2(f,[0.,0.,0.],[0.3, 0.3,0.3])


if __name__ == "__main__":
	main()
