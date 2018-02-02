# allRhovalues.py
import os, sys
import ROOT

from math import exp,pi,isnan,atan,log
import scipy.optimize

from utils import gaus2D, twoDgausFunc, twoDgaussFuncLike, oneDgaussFuncLike

from random import random

def isRegularError(err):
	if isnan(err) or err <= 0.:
		return False
	return True

def getValues(inFileName):
	vals = []
	with open(inFileName, 'r') as inFile:
		for line in inFile.readlines():
			flts = [float(v) for v in line.split()]
			val  = ((flts[2],flts[3]),(flts[4],flts[5]),flts[6])
			vals.append(val)
	return vals

def avgAndSig(vals):
	avg    = 0.
	weight = 0.
	for v in vals:
		if len(v) == 2:
			w = v[1]
		else:
			w = 1.
		avg += v[0]*w
		weight += w
	var


def avgAndComa(vals):
	avgs = [0.,0.]
	weight = 0.
	for v in vals:
		if len(v) == 3:
			w = v[2]
		else:
			w = 1.
		avgs[0] += w*v[0]
		avgs[1] += w*v[1]
		weight  += w
	avgs[0] /= weight
	avgs[1] /= weight
	coma = [[0.,0.],[0.,0.]]
	for v in vals:
		dels = [v[0] - avgs[0], v[1] - avgs[1]]
		if len(v) == 3:
			w = v[2]
		else:
			w = 1.
		for i in range(2):
			for j in range(2):
				coma[i][j] += dels[i]*dels[j]/weight
	return avgs, coma

def getTbin(fn):
	if 'global' in fn:
		raise ValueError("Global file encountered")
	for i in range(4):
		if "_"+str(i)+"_" in fn or "_"+str(i)+".dat" in fn:
			return i
	raise ValueError("Could't determine tBin for '" + fn + "'")

def main():

	folder      = "./rhoMassesAndWidthForScript"
	outFileName = "./allRhoFitResults.dat"
	mMin        =    .7
	mMax        =    .85
	GMin        =    .08
	GMax        =    .25
	chiMin      =   0.
	chiMax      =  20.
	nBins       = 100
	nBins2D     =  50

	massHist  = ROOT.TH1D("masses", "masses", nBins, mMin  , mMax  )
	widthHist = ROOT.TH1D("widths", "widths", nBins, GMin  , GMax  )
	chi2Hist  = ROOT.TH1D("chi2s" , "chi2s" , nBins, chiMin, chiMax)

	massWidthHist = ROOT.TH2D("2D", "2D", nBins2D, mMin, mMax, nBins2D, GMin, GMax)
#	with open(outFileName, 'w') as outFile:
	nTotal  = 0.
	allVals = []
	mVals   = []
	GVals   = []
	
	tBins = [0,1,2,3]
#	tBins = [0]
#	tBins = [1]
#	tBins = [2]
#	tBins = [3]

	maxErr = float('inf')
	for i in range(1,len(sys.argv)):
		if sys.argv[i].startswith("me:"):
			maxErr = float(sys.argv[i][3:])
			print "maxErr is set to", maxErr
		else:
			raise RuntimeError("Unknown option '"+sys.argv[i]+"'")

	if True:
		for fn in os.listdir(folder):
			if "global" in fn or "range" in fn or "_2mp_" in fn: # exlude "_2mp_" sine P and F are used separately
				continue
			if not getTbin(fn) in tBins:
				continue
			vals = getValues(folder + os.sep + fn)
			for v in vals:
				if v[0][1] > maxErr or v[1][1] > maxErr:
					continue
				if  v[0][0] > mMin and v[0][0] < mMax and v[1][0] > GMin and v[1][0] < GMax:
					nTotal =+ 1.
#					outFile.write(str(v[0][0])+' '+str(v[0][1])+' '+str(v[1][0])+' '+str(v[1][1])+' '+str(v[2])+'\n')
					allVals.append((v[0][0], v[1][0], 1./v[2]))
#					allVals.append((v[0][0], v[1][0], 1./(v[0][1]**2+v[1][1]**2)**.5))
					if not isRegularError((v[0][1]**2+v[1][1]**2)**.5) or not isRegularError((v[0][1]**2+v[1][1]**2)**.5):
						print v
#					mVals.append((v[0][0], 1./v[0][1]))
#					GVals.append((v[1][0], 1./v[1][1]))

					mVals.append([v[0][0],1./v[0][1]])

					GVals.append([v[1][0],random()])

				massHist.Fill(v[0][0])
				widthHist.Fill(v[1][0])
				chi2Hist.Fill(v[2])
				massWidthHist.Fill(v[0][0], v[1][0], 1./v[2])
		print len(allVals),":::<<<><><>>>",tBins


#	print " --- mass fit --- "
#	massHist.Fit("gaus")
#	print " - - - - - - - - "
#	massHist.Draw()
#	raw_input()
#	print " --- width fit --- "
#	widthHist.Fit("gaus")
#	print " - - - - - - - - - "
#	widthHist.Draw()
#	raw_input()
#	chi2Hist.Draw()
#	raw_input()

#	a,c = avgAndComa(allVals)
#	print [a[0], c[0][0]**.5, a[1], c[1][1]**.5, c[0][1]/(c[0][0]*c[1][1])**.5]

#	par     = [7.69640e-01, 1.21366e-02, 1.43298e-01, 2.20290e-02, 0.]
	par     = [.77, .012, .14, .022, 0.]
#	fitFunc = twoDgausFunc(massWidthHist)
#	fitFunc.isLike = True

	ffm = oneDgaussFuncLike(mVals)
	ffG = oneDgaussFuncLike(GVals)
	likeFunc = twoDgaussFuncLike(allVals)
	with open('allRhoValuesCombined.dat', 'w') as outFile:
		for v in allVals:
			outFile.write(str(v[0]) + ' ' + str(v[1]) + ' ' + str(v[2]) + '\n')
		

#	print fitFunc(par)
	res  = scipy.optimize.minimize(likeFunc, par)
	print res.x,2/pi*atan(res.x[-1])
#	print [(2*res.hess_inv[i,i])**.5 for i in range(len(res.x))]
#	resM = scipy.optimize.minimize(ffm, par[:2])
#	print resM.x
#	resG = scipy.optimize.minimize(ffG, par[2:4])
#	print resG.x
#	print ffm(par[ :2])
#	print ffG(par[2:4])
#	print likeFunc(res.x)
	
#	nh = fitFunc.makeDrawHist(res.x)
#	nh.Draw("COLZ")
#	massWidthHist.Draw("colz")
#	raw_input()


if __name__ == "__main__":
	main()
