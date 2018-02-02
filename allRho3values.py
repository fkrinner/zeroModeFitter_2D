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

def main():

	folder      = "./rho3MassesAndWidthForScript"
	outFileName = "./allRhoFitResults.dat"
	mMin        =    1.2
	mMax        =    2.2
	GMin        =    0.1
	GMax        =    0.3
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
		print len(allVals),":::<<<><><>>>"


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

#	par     = [7.69640e-01, 1.21366e-02, 1.43298e-01, 2.20290e-02, 0.]
	par     = [1.6, .012, .2, .022, 0.]
#	fitFunc = twoDgausFunc(massWidthHist)
#	fitFunc.isLike = True

	ffm = oneDgaussFuncLike(mVals)
	ffG = oneDgaussFuncLike(GVals)
	likeFunc = twoDgaussFuncLike(allVals)
#	print fitFunc(par)
	res  = scipy.optimize.minimize(likeFunc, par)
	print res.x,2/pi*atan(res.x[-1])
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
