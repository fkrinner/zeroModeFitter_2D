# allRhovalues.py
import os, sys
import ROOT

from math import exp,pi,isnan,atan,log
import scipy.optimize

from utils import gaus2D, twoDgausFunc, twoDgaussFuncLike, oneDgaussFuncLike, threeDgaussFuncLike

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
			val  = ((flts[2],flts[3]),(flts[4],flts[5]),(flts[6],flts[7]),flts[8])
			vals.append(val)
	return vals

def main():

	folder      = "./f0MassesAndWidthForScript"
	mMin        =   .8
	mMax        =   1.2
	g1Min       =    .0
	g1Max       =    .2
	g2Min       =   -.2
	g2Max       =    .5
	chiMin      =   0.
	chiMax      =  20.
	nBins       = 100
	nBins2D     =  50

	massHist  = ROOT.TH1D("masses", "masses", nBins, mMin  , mMax  )
	g1Hist    = ROOT.TH1D("g1", "g1", nBins, g1Min  , g1Max  )
	g2Hist    = ROOT.TH1D("g2", "g2", nBins, g2Min  , g2Max  )
	chi2Hist  = ROOT.TH1D("chi2s" , "chi2s" , nBins, chiMin, chiMax)

#	with open(outFileName, 'w') as outFile:
	nTotal  = 0.
	allVals = []
	mVals   = []
	GVals   = []

	fit = '012'

	maxErr = float('inf')
	for i in range(1,len(sys.argv)):
		if sys.argv[i].startswith("me:"):
			maxErr = float(sys.argv[i][3:])
			print "maxErr is set to", maxErr
		else:
			raise RuntimeError("Unknown option '"+sys.argv[i]+"'")

	if True:
		for fn in os.listdir(folder):
			if "global" in fn or "range" in fn: # exlude "_2mp_" sine P and F are used separately
				continue
			vals = getValues(folder + os.sep + fn)
			for v in vals:
				if v[0][1] > maxErr or v[1][1] > maxErr or v[2][1] > maxErr:
					continue

				if v[0][0] ==0.965 and v[1][0] == 0.165 and v[2][0] == 0.69465:
					print ":::::: ERRRRRRORRRRR"
					continue

				if  v[0][0] > mMin and v[0][0] < mMax and v[1][0] > g1Min and v[1][0] < g1Max and v[2][0] > g2Min and v[2][0] < g2Max:
					nTotal =+ 1.
					if fit == '01':
						allVals.append((v[0][0], v[1][0], 1./v[3]))
					elif fit == '02':
						allVals.append((v[0][0], v[2][0], 1./v[3]))
					elif fit == '12':
						allVals.append((v[1][0], v[2][0], 1./v[3]))
					elif fit == '012':
						allVals.append((v[0][0], v[1][0], v[2][0], 1./v[3]))
#					allVals.append((v[0][0], v[1][0], 1./v[3]))

				massHist.Fill(v[0][0])
				g1Hist.Fill(v[1][0])
				g2Hist.Fill(v[2][0])
				chi2Hist.Fill(v[3])
#				massWidthHist.Fill(v[0][0], v[1][0], 1./v[2])
		print len(allVals),":::<<<><><>>>"

#	par     = [9.77354e-01, 1.02718e-02, 6.19730e-02, 2.15549e-02, 4.40394e-02, 1.06449e-01, 0., 0., 0.]
	if fit == '01':
		par = [0.97960044,  0.01150727,  0.06686056,  0.02862707, -0.25047197]
	elif fit == '02':
		par = [0.97960044,  0.01150728,  0.0727524,   0.13166826,  0.81707373]
	elif fit == '12':
		par = [ 0.06686055,  0.02862706,  0.07275247,  0.13166824,  0.68740522]
	elif fit == '012':
#		par = [0.97960044,  0.01150727,  0.06686056,  0.02862707, 0.0727524 , 0.13166826 , -0.25047197, 0.81707373, 0.68740522]
		par = [0.97960044,  0.01150727,  0.06686055,  0.02862706, 0.07275241,  0.13166822, -0.25047208, 0.81707416,  0.68740525]

#	fitFunc = twoDgausFunc(massWidthHist)
#	fitFunc.isLike = True

#	with open("f0_parameters.dat", 'w') as outFile:
#		for v in allVals:
#			outFile.write(str(v[0]) + ' ' + str(v[1]) + ' ' + str(v[2]) + '\n')
#	return

	if fit in ['01','02','12']:
		likeFunc = twoDgaussFuncLike(allVals)
	elif fit == '012':
		likeFunc = threeDgaussFuncLike(allVals)
#	likeFunc = twoDgaussFuncLike(allVals)
	
	res  = scipy.optimize.minimize(likeFunc, par)
	print res.x
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
