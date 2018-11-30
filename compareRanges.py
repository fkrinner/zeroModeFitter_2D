#!/usr/bin/python
# compareRanges.py
# Created: 2018-11-07 15:24:43.786462
# Author: Fabian Krinner
import os, sys
from rootfabi import root_open
import numpy as np
import numpy.linalg as la

def getHistos(inFileName):
	with root_open(inFileName,"READ") as inFile:

		histRe = inFile.Get("real0")
		if not histRe:
			raise IOError("could not get real part from '" + inFileName + "'")
		histRe.SetDirectory(0)

		histIm = inFile.Get("imag0")
		if not histIm:
			raise IOError("Could not get imag part from '" + inFileName + "'")
		histIm.SetDirectory(0)

		histCr = inFile.Get("reImCorrel")
		if not histCr:
			raise IOError("Could not get correlation from '" + inFileName + "'")
		histCr.SetDirectory(0)

	return histRe, histIm, histCr

def getComas(re, im, correl):
	comaInvs = []
	for i in xrange(re.GetNbinsY()):
		for j in xrange(im.GetNbinsY()):
#			coma      = np.zeros((2,2))
			reim      = np.asarray([re.GetBinError(i+1, j+1),im.GetBinError(i+1, j+1)])
			coma      = np.multiply.outer(reim, reim)
			coma[0,1] = correl.GetBinContent(i+1, j+1)
			coma[1,0] = correl.GetBinContent(i+1, j+1)
			comaInvs.append(la.pinv(coma))
	return comaInvs

def getData(re,im,_):
	vals = []
	for i in xrange(re.GetNbinsX()):
		for j in xrange(im.GetNbinsY()):
			vals.append([re.GetBinContent(i+1, j+1), im.GetBinContent(i+1, j+1)])
	return np.asarray(vals)

def getChi2(vals1, vals2, comaInvs):
	c2    = 0.
	delta = vals1 - vals2
	for i in xrange(len(vals1)):
		c2 += np.dot(delta[i], np.dot(comaInvs[i], delta[i]))
	return c2

def getRangeMax(fileName):
	return float(fileName.split("_")[~0][8:~4])

def main():
	comaFileName = "./1mp1p1mmP_rangeMin0.0_rangeMax2.28.root"
	comaInvs     = getComas(*getHistos(comaFileName))
	fileList     = ["1mp1p1mmP_rangeMin0.0_rangeMax2.28.root","1mp1p1mmP_rangeMin0.0_rangeMax2.24.root","1mp1p1mmP_rangeMin0.0_rangeMax2.16.root"]
	vals1        = getData(*getHistos(fileList[0]))
	for i in xrange(len(fileList)-1):
		vals2 = getData(*getHistos(fileList[i+1]))
		print getChi2(vals1, vals2, comaInvs)

if __name__ == "__main__":
	main()
