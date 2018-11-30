#!/usr/bin/python
# makePhasePlots.py
# Created: 2018-10-23 19:19:18.969641
# Author: Fabian Krinner
import matplotlib
matplotlib.use("Agg")

import os, sys
from cmath import phase, pi

import ROOT

os.system(". add_to_front PATH /nfs/mnemosyne/sys/slc6/contrib/texlive/2013/bin/x86_64-linux")
import modernplotting.root
import modernplotting.mpplot
import modernplotting.colors
import modernplotting.toolkit
import modernplotting.specialPlots

from LaTeX_strings import tBins

def parseFile(inFileName):
	retVal = []
	with open(inFileName, 'r') as inFile:
		for line in inFile.readlines():
			vals = [float(v) for v in line.split()]
			retVal.append(vals)
	return retVal

def getPhase(vals, delta = 1.e-5):
	coma    = [[vals[5], vals[7]],[vals[7],vals[6]]]
	ampl    =  vals[3] + 1.j*vals[4]
	pha     = phase(ampl)
	jacPlus = [(phase(ampl+delta)-pha)/delta, (phase(ampl+1.j*delta)-pha)/delta]
	jacMinu = [(pha-phase(ampl-delta))/delta, (pha-phase(ampl-1.j*delta))/delta]
	jac     = [(jacPlus[0]+jacMinu[0])/2,(jacPlus[1]+jacMinu[1])/2]

	err = 0.
	for i in range(2):
		for j in range(2):
			err += jac[i]*jac[j]*coma[i][j]
	try:
		err **= .5
	except ValueError:
		err = 3.
	return pha, err

def getPhases(inFileName):
	dat = parseFile(inFileName)
	retVal = []
	for vals in dat:
		try:
			pha,err = getPhase(vals)
		except ValueError:
			print inFileName
			raise ValueError
		retVal.append([vals[0], pha,err])
	return retVal

def getPhaseHist(name, inFileName, leastDiff = True):
	phases = getPhases(inFileName)
	if leastDiff:
		for i in range(len(phases)-1):
			if phases[i+1][1] > phases[i][1] + pi:
				phases[i+1][1] -= 2*pi
			if phases[i+1][1] < phases[i][1] - pi:
				phases[i+1][1] += 2*pi
				
	hist   = ROOT.TH1D(name, name, 50, .5, 2.5)
	for vals in phases:
		nBin = hist.GetXaxis().FindBin(vals[0])
		hist.SetBinContent(nBin, vals[1])
		hist.SetBinError(nBin, vals[2])
	return hist

def getRangeRangeList(folderName):
	retVal = []
	for fn in os.listdir(folderName):
		if fn.startswith("1mp_rho_cpls_"):
			tBin     = int(fn[~4])
			chunks   = fn.split("_")
			rangeMin = float(chunks[~2][8:])
			rangeMax = float(chunks[~1][8:])
			retVal.append((folderName + os.sep + fn, tBin, rangeMin, rangeMax))
	return retVal

def main():
	style                  = modernplotting.mpplot.PlotterStyle()
	style.p1dDefaultMarker = '.'
	tStringXpos            = 0.017
	tStringYpos            = 0.93
	titleFontSize          = 13

	title = r"arg$(1^{-+}1^+\rho(770)\pi$P$)-$arg$(4^{-+}0^+\rho(770)\pi$F$)$"

	folderName = "./cplFiles"
	fileRangeList = getRangeRangeList(folderName)

	ofns = []
	for rangeFn in fileRangeList:
		tBin       = rangeFn[1]
		tString    = tBins[tBin]
		style.titleRight = title

		outFileName = "./cplPlots/tBin"+str(tBin)+"_rangeMin"+str(rangeFn[2])+"_rangeMax"+str(rangeFn[3])+".pdf"

		if not outFileName in ofns:
			ofns.append(outFileName)
		else:
			raise IOError("File '" + outFileName + "' found twice")

		inFileName = rangeFn[0]
		hist       = getPhaseHist("tBin"+str(tBin), inFileName)
#		if tBin == 3:
#			for i in range(50):
#				if hist.GetBinContent(i+1) > 4:
#					hist.SetBinContent(i+1, hist.GetBinContent(i+1) - 2*pi)
		with modernplotting.toolkit.PdfWriter(outFileName) as pdfOutput:
			style.titleRight = "$" + str(rangeFn[2]) + r" < m_{\pi^+\pi^-} < " + str(rangeFn[3]) + r"$\,GeV"
			plot = style.getPlot1D()
			modernplotting.root.plotTH1D(hist, plot, yErrors = True, maskValue = 0.)
			plot.setXlabel(r"$m_{3\pi}\,[\SI{}{GeV/c^2}]$")
			plot.setYlabel(r"Phase [rad]")
			plot.axes.text(tStringXpos,tStringYpos, tString, transform = plot.axes.transAxes, size = titleFontSize)
#			if tBin in [0,2,3]:
#				plot.setYlim(hist.GetMinimum(), hist.GetMaximum() + .6)
			print "I am here"
			plot.setYlim(2.,4.)
			plot.finishAndSaveAndClose(pdfOutput)

if __name__ == "__main__":
	main()
