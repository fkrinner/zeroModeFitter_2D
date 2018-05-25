import matplotlib
matplotlib.use("Agg")

from fabi import getch
import plotfabi as pf
import pyRootPwa
import numpy as np
import os, sys
from modes import INTENS, PHASE, REAL, IMAG
from parseFiles import parseTH1D, parseTGraph, parseArgand
os.system(". add_to_front PATH /nfs/mnemosyne/sys/slc6/contrib/texlive/2013/bin/x86_64-linux")
import modernplotting.root
import modernplotting.mpplot
import modernplotting.colors
import modernplotting.toolkit
import modernplotting.specialPlots

from rootfabi import root_open

import LaTeX_strings
from cmath import phase

def findMinBin(hist):
	for i in range(0,hist.GetNbinsX()):
		if not hist.GetBinContent(i+1,1) == 0.:
			maxZeroBin = i-1
			break
	while True:
		allZero = True
		for i in range(hist.GetNbinsY()):
			if not hist.GetBinContent(maxZeroBin + 1, i+1) == 0.:
				allZero = False
				break
		if allZero:
			return maxZeroBin + 1
		maxZeroBin -= 1
		if maxZeroBin == -1:
			return maxZeroBin + 1

def getMaximumRanges(argands):
	inf  =  float('inf')
	xMin =  inf
	xMax = -inf
	yMin =  inf
	yMax = -inf
	for argand in argands:
		xMin = min(xMin, argand.GetXaxis().GetXmin())
		xMax = max(xMax, argand.GetXaxis().GetXmax())
		yMin = min(yMin, argand.GetYaxis().GetXmin())
		yMax = max(yMax, argand.GetYaxis().GetXmax())
	dX = xMax - xMin
	dY = yMax - yMin
	dd = dX - dY
	if dd < 0.:
		yMin -= dd/2
		yMax += dd/2
	else:
		xMin += dd/2
		xMax -= dd/2
	return (xMin, xMax), (yMin, yMax)

def setAxesRange(graph):
	hist = graph.GetHistogram()
	xMin   = hist.GetXaxis().GetXmin()
	xMax   = hist.GetXaxis().GetXmax()
	xRange = xMax - xMin
	yMin   = hist.GetYaxis().GetXmin()
	yMax   = hist.GetYaxis().GetXmax()
	yRange = yMax - yMin
	if xRange >= yRange:
		delRange = (xRange - yRange)/2
		hist.GetYaxis().SetRangeUser(yMin - delRange, yMax + delRange)
		ranges = ((xMin, xMax),(yMin - delRange, yMax + delRange))
	else:
		delRange = (yRange - xRange)/2
		hist.GetXaxis().SetLimits(xMin - delRange, xMax + delRange)
		ranges = ((xMin - delRange, xMax + delRange), (yMin, yMax))
	return ranges

def mergeRanges(ranges1, ranges2):
	xMin   = min(ranges1[0][0], ranges2[0][0])
	xMax   = max(ranges1[0][1], ranges2[0][1])
	xRange = xMax - xMin
	yMin   = min(ranges1[1][0], ranges2[1][0])
	yMax   = max(ranges1[1][1], ranges2[1][1])
	yRange = yMax - yMin
	if xRange >= yRange:
		delRange = (xRange - yRange)/2
		ranges = ((xMin, xMax),(yMin - delRange, yMax + delRange))
	else:
		delRange = (yRange - xRange)/2
		ranges = ((xMin - delRange, xMax + delRange), (yMin, yMax))
	return ranges

def getStdCmd():
	print "00 : 0mp0pp"
	print "01 : 0mp1mm"
	print "------------"
	print "10 : 1pp0pp"
	print "11 : 1pp1mm"
	print "------------"
	print "20 : 2mp0pp"
	print "21P: 2mp1mmP"
	print "21F: 2mp1mmF"
	print "22 : 2mp2pp"
	print "------------"
	print "f  : Folder"
	mode = raw_input()
	if mode == '00':
		return ["", "./results/0mp0ppIntens.pdf", ["./results/0mp0pp_only0pp.intens", "./results/0mp0pp_only1mm.intens"],
		            "./results/0mp0ppArgand.pdf", ["./results/0mp0pp_only0pp.argand", "./results/0mp0pp_only1mm.argand"]]
	if mode == '00-':
		return ["", "./results/0mp0ppIntens.pdf", ["./results/0mp0pp_only0pp.intens"],
		            "./results/0mp0ppArgand.pdf", ["./results/0mp0pp_only0pp.argand"]]
	if mode == '01':
		return ["", "./results/0mp1mmIntens.pdf", ["./results/0mp1mm_only0pp.intens", "./results/0mp1mm_only1mm.intens"],
		            "./results/0mp1mmArgand.pdf", ["./results/0mp1mm_only0pp.argand", "./results/0mp1mm_only1mm.argand"]]
	if mode == '01-':
		return ["", "./results/0mp1mmIntens.pdf", ["./results/0mp1mm_only0pp.intens"],
		            "./results/0mp1mmArgand.pdf", ["./results/0mp1mm_only0pp.argand"]]
	if mode == '10':
		return ["", "./results/1pp0ppIntens.pdf", ["./results/1pp0pp_only0pp.intens", "./results/1pp0pp_only1mm.intens"],
		            "./results/1pp0ppArgand.pdf", ["./results/1pp0pp_only0pp.argand", "./results/1pp0pp_only1mm.argand"]]
	if mode == '11':
		return ["", "./results/1pp1mmIntens.pdf", ["./results/1pp1mm_only0pp.intens","./results/1pp1mm_only0pp.intens"],
		            "./results/1pp1mmArgand.pdf", ["./results/1pp1mm_only0pp.argand", "./results/1pp1mm_only1mm.argand"]]
	if mode == "20":
		return ["", "./results/2mp0ppIntens.pdf", ["./results/2mp0pp_only0pp.intens","./results/2mp0pp_only1mmP.intens","./results/2mp0pp_only1mmF.intens","./results/2mp0pp_only2pp.intens"],
		            "./results/2mp0ppArgand.pdf", ["./results/2mp0pp_only0pp.argand","./results/2mp0pp_only1mmP.argand","./results/2mp0pp_only1mmF.argand","./results/2mp0pp_only2pp.argand"]]
	if mode == "21P":
		return ["", "./results/2mp1mmPIntens.pdf", ["./results/2mp1mmP_only0pp.intens","./results/2mp1mmP_only1mmP.intens","./results/2mp1mmP_only1mmF.intens","./results/2mp1mmP_only2pp.intens"],
		            "./results/2mp1mmPArgand.pdf", ["./results/2mp1mmP_only0pp.argand","./results/2mp1mmP_only1mmP.argand","./results/2mp1mmP_only1mmF.argand","./results/2mp1mmP_only2pp.argand"]]
	if mode == "21F":
		return ["", "./results/2mp1mmFIntens.pdf", ["./results/2mp1mmF_only0pp.intens","./results/2mp1mmF_only1mmP.intens","./results/2mp1mmF_only1mmF.intens","./results/2mp1mmF_only2pp.intens"],
		            "./results/2mp1mmFArgand.pdf", ["./results/2mp1mmF_only0pp.argand","./results/2mp1mmF_only1mmP.argand","./results/2mp1mmF_only1mmF.argand","./results/2mp1mmF_only2pp.argand"]]
	if mode == "22":
		return ["", "./results/2mp2ppIntens.pdf", ["./results/2mp2pp_only0pp.intens","./results/2mp2pp_only1mmP.intens","./results/2mp2pp_only1mmF.intens","./results/2mp2pp_only2pp.intens"],
		            "./results/2mp2ppArgand.pdf", ["./results/2mp2pp_only0pp.argand","./results/2mp2pp_only1mmP.argand","./results/2mp2pp_only1mmF.argand","./results/2mp2pp_only2pp.argand"]]
	if mode == 'f':
		folder   = raw_input("Folder:")
		suffix   = raw_input("Suffix:")
		intenses = []
		argands  = []
		for fn in os.listdir(folder):
			ffn = folder + os.sep + fn
			if fn.endswith(suffix+".intens"):
				intenses.append(ffn)
			if fn.endswith(suffix+".argand"):
				argands.append(ffn)
		return ["", "intens_"+folder+suffix+".pdf", intenses, "argands_"+folder+suffix+".pdf", argands]

def scaleHist(hist, factor):
	for x in range(hist.GetNbinsX()):
		for y in range(hist.GetNbinsY()):
			hist.SetBinContent(x+1,y+1,hist.GetBinContent(x+1,y+1)*factor)
			hist.SetBinError(x+1,y+1,hist.GetBinError(x+1,y+1)*factor)

class resultViewer:
	def __init__(self, intensHists, realHists, imagHists, phaseHists, startBin = 34, startCommand = None, reImCorrel = None, noRun = False, showPlots = None):
		pyRootPwa.ROOT.gStyle.SetOptStat(0)
		self.nHists       = len(intensHists)
		if self.nHists < 1:
			raise ValueError("No histograms given")
		self.bin          = startBin
		self.binMin       = findMinBin(intensHists[0])
		self.binMax       = intensHists[0].GetNbinsX()
		self.binsY        = intensHists[0].GetNbinsY()
		self.intensHists  = intensHists
		self.realHists    = realHists
		self.imagHists    = imagHists
		self.phaseHists   = phaseHists
		self.reImCorrel   = reImCorrel
		if not len(self.realHists) ==self. nHists or not len(self.imagHists) == self.nHists or not len(self.phaseHists) == self.nHists:
			raise ValueError("Size of histograms does not match")

		self.corrColor = modernplotting.colors.colorScheme.blue
		self.theoColor = modernplotting.colors.makeColorLighter(modernplotting.colors.colorScheme.gray, .2)
		self.dataColor = modernplotting.colors.colorScheme.red
#		self.addiColor = modernplotting.colors.colorScheme.gray
		self.addiColor = modernplotting.colors.makeColorLighter(modernplotting.colors.colorScheme.gray, .2)

		self.lineArgs = {'marker': None, 'linestyle':'dashed', 'markersize':0, 'linewidth':.4, 'zorder' :0, 'color':'.5'}

		self.plotCorr     = True
		self.plotTheo     = True
		self.plotData     = True
		self.noEllipse    = False
		self.XcheckArgand = False

		self.addiXoffset = None

		self.mMin = 0.27
		self.mMax = 1.94

		self.intensLabel = LaTeX_strings.intens
		self.realLabel   = LaTeX_strings.real
		self.imagLabel   = LaTeX_strings.imag
		self.m2PiString  = LaTeX_strings.m2Pi
		self.m3PiString  = LaTeX_strings.m3Pi

		self.scaleArgandZeroLine = 1.
		self.titleFontSize       = 13
		self.showColorBar        = False

		self.makeLegend        = False
		self.legendCorrected   = r"Corrected zero mode"
		self.legendUncorrected = r"Uncorrected zero mode"
		self.legendFixed       = r"Fixed shape"
		self.legendMethods     = r"Single methods"
		self.twoDtitleLeft       = None

		self.noRun = noRun
		if showPlots is None:
			self.showPlots = [True, True, False, False, False]
		else:
			if not len(showPlots) == 5:
				raise IndexError("'showPlots' invalid")
			self.showPlots = showPlots[:]

		self.connectPoints      = []
		self.labelPoints        = []
		self.shiftMap           = {}
		self.yAxisShift         = 0.
		self.xticks             = None
		self.yticks             = None


		self.titleRight         = ""
		self.tString            = ""
		self.overrideMassString = ""

		self.tStringXpos = 0.017
		self.tStringYpos = 0.93

		if startCommand is not None:
			if startCommand.startswith("wq:"):
				fileName = startCommand[3:]
				self.writeAmplFiles(self.bin, fileName = fileName)
				self.noRun = True
			else:
				raise "Unknwon command '" + startCommand + "'"

		
		self.sliceCanvasses = [None,            None,    None,   None,  None  ]
		self.sliceNames     = ["IntensitySlice","Argand","Phase","Real","Imag"]
		self.sliceModes     = [ INTENS,          None,    PHASE,  REAL,  IMAG ]

		if not self.noRun:
			self.intensCanvas = pyRootPwa.ROOT.TCanvas("Intensity","Intensity")
			self.intensCanvas.SetWindowSize(500,500)
			self.initSliceCanvasses()

		self.scaleFakk       = 1.
		self.printLiminary   = False
		self.scaleTo         = "corr"
		self.topMarginIntens = 1.2

	def writeToRootFile(self, rootFileName):
		with root_open(rootFileName, "RECREATE") as outFile:
			for i,hist in enumerate(self.intensHists):
				hist.SetName("intens"+str(i))
				hist.Write()
			for i,hist in enumerate(self.realHists):
				hist.SetName("real"+str(i))
				hist.Write()
			for i,hist in enumerate(self.imagHists):
				hist.SetName("imag"+str(i))
				hist.Write()
			for i,hist in enumerate(self.phaseHists):
				hist.SetName("phase"+str(i))
				hist.Write()

	def replaceFromRootFile(self, rootFileName, index):
		with root_open(rootFileName, "READ") as inFile:
			histIntens = inFile.Get("intens"+str(index))
			if not histIntens:
				raise IOError("Could not load 'intens"+str(index)+ "'")
			histIntens.SetDirectory(0)
			self.intensHists[index] = histIntens

			histReal = inFile.Get("real"+str(index))
			if not histReal:
				raise IOError("Could not load 'real"+str(index)+ "'")
			histReal.SetDirectory(0)
			self.realHists[index] = histReal

			histImag = inFile.Get("imag"+str(index))
			if not histImag:
				raise IOError("Could not load 'imag"+str(index)+ "'")
			histImag.SetDirectory(0)
			self.imagHists[index] = histImag

			histPhase = inFile.Get("phase"+str(index))
			if not histPhase:
				raise IOError("Could not load 'phase"+str(index)+ "'")
			histPhase.SetDirectory(0)
			self.phaseHists[index] = histPhase

	def initSliceCanvasses(self):
		for i in range(5):
			if self.showPlots[i]:
				if not self.sliceCanvasses[i]:
					self.sliceCanvasses[i]  = pyRootPwa.ROOT.TCanvas(self.sliceNames[i],self.sliceNames[i])
					self.sliceCanvasses[i].SetWindowSize(500,500)
			else:
				if self.sliceCanvasses[i]:
					self.sliceCanvasses[i].Close()
					self.sliceCanvasses[i] = None	

	def setShowPlots(self):
		self.showPlots = []
		for i in range(5):
			while True:
				ans = raw_input("Show "+self.sliceNames[i]+" slice? (y/n)")
				if ans == "y":
					self.showPlots.append(True)
					break
				if ans == "n":
					self.showPlots.append(False)
					break
		self.initSliceCanvasses()

	def scale(self, factor):
		for hist in self.intensHists:
			scaleHist(hist, factor)
		for hist in self.realHists:
			scaleHist(hist, factor**.5)
		for hist in self.imagHists:
			scaleHist(hist, factor**.5)
		scaleHist(self.reImCorrel, factor)
		self.scaleFakk *= factor

	def getArgand(self, nBin, index = 0):
		if index >= self.nHists:
			raise IndexError("Bin index too large")
		reals = []
		imags = []
		reErr = []
		imErr = []
		for i in range(self.binsY):
			re  = self.realHists[index].GetBinContent(nBin + 1, i+1)
			im  = self.imagHists[index].GetBinContent(nBin + 1, i+1)
			reR = self.realHists[index].GetBinError(  nBin + 1, i+1)
			imR = self.imagHists[index].GetBinError(  nBin + 1, i+1)
			if re == 0. and im == 0.:
				continue
			reals.append(re)
			imags.append(im)
			reErr.append(reR)
			imErr.append(imR)
		argand = pyRootPwa.ROOT.TGraphErrors(len(reals), np.asarray(reals, dtype = np.float64), np.asarray(imags, dtype = np.float64), np.asarray(reErr, np.float64), np.asarray(imErr, dtype = np.float64))
		if index > 0:
			argand.SetLineColor(index + 1)
		else:
			setAxesRange(argand)
#		print "========================================================"
#		for i in range(len(reals)):
#			print index, i,'a',reals[i]**2 +imags[i]**2,reals[i],imags[i],phase(reals[i]+1.j*imags[i])
#		print "========================================================"
		return argand

	def getArgandData(self, nBin, index, getCOMA = False):
		if getCOMA and self.reImCorrel is None:
			raise RuntimeError("Can not get COMA with no reImCorrel-histogram set")
		reals = []
		imags = []
		comas = []
		for i in range(self.binsY):
			re  = self.realHists[index].GetBinContent(nBin + 1, i+1)
			im  = self.imagHists[index].GetBinContent(nBin + 1, i+1)
			if re == 0. and im == 0.:
				continue
			reals.append(re)
			imags.append(im)

			if getCOMA:
				reR    = self.realHists[index].GetBinError(  nBin + 1, i+1)
				imR    = self.imagHists[index].GetBinError(  nBin + 1, i+1)
				correl = self.reImCorrel.GetBinContent(      nBin + 1, i+1)
				coma = [ [reR**2, correl],[correl, imR**2] ]
				comas.append(coma)
		return reals, imags, comas

	def getSlice(self, nBin, index = 0, mode = INTENS):
		if index >= self.nHists:
			raise IndexError("Bin index too large")
		if mode == INTENS:
#			b = 'I'
			hist = self.intensHists[index].ProjectionY("intens_"+str(index), nBin+1, nBin +1)
		elif mode == REAL:
			b = 'r'
			hist = self.realHists[index].ProjectionY("real_"+str(index), nBin+1, nBin +1)
		elif mode == IMAG:			
			b = 'i'
			hist = self.imagHists[index].ProjectionY("imag_"+str(index), nBin+1, nBin +1)
		elif mode == PHASE:
#			b = 'p'
			hist = self.phaseHists[index].ProjectionY("phase_"+str(index), nBin+1, nBin +1)
		if index > 0:
			hist.SetLineColor(index + 1)
#		print "========================================================"
#		for i in range(hist.GetNbinsX()):
#			if not hist.GetBinContent(i+1) == 0.:
#				print index,i,b, hist.GetBinContent(i+1)
#		print "========================================================"
		hist.SetTitle(hist.GetName())
		return hist

	def getMarkerLines(self, nBin):
		yMin = self.intensHists[0].GetYaxis().GetXmin()
		yMax = self.intensHists[0].GetYaxis().GetXmax()
		xMin = self.intensHists[0].GetXaxis().GetBinLowEdge(nBin+1)
		xMax = self.intensHists[0].GetXaxis().GetBinUpEdge( nBin+1)
		minLine = pyRootPwa.ROOT.TLine(xMin, yMin, xMin, yMax)
		maxLine = pyRootPwa.ROOT.TLine(xMax, yMin, xMax, yMax)
		minLine.SetLineWidth(3)
		maxLine.SetLineWidth(3)
		return maxLine, minLine

	def drawBin(self, nBin):
		self.intensCanvas.cd()
		self.intensCanvas.Clear()
		maxLine, minLine = self.getMarkerLines(nBin)
		self.intensHists[0].Draw("COL")
		maxLine.Draw()
		minLine.Draw()
		self.intensCanvas.Update()

		for i in range(5):
			if self.showPlots[i] and not i == 1:
				self.sliceCanvasses[i].cd()
				self.sliceCanvasses[i].Clear()
				slices = [self.getSlice(nBin, mode = self.sliceModes[i])]
				slices[0].Draw()
				for h in range(1, self.nHists):
					slices.append(self.getSlice(nBin, h, mode = self.sliceModes[i]))
					slices[-1].Draw("SAME")
				self.sliceCanvasses[i].Update()
			elif i == 1 and self.showPlots[i]:	
				self.sliceCanvasses[i].cd()
				self.sliceCanvasses[i].Clear()
				argands = [self.getArgand(nBin)]
				argands[0].Draw()
				for h in range(1, self.nHists):
					argands.append(self.getArgand(nBin, h))
					argands[-1].Draw("SAME")
				self.sliceCanvasses[i].Update()

	def run(self):
		while True:
			if self.noRun:
				break
			self.drawBin(self.bin)
			cmd = getch()
			if cmd == 'q':
				break
			if cmd == "Q":
				sys.exit(0)
			elif cmd == 's':
				self.bin += 1
				if self.bin == self.binMax:
					self.bin -= 1
			elif cmd == 'a':
				self.bin -= 1
				if self.bin < self.binMin:
					self.bin += 1
			elif cmd == 'p':
				self.intensHists[0].SetMaximum(self.intensHists[0].GetMaximum()/1.2)
			elif cmd == 'm':
				self.intensHists[0].SetMaximum(self.intensHists[0].GetMaximum()*1.2)
			elif cmd == 'w':
				self.writeBinToPdf(self.bin)
			elif cmd == 'r':
				self.wiriteReImToPdf(self.bin)
			elif cmd == 't':
				self.writeAmplFiles(self.bin)
			elif cmd == 'e':
				self.executeCommad()
			elif cmd == 'l':
				self.setShowPlots()
			elif cmd == '*':
				stdCmd = getStdCmd()
				self.writeBinToPdf(self.bin, stdCmd)
			else:
				print "Unknown command '" + cmd + "'"

	def executeCommad(self):
		commad = raw_input(">>> ")
		try:
			exec(commad)
		except:
			print "Could not execute '" + command + "'"

	def writeAmplFiles(self, nBin, index = 0, fileName = ""):
		if fileName == "":
			name = raw_input("outputFileName:")
		else:
			name = fileName
		if name == "":
			print "no name given"
			return
		with open(name + ".intens", 'w') as out:
			hist = self.getSlice(nBin, index)
			for i in range(hist.GetNbinsX()):
				out.write(str(hist.GetXaxis().GetBinLowEdge(i+1)) + ' ' + str(hist.GetXaxis().GetBinUpEdge(i+1)) + ' ' + str(hist.GetBinContent(i+1)) + ' ' + str(hist.GetBinError(i+1)) + '\n')
		with open(name + ".argand", 'w') as out:
			argand = self.getArgand(nBin, index)
			for i in range(argand.GetN()):
				X  = argand.GetX()[i]
				Y  = argand.GetY()[i]
				XE = argand.GetErrorX(i)
				YE = argand.GetErrorY(i)
				out.write(str(X) + ' ' + str(XE) + ' ' + str(Y) + ' ' + str(YE) + '\n')
#		print "Written."

	def getLaTeXMassString(self, nBin):
		if not self.overrideMassString == "":
			return self.overrideMassString
		mMin = self.intensHists[0].GetXaxis().GetBinLowEdge(nBin+1)
		mMax = self.intensHists[0].GetXaxis().GetBinUpEdge( nBin+1)
		retVal = LaTeX_strings.getMassString(mMin, mMax)
		return retVal

	def writeBinToPdf(self, nBin, stdCmd = None):
		style = modernplotting.mpplot.PlotterStyle()

		xSize1 = 4.
		xSize2 = 5.
		ySize  = 4.

		style.p1dFigSize = (xSize1,ySize)
		style.p2dFigSize = (xSize2,ySize)

		style.setFixed1DAxisPos(0.19, 0.16, 0.77, 0.77)
		style.setFixed2DAxisPos(0.19*xSize1/xSize2, 0.16, 0.655, 0.77)

		style.setFontSize(self.titleFontSize)

#		if self.showColorBar:
#			style.p2dFigSize = np.array([style.p2dFigSize[0]*1.25, style.p2dFigSize[1]])
#			style.p2dAxisPos[0] /= 1.25

		style.errorEllipsesEdgeColor = modernplotting.colors.makeColorLighter(self.corrColor, .5)
		style.errorEllipsesFaceColor = modernplotting.colors.makeColorLighter(self.corrColor, .5)
		style.finishPreliminary      = self.printLiminary

		if self.twoDtitleLeft is not None:
			style.titleLeft	= self.twoDtitleLeft

		if not self.titleRight == "":
			style.titleRight = self.titleRight

		if stdCmd is not None:
			if not len(stdCmd) == 5:
				raise IndexError("Wrong number of standard commands")
		if stdCmd is not None:
			twoDimPlotName = stdCmd[0]
		else:
			twoDimPlotName = raw_input("Name of the 2D plot:")
		if twoDimPlotName == "":
			print "No name given, skipping the 2D plot"
		elif not twoDimPlotName.endswith(".pdf") and not twoDimPlotName.endswith(".png"):
			print "Invalid file name extension in file name: '" + twoDimPlotName + "' skipping 2d plot"
		else:
			with modernplotting.toolkit.PdfWriter(twoDimPlotName) as pdfOutput:
				plot = style.getPlot2D()
				modernplotting.root.plotTH2D(self.intensHists[0], plot, maskValue = 0.)
				if self.showColorBar:
					plot.setZshowColorBar()
					plot.colorbar.ax.yaxis.offsetText.set_position((4.,1.))
					plot.colorbar.formatter.set_powerlimits((2, 6))
				plot.setZlim((0., self.intensHists[0].GetMaximum()))
				plot.setXlabel(self.m3PiString)
				plot.setYlabel(self.m2PiString)
				plot.axes.text(self.tStringXpos,self.tStringYpos, self.tString, transform = plot.axes.transAxes, size = self.titleFontSize)

				plot.finishAndSaveAndClose(pdfOutput)
				
		if stdCmd is not None:
			slicePlotName = stdCmd[1]
			addFiles      = stdCmd[2]
		else:
			addFiles = []
			while True:
				slicePlotName = raw_input("Name of the intensity slice:")
				if slicePlotName.endswith(".intens"):
					fileName = slicePlotName
					if not os.path.isfile(fileName):
						print "file '" + fileName + "' does not exist"
					else:
						addFiles.append(fileName)
						print "Added '" + fileName + "' as additional plot"
				elif slicePlotName.endswith(".pdf") or slicePlotName.endswith(".png") or slicePlotName == "":
					break
				else:
					print "Invalid file name given: '" + slicePlotName + "'"

		style.titleLeft = self.getLaTeXMassString(nBin)

		handlers = []
		legends  = []
		if self.plotCorr:
			handlers.append(matplotlib.patches.Patch(color = self.corrColor))
			legends.append(self.legendCorrected)
		if self.plotData:
			handlers.append(matplotlib.patches.Patch(color = self.dataColor))
			legends.append(self.legendUncorrected)
		if self.plotTheo:
			handlers.append(matplotlib.patches.Patch(color = self.theoColor))
			legends.append(self.legendFixed)
		if len(addFiles) > 0:
			handlers.append(matplotlib.patches.Patch(color = self.addiColor))
			legends.append(self.legendMethods)

		pointLabels     = {}
		for i in self.labelPoints:
			mass = self.intensHists[0].GetYaxis().GetBinCenter(i+1)
			pointLabels[i] = "{:1.2f}".format(mass)
		if slicePlotName == "":
			print "No name given, skipping the intensity slice"
		else:
			with modernplotting.toolkit.PdfWriter(slicePlotName) as pdfOutput:
				plot = style.getPlot1D()
				hists = [self.getSlice(nBin, index) for index in range(self.nHists)]
				
				nnBins = hists[0].GetNbinsX()
				firstUpperBinBorder = -1.
				for i in range(nnBins):
					if not hists[0].GetBinContent(nnBins-i) == 0.:
						firstUpperBinBorder = hists[0].GetXaxis().GetBinUpEdge(nnBins-i)
						break
				if firstUpperBinBorder == -1.:
					raise ValueError("Could not determine upper limit")
				if self.plotData:
					modernplotting.root.plotTH1D(hists[1], plot, yErrors = True, maskValue = 0., markerDefinitions = { 'zorder':1, 'color':self.dataColor})
				if len(hists) > 2 and self.plotTheo:
					modernplotting.root.plotTH1D(hists[2], plot, markerDefinitions = { 'zorder':2, 'color': self.theoColor})
				if self.plotCorr:
					modernplotting.root.plotTH1D(hists[0], plot, yErrors = True, maskValue = 0., markerDefinitions = { 'zorder':3, 'color':self.corrColor})
				for fn in addFiles:
					hist = parseTH1D(fn, self.scaleFakk, addX = self.addiXoffset)
					modernplotting.root.plotTH1D(hist, plot, yErrors = True, maskValue = 0., markerDefinitions = { 'zorder':0, 'color':self.addiColor})
				plot.setXlabel(self.m2PiString)
#				plot.setXlim(self.mMin,self.mMax)
				plot.setXlim(self.mMin,firstUpperBinBorder)
				plot.setYlabel(self.intensLabel)
				plot.axes.yaxis.offsetText.set_position((-.15,1.))
				if self.scaleTo == "corr":
					plot.setYlim(0.,hists[0].GetMaximum()*self.topMarginIntens)
				elif self.scaleTo == "maxCorrData":
					plot.setYlim(0.,max(hists[0].GetMaximum(),hists[1].GetMaximum())*self.topMarginIntens)
				else:
					raise RuntimeError("Unknwons scale option '" + self.scaleTo + "'")
				
				plot.axes.text(self.tStringXpos,self.tStringYpos, self.tString, transform = plot.axes.transAxes, size = self.titleFontSize)

				if self.makeLegend:
					plot.fig.legend(handlers, legends, fontsize = "x-small",loc =0 , mode = "expand", ncol = 2, borderaxespad=0., bbox_to_anchor=style.p1dAxisPos)	

				plot.finishAndSaveAndClose(pdfOutput)
		if stdCmd is not None:
			argandPlotName = stdCmd[3]
			addFiles       = stdCmd[4]
			for fn in addFiles:
				if not os.path.isfile(fn):
					raise IOError("File '" + fn + "' does not exist")
		else:
			addFiles = []
			while True:
				argandPlotName = raw_input("Name of the argand plot:")
				if argandPlotName.endswith(".argand"):
					fileName = argandPlotName
					if not os.path.isfile(fileName):
						print "file '" + fileName + "' does not exist"
					else:
						addFiles.append(fileName)
						print "Added '" + fileName + "' as additional plot"
				elif argandPlotName.endswith(".pdf") or argandPlotName.endswith(".png") or argandPlotName == "":
					break
				else:
					print "Invalid file name: '" + argandPlotName + "'"

		if argandPlotName == "":
			print "No name given, skipping the argand plot"
		else:
			with modernplotting.toolkit.PdfWriter(argandPlotName) as pdfOutput:
				plot = style.getPlot1D()
				if self.noEllipse or not self.reImCorrel:
					argands   = [self.getArgand(nBin, index) for index in range(self.nHists)]
					addGraphs = []
					for fn in addFiles:
						graph = parseTGraph(fn, self.scaleFacc**.5)
						addGraphs.append(graph)
						modernplotting.root.plotTH1D(graph,plot,yErrors=True,xErrors=True,maskValue=0.,markerDefinitions={'linestyle':'solid','linewidth':.2,'zorder':0,'color':self.addiColor})
					if self.plotData:
						modernplotting.root.plotTH1D(argands[1],plot,yErrors=True,xErrors=True,maskValue=0.,markerDefinitions={'linestyle':'solid','linewidth':.2,'zorder':2,'color':self.dataColor})
					if len(argands) > 2 and self.plotTheo:
						modernplotting.root.plotTH1D(argands[2],plot,maskValue=0.,markerDefinitions={'marker':None,'linestyle':'solid','zorder':1,'color':self.theoColor,'linewidth':1.})
					if self.plotCorr:
						modernplotting.root.plotTH1D(argands[0],plot,yErrors=True,xErrors=True,maskValue=0.,markerDefinitions={'marke':None,'linestyle':'solid','linewidth':.2,'zorder':3,'color':self.corrColor})
					ranges   = setAxesRange(argands[0])
				else:
					ranges   = setAxesRange(self.getArgand(nBin, 0))
					for fn in addFiles:
						X,EX,Y,EY = parseArgand(fn, skipZero = True, fakk = self.scaleFakk**.5)
						hasAdd = True
						if not self.XcheckArgand:
							plot.plot(X,Y,**{'linestyle':'solid','linewidth':.7,'zorder':2,'color':self.addiColor})
						else:
							plot.plot(X,Y,**{'linestyle':'solid','linewidth':.0,'zorder':5,'color':'red','markersize':2.,'markeredgecolor':'none'}) # here here ololo
					if self.plotData:
						X,Y,COMA = self.getArgandData(nBin, 1, getCOMA = False)
#						print "data"
#						print X
#						print Y
#						print "/data"
						plot.plot(X, Y,**{'linestyle':'solid','linewidth':.7,'zorder':3,'color':self.dataColor})
					if self.nHists > 2 and self.plotTheo:
						X,Y,COMA = self.getArgandData(nBin,2,getCOMA=False)
						plot.plot(X, Y, **{'marker':None,'markersize':0.,'linestyle':'solid','linewidth':1.,'zorder':2,'color':self.theoColor})
						if len(self.connectPoints) > 0:
							XD, YD, CD = self.getArgandData(nBin, 0, getCOMA = True)
							for p in self.connectPoints:
								if len(X) <= p:
									print "WARNING: Connect point",p,"not possible. Index out of range"
									continue
								plot.plot([X[p],XD[p]],[Y[p],YD[p]], markersize = 0., color = 'k', linewidth = 0.1 )
					if self.plotCorr:
						X,Y,COMA = self.getArgandData(nBin, 0, getCOMA = True)
#						with open("coma_from_the_rv_C.dat",'w') as outFile:
#							for i in range(len(X)):
#								outFile.write(str(COMA[i][0][0])+' '+str(COMA[i][0][1])+' '+str(COMA[i][1][1])+'\n')

						modernplotting.specialPlots.plotErrorEllipses(plot,X,Y,COMA,markerDefinitions={'linestyle':'solid','linewidth':1.,'zorder':1,'color':self.corrColor,'markersize':0.})
						plot.plot(X,Y, **{'linestyle' : 'solid', 'linewidth' : 1., 'zorder' : 4, 'color' : self.corrColor})
						for i in pointLabels:
							if i < len(X):
								x = X[i]
								y = Y[i]
								if i in self.shiftMap:
									x += self.shiftMap[i][0]
									y += self.shiftMap[i][1]
								plot.axes.text(x,y,pointLabels[i],size = self.titleFontSize)


				if self.scaleTo == "corr":
					pass
				elif self.scaleTo == "maxCorrData":
					ranges    = mergeRanges(ranges, setAxesRange(self.getArgand(nBin, 1)))
				else:
					raise RuntimeError("Unknwons scale option '" + self.scaleTo + "'")

				if self.makeLegend:
					plot.fig.legend(handlers, legends, fontsize = "x-small",loc = 0 , mode = "expand", ncol = 2, borderaxespad=0., bbox_to_anchor=style.p1dAxisPos)	
				plot.setXlabel(self.realLabel)
				plot.setYlabel(self.imagLabel)
				plot.axes.yaxis.offsetText.set_position((-.15,1.))
				plot.axes.xaxis.offsetText.set_position((1.05,0.))
				plot.axes.text(self.tStringXpos,self.tStringYpos, self.tString, transform = plot.axes.transAxes,size = self.titleFontSize)
				fakk   = .1
				yRange = (ranges[1][0], ranges[1][1] +  fakk*(ranges[1][1] - ranges[1][0]))

				if self.xticks:
					plot.axes.xaxis.set_ticks(self.xticks)

				if self.yticks:
					plot.axes.yaxis.set_ticks(self.yticks)


				plot.setXlim(ranges[0])
				yRange = (yRange[0] + self.yAxisShift,yRange[1] + self.yAxisShift)
				plot.setYlim(yRange)

				plot.plot([ranges[0][0],ranges[0][1]],[0.,0.], **self.lineArgs)
				plot.plot([0.,0.],[yRange[0],yRange[0] + (yRange[1]-yRange[0])*self.scaleArgandZeroLine],**self.lineArgs)

				plot.finishAndSaveAndClose(pdfOutput)
		print "Finished with creating .pdf-files"

	def wiriteReImToPdf(self, nBin, outFileName = None):
		style = modernplotting.mpplot.PlotterStyle()
		xSize1 = 4.
		xSize2 = 5.
		ySize  = 4.
		style.setFixed1DAxisPos(0.19, 0.16, 0.77, 0.77)
		style.setFixed2DAxisPos(0.19*xSize2/xSize1, 0.16, 0.77, 0.77)
		style.finishPreliminary = self.printLiminary
		style.titleLeft = self.getLaTeXMassString(nBin)
		if not self.titleRight == "":
			style.titleRight = self.titleRight

		if outFileName is None:
			while True:
				outFileNameBase = raw_input("OutFileName (must contain <ri>):")
				if not outFileNameBase.endswith(".pdf"):
					print "File name must end with '.pdf'"
					continue
				if not "<ri>" in outFileNameBase:
					print "File name must contain '<ri>'"
					continue
				break
		else:
			outFileNameBase = outFileName

		if outFileNameBase == "":
			print "No name given, skipping"
			return False
		addFiles = []

		handlers = []
		legends  = []
		if self.plotCorr:
			handlers.append(matplotlib.patches.Patch(color = self.corrColor))
			legends.append(self.legendCorrected)
		if self.plotData:
			handlers.append(matplotlib.patches.Patch(color = self.dataColor))
			legends.append(self.legendUncorrected)
		if self.plotTheo:
			handlers.append(matplotlib.patches.Patch(color = self.theoColor))
			legends.append(self.legendFixed)
		if len(addFiles) > 0:
			handlers.append(matplotlib.patches.Patch(color = self.addiColor))
			legends.append(self.legendMethods)

		for mode in [(REAL,'real'),(IMAG,'imag')]:
			with modernplotting.toolkit.PdfWriter(outFileNameBase.replace("<ri>",mode[1])) as pdfOutput:
				plot = style.getPlot1D()
				hists = [self.getSlice(nBin, index, mode = mode[0]) for index in range(self.nHists)]
				nnBins = hists[0].GetNbinsX()
				firstUpperBinBorder = -1.
				for i in range(nnBins):
					if not hists[0].GetBinContent(nnBins-i) == 0.:
						firstUpperBinBorder = hists[0].GetXaxis().GetBinUpEdge(nnBins-i)
						break
				if firstUpperBinBorder == -1.:
					raise ValueError("Could not determine upper limit")
				for fn in addFiles:
					hist = parseTH1D(fn, addX = self.addiXoffset)
					modernplotting.root.plotTH1D(hist, plot, yErrors = True, maskValue = 0., markerDefinitions = { 'zorder':0, 'color':self.addiColor})
				if self.plotData:
					modernplotting.root.plotTH1D(hists[1], plot, yErrors = True, maskValue = 0., markerDefinitions = { 'zorder':1, 'color':self.dataColor})
				if len(hists) > 2 and self.plotTheo:
					modernplotting.root.plotTH1D(hists[2], plot, markerDefinitions = { 'zorder':2, 'color': self.theoColor})
				if self.plotCorr:
					modernplotting.root.plotTH1D(hists[0], plot, yErrors = True, maskValue = 0., markerDefinitions = { 'zorder':3, 'color':self.corrColor})
				plot.setXlabel(self.m2PiString)
#				plot.setXlim(self.mMin,self.mMax)
				plot.setXlim(self.mMin,firstUpperBinBorder)
				if mode[0] == REAL:
					plot.setYlabel(self.realLabel)
				elif mode[0] == IMAG:
					plot.setYlabel(self.imagLabel)
				plot.axes.yaxis.offsetText.set_position((-.15,1.))
				if self.scaleTo == "corr":
					maxx = hists[0].GetMaximum()
					minn = hists[0].GetMinimum()
					rang = maxx - minn

					plot.setYlim(minn,minn + rang*self.topMarginIntens)
				elif self.scaleTo == "maxCorrData":
					maxx = max(hists[0].GetMaximum(),hists[1].GetMaximum())
					minn = min(hists[0].GetMinimum(),hists[1].GetMinimum())
					rang = maxx - minn
					plot.setYlim(minn,minn + rang*self.topMarginIntens)
				else:
					raise RuntimeError("Unknwons scale option '" + self.scaleTo + "'")
			
				plot.axes.text(self.tStringXpos,self.tStringYpos, self.tString, transform = plot.axes.transAxes,size = self.titleFontSize)

				if self.makeLegend:
					plot.fig.legend(handlers, legends, fontsize = "x-small",loc =0 , mode = "expand", ncol = 2, borderaxespad=0., bbox_to_anchor=style.p1dAxisPos)	

				plot.finishAndSaveAndClose(pdfOutput)

		print "Writing to .pdf finished."
		return True

