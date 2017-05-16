import matplotlib
matplotlib.use("Agg")

from fabi import getch
import plotfabi as pf
import pyRootPwa
import numpy as np
import os, sys
from modes import INTENS, PHASE
from parseFiles import parseTH1D, parseTGraph, parseArgand
os.system(". add_to_front PATH /nfs/mnemosyne/sys/slc6/contrib/texlive/2013/bin/x86_64-linux")
import modernplotting.root
import modernplotting.mpplot
import modernplotting.colors
import modernplotting.toolkit
import modernplotting.specialPlots

import LaTeX_strings


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
	if mode == '01':
		return ["", "./results/0mp1mmIntens.pdf", ["./results/0mp1mm_only0pp.intens", "./results/0mp1mm_only1mm.intens"],
		            "./results/0mp1mmArgand.pdf", ["./results/0mp1mm_only0pp.argand", "./results/0mp1mm_only1mm.argand"]]
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

class resultViewer:
	def __init__(self, intensHists, realHists, imagHists, phaseHists, startBin = 34, startCommand = "", reImCorrel = None, noRun = False):
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
		self.theoColor = modernplotting.colors.makeColorLighter(modernplotting.colors.colorScheme.blue, .5)
		self.dataColor = modernplotting.colors.colorScheme.red
		self.addiColor = modernplotting.colors.colorScheme.gray

		self.lineArgs = {'marker': None, 'linestyle':'dashed', 'markersize':0, 'linewidth':.4, 'zorder' :0, 'color':'.5'}

		self.plotCorr  = True
		self.plotTheo  = True
		self.plotData  = True
		self.noEllipse = False

		self.mMin = 0.27
		self.mMax = 1.94

		self.noRun = noRun

		self.titleRight = ""
		self.tString  = ""

		self.tStringXpos = 0.015
		self.tStringYpos = 0.93

		if not startCommand == "":
			if startCommand.startswith("wq:"):
				fileName = startCommand[3:]
				self.writeAmplFiles(self.bin, fileName = fileName)
				self.noRun = True
			else:
				raise "Unknwon command '" + startCommand + "'"

		if not  self.noRun:
			self.intensCanvas = pyRootPwa.ROOT.TCanvas("Intensity","Intensity")
			self.sliceCanvas  = pyRootPwa.ROOT.TCanvas("IntensitySlice","IntensitySlice")
			self.argandCanvas = pyRootPwa.ROOT.TCanvas("Argand"        ,"Argand"        )
			self.phaseCanvas  = pyRootPwa.ROOT.TCanvas("Phase"        ,"Phase"        )
			self.intensCanvas.SetWindowSize(500,500)
			self.sliceCanvas.SetWindowSize( 500,500)
			self.argandCanvas.SetWindowSize(500,500)
			self.phaseCanvas.SetWindowSize( 500,500)


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
		return argand

	def getArgandData(self, nBin, index, getCOMA = False):
		if not self.reImCorrel and getCOMA:
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
			hist = self.intensHists[index].ProjectionY("_intens_"+str(index), nBin+1, nBin +1)
		elif mode == PHASE:
			hist = self.phaseHists[index].ProjectionY("_phase_"+str(index), nBin+1, nBin +1)
		if index > 0:
			hist.SetLineColor(index + 1)
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

		self.sliceCanvas.cd()
		self.sliceCanvas.Clear()
		slices = [self.getSlice(nBin)]
		slices[0].Draw()
		for i in range(1, self.nHists):
			slices.append(self.getSlice(nBin, i))
			slices[-1].Draw("SAME")
		self.sliceCanvas.Update()

		self.phaseCanvas.cd()
		self.phaseCanvas.Clear()
		phases = [self.getSlice(nBin, mode = PHASE)]
		phases[0].Draw()
		for i in range(1, self.nHists):
			slices.append(self.getSlice(nBin, index = i, mode = PHASE))
			slices[-1].Draw("SAME")
		self.phaseCanvas.Update()

		self.argandCanvas.cd()
		self.argandCanvas.Clear()
		argands = [self.getArgand(nBin)]
		argands[0].Draw()
		for i in range(1, self.nHists):
			argands.append(self.getArgand(nBin, i))
			argands[-1].Draw("SAME")
		self.argandCanvas.Update()

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
			elif cmd == 't':
				self.writeAmplFiles(self.bin)
			elif cmd == 'e':
				self.executeCommad()
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
		mMin = self.intensHists[0].GetXaxis().GetBinLowEdge(nBin+1)
		mMax = self.intensHists[0].GetXaxis().GetBinUpEdge(nBin+1)
		retVal = LaTeX_strings.getMassString(mMin, mMax)
		return retVal

	def writeBinToPdf(self, nBin, stdCmd = None):
		style = modernplotting.mpplot.PlotterStyle()
		if not self.titleRight == "":
			style.titleRight = self.titleRight

		if stdCmd:
			if not len(stdCmd) == 5:
				raise IndexError("Wrond number of standard commands")
		if stdCmd:
			twoDimPlotName = stdCmd[0]
		else:
			twoDimPlotName = raw_input("Name of the 2D plot:")
		if twoDimPlotName == "":
			print "No name given, skipping the 2D plot"
		elif not twoDimPlotName.endswith(".pdf"):
			print "Invalid file name extension in file name: '" + twoDimPlotName + "' skipping 2d plot"
		else:
			with modernplotting.toolkit.PdfWriter(twoDimPlotName) as pdfOutput:
				plot = style.getPlot2D()
				modernplotting.root.plotTH2D(self.intensHists[0], plot, maskValue = 0.)
				plot.setZlim((0., self.intensHists[0].GetMaximum()))
				plot.setXlabel(LaTeX_strings.m3Pi)
				plot.setYlabel(LaTeX_strings.m2Pi)
				plot.axes.text(self.tStringXpos,self.tStringYpos, self.tString, transform = plot.axes.transAxes)
				pdfOutput.savefigAndClose()
		if stdCmd:
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
				elif slicePlotName.endswith(".pdf") or slicePlotName == "":
					break
				else:
					print "Invalid file name given: '" + slicePlotName + "'"

		style.titleLeft = self.getLaTeXMassString(nBin)

		if slicePlotName == "":
			print "No name given, skipping the intensity slice"
		else:
			with modernplotting.toolkit.PdfWriter(slicePlotName) as pdfOutput:
				plot = style.getPlot1D()
				hists = [self.getSlice(nBin, index) for index in range(self.nHists)]
				for fn in addFiles:
					hist = parseTH1D(fn)
					modernplotting.root.plotTH1D(hist, plot, yErrors = True, maskValue = 0., markerDefinitions = { 'zorder':0, 'color':self.addiColor})
				if self.plotData:
					modernplotting.root.plotTH1D(hists[1], plot, yErrors = True, maskValue = 0., markerDefinitions = { 'zorder':1, 'color':self.dataColor})
				if len(hists) > 2 and self.plotTheo:
					modernplotting.root.plotTH1D(hists[2], plot, markerDefinitions = { 'zorder':2, 'color': self.theoColor})
				if self.plotCorr:
					modernplotting.root.plotTH1D(hists[0], plot, yErrors = True, maskValue = 0., markerDefinitions = { 'zorder':3, 'color':self.corrColor})
				plot.setXlabel(LaTeX_strings.m3Pi)
				plot.setXlim(self.mMin,self.mMax)
				plot.setYlabel(LaTeX_strings.intens)
				plot.axes.yaxis.offsetText.set_position((-.15,1.))
				plot.setYlim(0.,hists[0].GetMaximum()*1.2)
				plot.axes.text(self.tStringXpos,self.tStringYpos, self.tString, transform = plot.axes.transAxes)

				pdfOutput.savefigAndClose()
		if stdCmd:
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
				elif argandPlotName.endswith(".pdf") or argandPlotName == "":
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
						graph = parseTGraph(fn)
						addGraphs.append(graph)
						modernplotting.root.plotTH1D(graph, plot, yErrors = True, xErrors = True, maskValue = 0., markerDefinitions = {'linestyle' : 'solid', 'linewidth' : .2, 'zorder' : 0, 'color':self.addiColor})
					if self.plotData:
						modernplotting.root.plotTH1D(argands[1], plot, yErrors = True, xErrors = True, maskValue = 0., markerDefinitions = {'linestyle' : 'solid', 'linewidth' : .2, 'zorder' : 1, 'color':self.dataColor})
					if len(argands) > 2 and self.plotTheo:
						modernplotting.root.plotTH1D(argands[2], plot, maskValue = 0.,markerDefinitions = {'marker' : None, 'linestyle' : 'solid', 'zorder' : 2, 'color': self.theoColor})
					if self.plotCorr:
						modernplotting.root.plotTH1D(argands[0], plot, yErrors = True, xErrors = True, maskValue = 0., markerDefinitions = {'linestyle' : 'solid', 'linewidth' : .2, 'zorder' : 3, 'color' : self.corrColor})
					ranges = setAxesRange(argands[0])
	#				ranges = getMaximumRanges(argands + addGraphs)
	#				ranges = ((-500., 500.), (-500., 500.))
					plot.setXlabel(LaTeX_strings.real)
					plot.setYlabel(LaTeX_strings.imag)
					plot.axes.text(self.tStringXpos,self.tStringYpos, self.tString, transform = plot.axes.transAxes)
					plot.setXlim(ranges[0])
					plot.setYlim(ranges[1])
				else:
					ranges = setAxesRange(self.getArgand(nBin, 0))
					for fn in addFiles:
						X,EX,Y,EY = parseArgand(fn, skipZero = True)
						plot.plot(X, Y, **{'linestyle' : 'solid', 'linewidth' : .7, 'zorder' : 2, 'color':self.addiColor})
					if self.plotData:
						X,Y,COMA = self.getArgandData(nBin, 1, getCOMA = False)
						plot.plot(X, Y, **{'linestyle' : 'solid', 'linewidth' : .7, 'zorder' : 2, 'color':self.dataColor})
					if self.nHists > 2 and self.plotTheo:
						X,Y,COMA = self.getArgandData(nBin, 2, getCOMA = False)
						plot.plot(X, Y, **{'marker' : None, 'markersize' : 0., 'linestyle' : 'solid', 'linewidth' : .7,  'zorder' : 3, 'color': self.theoColor})
					if self.plotCorr:
						X,Y,COMA = self.getArgandData(nBin, 0, getCOMA = True)
						modernplotting.specialPlots.plotErrorEllipses(plot, X, Y, COMA, markerDefinitions = {'linestyle' : 'solid', 'linewidth' : 1., 'zorder' : 1, 'color' : self.corrColor})
						plot.plot(X,Y, **{'linestyle' : 'solid', 'linewidth' : 1., 'zorder' : 4, 'color' : self.corrColor})
					plot.setXlabel(LaTeX_strings.real)
					plot.setYlabel(LaTeX_strings.imag)
					plot.axes.text(self.tStringXpos,self.tStringYpos, self.tString, transform = plot.axes.transAxes)
					plot.setXlim(ranges[0])
					plot.setYlim(ranges[1])

				plot.plot([ranges[0][0],ranges[0][1]],[0.,0.], **self.lineArgs)
				plot.plot([0.,0.],[ranges[1][0],ranges[1][1]], **self.lineArgs)

				pdfOutput.savefigAndClose()
		print "Writing to .pdf finished."

