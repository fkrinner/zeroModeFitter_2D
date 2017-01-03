from fabi import getch
import plotfabi as pf
import pyRootPwa
import numpy as np
import os
from parseFiles import parseTH1D, parseTGraph
os.system(". add_to_front PATH /nfs/mnemosyne/sys/slc6/contrib/texlive/2013/bin/x86_64-linux")
import modernplotting.root
import modernplotting.mpplot
import modernplotting.colors
import modernplotting.toolkit

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

class resultViewer:
	def __init__(self, intensHists, realHists, imagHists, startBin = 30):
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
		if not len(self.realHists) ==self. nHists or not len(self.imagHists) == self.nHists:
			raise ValueError("Size of histograms does not match")
		self.intensCanvas = pyRootPwa.ROOT.TCanvas("Intensity","Intensity")
		self.sliceCanvas  = pyRootPwa.ROOT.TCanvas("IntensitySlice","IntensitySlice")
		self.argandCanvas = pyRootPwa.ROOT.TCanvas("Argand"        ,"Argand"        )
		self.intensCanvas.SetWindowSize(500,500)
		self.sliceCanvas.SetWindowSize( 500,500)
		self.argandCanvas.SetWindowSize(500,500)

		self.corrColor = modernplotting.colors.colorScheme.blue
		self.theoColor = modernplotting.colors.makeColorLighter(modernplotting.colors.colorScheme.blue, .5)
		self.dataColor = modernplotting.colors.colorScheme.gray
		self.addiColor = modernplotting.colors.colorScheme.red

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
				break
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

	def getSlice(self, nBin, index = 0):
		if index >= self.nHists:
			raise IndexError("Bin index too large")
		intensHist = self.intensHists[index].ProjectionY("_py_"+str(index), nBin+1, nBin +1)
		if index > 0:
			intensHist.SetLineColor(index + 1)
		return intensHist

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
			self.drawBin(self.bin)
			cmd = getch()
			if cmd == 'q':
				break
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
			else:
				print "Unknown command '" + cmd + "'"

	def writeAmplFiles(self, nBin, index = 0):
		name = raw_input("outputFileName:")
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

	def writeBinToPdf(self, nBin):
		style = modernplotting.mpplot.PlotterStyle()
		twoDimPlotName = raw_input("Name of the 2D plot:")
		if twoDimPlotName == "":
			print "No name given, skipping the 2D plot"
		else:
			with modernplotting.toolkit.PdfWriter(twoDimPlotName) as pdfOutput:
				plot = style.getPlot2D()
				modernplotting.root.plotTH2D(self.intensHists[0], plot, maskValue = 0.)
				plot.setZlim((0., self.intensHists[0].GetMaximum()))
				plot.setXlabel(pf.mPiPiPi + ' ' + pf.MeVc2 )
				plot.setYlabel(pf.mPiPi   + ' ' + pf.MeVc2)
				pdfOutput.savefigAndClose()
		addFiles = []
		while True:
			slicePlotName = raw_input("Name of the intensity slice:")
			if slicePlotName.startswith(":"):
				fileName = slicePlotName[1:]
				if not os.path.isfile(fileName):
					print "file '" + fileName + "' does not exist"
				else:
					addFiles.append(fileName)
					print "Added '" + fileName + "' as additional plot"
			else:
				break
		if slicePlotName == "":
			print "No name given, skipping the intensity slice"
		else:
			with modernplotting.toolkit.PdfWriter(slicePlotName) as pdfOutput:
				plot = style.getPlot1D()
				hists = [self.getSlice(nBin, index) for index in range(self.nHists)]
				for fn in addFiles:
					hist = parseTH1D(fn)
					modernplotting.root.plotTH1D(hist, plot, yErrors = True, maskValue = 0., markerDefinitions = { 'zorder':0, 'color':self.addiColor})
				modernplotting.root.plotTH1D(hists[1], plot, yErrors = True, maskValue = 0., markerDefinitions = { 'zorder':1, 'color':self.dataColor})
				if len(hists) > 2:
					modernplotting.root.plotTH1D(hists[2], plot, markerDefinitions = { 'zorder':2, 'color': self.theoColor})
				modernplotting.root.plotTH1D(hists[0], plot, yErrors = True, maskValue = 0., markerDefinitions = { 'zorder':3, 'color':self.corrColor})
				plot.setXlabel(pf.mPiPi + ' ' + pf.MeVc2)
				plot.setYlabel(pf.intens + ' ' + pf.AU)
				plot.setYlim(0.,hists[0].GetMaximum()*1.2)
				pdfOutput.savefigAndClose()

		addFiles = []
		while True:
			argandPlotName = raw_input("Name of the argand plot:")
			if argandPlotName.startswith(":"):
				fileName = argandPlotName[1:]
				if not os.path.isfile(fileName):
					print "file '" + fileName + "' does not exist"
				else:
					addFiles.append(fileName)
					print "Added '" + fileName + "' as additional plot"
			else:
				break

		if argandPlotName == "":
			print "No name given, skipping the argand plot"
		else:
			with modernplotting.toolkit.PdfWriter(argandPlotName) as pdfOutput:
				plot = style.getPlot1D()
				argands = [self.getArgand(nBin, index) for index in range(self.nHists)]
				for fn in addFiles:
					graph = parseTGraph(fn)
					modernplotting.root.plotTH1D(graph, plot, yErrors = True, xErrors = True, maskValue = 0., markerDefinitions = {'linestyle' : 'solid', 'linewidth' : .2, 'zorder' : 0, 'color':self.addiColor})
				modernplotting.root.plotTH1D(argands[1], plot, yErrors = True, xErrors = True, maskValue = 0., markerDefinitions = {'linestyle' : 'solid', 'linewidth' : .2, 'zorder' : 1, 'color':self.dataColor})
				if len(argands) > 2:
					modernplotting.root.plotTH1D(argands[2], plot, maskValue = 0.,markerDefinitions = {'marker' : None, 'linestyle' : 'solid', 'zorder' : 2, 'color': self.theoColor})
				modernplotting.root.plotTH1D(argands[0], plot, yErrors = True, xErrors = True, maskValue = 0., markerDefinitions = {'linestyle' : 'solid', 'linewidth' : .2, 'zorder' : 3, 'color' : self.corrColor})
				ranges = setAxesRange(argands[0])
				plot.setXlabel(pf.realPart + ' ' + pf.AU)
				plot.setYlabel(pf.imagPart + ' ' + pf.AU)
				plot.setXlim(ranges[0])
				plot.setYlim(ranges[1])
				pdfOutput.savefigAndClose()

