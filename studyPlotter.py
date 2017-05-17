import matplotlib
matplotlib.use("Agg")

import modernplotting.root
import modernplotting.mpplot
import modernplotting.colors
import modernplotting.toolkit
import modernplotting.specialPlots

def makeStringArray(hist2D):
	retVal = []
	for i in range(hist2D.GetNbinsX()):
		line = []
		for j in range(hist2D.GetNbinsY()):
			ss = "${0:.2f}$".format(hist2D.GetBinContent(i+1,j+1))
			line.append(ss)
		retVal.append(line)	
	return retVal

def makeValuePlot(plot, hist2D):
	modernplotting.root.plotTH2D(hist2D, plot)#, maskValue = 0.)
	stringArray = makeStringArray(hist2D)
#	modernplotting.specialPlots.plotValueString(plot, stringArray, fontsize = 8, fontweight='bold', color = 'w')
	modernplotting.specialPlots.plotValueString(plot, stringArray, fontsize = 8)

def main():
	import pyRootPwa
	from random import random
	
	nBins = 6
	hist = pyRootPwa.ROOT.TH2D("hist", "hist", nBins, 0., 1., nBins, 0.,1.)
	for i in range(nBins):
		for j in range(nBins):
			hist.SetBinContent(i+1, j+1, random()*2-1.)
	style = modernplotting.mpplot.PlotterStyle()

	with modernplotting.toolkit.PdfWriter("testValuePlot.pdf") as pdfOutput:
		plot = style.getPlot2D()
		makeValuePlot(plot, hist)
		pdfOutput.savefigAndClose()
if __name__ == "__main__":
	main()
