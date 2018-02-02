import matplotlib
matplotlib.use("Agg")

import modernplotting.mpplot
import modernplotting.toolkit
import modernplotting.specialPlots as mpsp
import studyPlotter

import pyRootPwa
class comparisonPlot:
	def __init__(self, inFileName):
		first = True
		self.vals = []
		self.nY = 0
		with open(inFileName, 'r') as inin:
			for line in inin.readlines():
				if first:
					self.markers = line.split()
					self.nX = len(self.markers)
					first = False
					continue
				self.vals.append([float(v) for v in line.split()])
				if self.nY == 0:
					self.nY = len(self.vals[-1])
				else:
					if not len(self.vals[-1]) == self.nY:	
						raise ValueError("Dimension mismatch Y")
		if not self.nX == len(self.vals):
			raise ValueError("Dimension mismatch X")

	def writePDF(self, outFileName):
		hist = pyRootPwa.ROOT.TH2D("hist", "hist", self.nX, 0. ,self.nX, self.nY, 0.,self.nY)
		for i in range(self.nX):
			for j in range(self.nY):
				hist.SetBinContent(i+1, j+1, self.vals[i][j])

		style = self.getPlotStyle()
		with modernplotting.toolkit.PdfWriter(outFileName) as pdfOutput:
			plot = style.getPlot2D()
			plot.axes.get_xaxis().set_ticks([(i + 0.5) for i in range(self.nX)])
			plot.axes.get_yaxis().set_ticks([(i + 0.5) for i in range(self.nY)])
			studyPlotter.makeValuePlot(plot, hist)
			
			plot.axes.set_yticklabels(self.markers[:self.nY])
			plot.axes.set_xticklabels(self.markers, rotation = 90)
			plot.setZlim((0.,10.))

			pdfOutput.savefigAndClose()

	def getPlotStyle(self):
		style = modernplotting.mpplot.PlotterStyle()
		style.p2dColorMap = 'Reds'
		
		return style


def main():
	inFileName  = "./studies_0mp.txt"
	outFileName =  "./studies_0mp_newScript.pdf"
	lePlot      = comparisonPlot(inFileName)
	lePlot.writePDF(outFileName)

if __name__ == "__main__":
	main()
