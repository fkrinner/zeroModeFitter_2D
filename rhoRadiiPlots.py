import matplotlib
matplotlib.use("Agg")
import modernplotting.mpplot
import modernplotting.toolkit
import modernplotting.specialPlots as mpsp

import LaTeX_strings
from utils import checkLaTeX

from math import isnan

from globalDefinitions import Pr0

def parseRadiusFile(inFileName, Pr0 = Pr0):
	data  = []
	nVals = None
	with open(inFileName, 'r') as inFile:
		for line in inFile.readlines():
			chunks = line.split()
			if nVals is None:
				nVals = len(chunks)
			else:
				if not nVals == len(chunks):
					raise ValueError("Numbers of values do not match")
			if nVals == 0:
				raise ValueError("No values found")
			vals = [float(chunks[1]), Pr0/float(chunks[2]), Pr0/float(chunks[2])**2 * float(chunks[3])] # Masses, values and errors: r = Pr0/Pr -> Dr = Pr0/Pr^2 DPr
			data.append(vals)
	return data

def makeRadiusPlots(inFileNameBase, outFileName, LaTeXwaveName = "", isobLaTeX = "", mMin = 1.0, mMax = 2.5, minn = 0.0001, maxx = 1., maxErrVal = None, dolog = True):
	style       = modernplotting.mpplot.PlotterStyle()
	style.titleRight = LaTeXwaveName

	color1 = modernplotting.colors.colorScheme.blue
	color2 = modernplotting.colors.colorScheme.red
	color3 = modernplotting.colors.colorScheme.green
	color4 = modernplotting.colors.colorScheme.orange

	colors      = [color1,color2,color3,color4]
	with  modernplotting.toolkit.PdfWriter(outFileName) as pdfOutput:
		plot = style.getPlot1D()
		plot.setXlabel(LaTeX_strings.m3Pi)
		plot.setYlim(minn, maxx)
		plot.setXlim(mMin, mMax)
		plot.setYlabel(r"$r_{"+isobLaTeX+r"}\,[\text{fm}]$")
		if dolog:
			plot.setYlog()
			plot.setYlim(minn, maxx)
		for t in range(4):
			inFileName   = inFileNameBase.replace("<tBin>", str(t))
			tBinVals     = parseRadiusFile(inFileName)
			x = [tBinVals[j][0] for j in range(len(tBinVals))]
			y = [tBinVals[j][1] for j in range(len(tBinVals))]
			e = [tBinVals[j][2] for j in range(len(tBinVals))]
			if maxErrVal is not None:
				nX = []
				nY = []
				nE = []
				for iE,err in enumerate(e):
					if err < maxErrVal:
						nX.append(x[iE])
						nY.append(y[iE])
						nE.append(err)
				x = nX
				y = nY
				e = nE
			if dolog:
				nX = []
				nY = []
				nE = []
				for iY, Y in enumerate(y):
					if Y >= minn:
						nX.append(x[iY])
						nY.append(y[iY])
						nE.append(e[iY])
				x = nX
				y = nY
				e = nE
			print y
			plot.axes.errorbar(x,y, yerr = e, color = colors[t], linestyle =  " ")
			plot.plot(x,y, color = colors[t], markeredgecolor = 'none', linestyle =  " ", markersize = 2. )

		plot.finishAndSaveAndClose(pdfOutput)

def main():
	checkLaTeX()
	toDo = "pi1"
	LaTeXwaveName  = r"$1^{-+}1^+[\pi\pi]_{1^{--}}\pi$P"
	if toDo == "pi1":
		inFileNameBase = "./pi1Radii_exotic_<tBin>.dat"
		outFileName    = "./pi1RadiiFromExotic.pdf"
		isobLaTeX      = r"\pi_1"
	elif toDo == "rho":
		inFileNameBase = "./rhoRadii_exotic_<tBin>.dat"
		outFileName    = "./rhoRadiiFromExotic.pdf"
		isobLaTeX      = r"\rho"
	makeRadiusPlots(inFileNameBase, outFileName, LaTeXwaveName, isobLaTeX)

if __name__ == "__main__":
	main()
