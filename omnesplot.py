#!/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/Python_ultra/Python-2.7.10/bin/python
# omnesplot.py
# Created: 2018-07-09 16:14:34.144663
# Author: Fabian Krinner
import matplotlib
matplotlib.use("Agg")
import os, sys
os.system(". add_to_front PATH /nfs/mnemosyne/sys/slc6/contrib/texlive/2013/bin/x86_64-linux")
import modernplotting.root
import modernplotting.mpplot
import modernplotting.colors
import modernplotting.toolkit
import modernplotting.specialPlots


def getData():
	inFileName = "./Omnes11_new.dat"
	data = []
	with open(inFileName, "r") as inFile:
		for line in inFile.readlines():
			vals = [float(v) for v in line.split()]
			data.append(vals)
	return data

def main():
	data        = getData()
	outFileName = "./omnes.pdf"	
	style = modernplotting.mpplot.PlotterStyle()
	xSize1 = 4.
	xSize2 = 4.
	ySize  = 4.

	style.p2dFigSize = (xSize2,ySize)
	style.setFixed1DAxisPos(0.16, 0.16, 0.77, 0.77)
	style.setFixed2DAxisPos(0.16*xSize2/xSize1, 0.16, 0.77, 0.77)

	style.titleRight = r"$[\pi\pi]_{1^{--}}$ Omn\`es function "
	with modernplotting.toolkit.PdfWriter(outFileName) as pdfOutput:
		plot = style.getPlot1D()
		plot.plot([v[0] for v in data],[v[1]for v in data], markersize = 0.)
		plot.plot([v[0] for v in data],[v[2]for v in data], markersize = 0.)
		plot.setXlim(0.,3.)
		plot.setYlim(-3.,7.)
		plot.setXlabel(r"$s\,[(\text{GeV}/c^2)^2]$")
		plot.setYlabel(r"$\Re/\Im(\Omega)$")
		plot.finishAndSaveAndClose(pdfOutput)


if __name__ == "__main__":
	main()
