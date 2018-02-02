import matplotlib
matplotlib.use("Agg")
import modernplotting.mpplot
import modernplotting.toolkit
import modernplotting.specialPlots as mpsp

import LaTeX_strings
from utils import checkLaTeX

from makeResonanceParameterPlots import makeLineFit

def parseFile(inFileName, fakk = 1):
	data  = []
	nVals = None
	with open(inFileName, 'r') as inFile:
		for line in inFile.readlines():
			chunks = line.split()
			if not nVals:
				nVals = len(chunks)
			else:
				if not nVals == len(chunks):
					raise ValueError("Numbers of values do not match")
			if nVals == 0:
				raise ValueError("No values found")
			vals = [int(chunks[0]), float(chunks[1])]
			vals.append(float(chunks[2])*fakk)
			vals.append(float(chunks[3])*fakk)
			for i in range(4,len(chunks)):
				vals.append(float(chunks[i]))
			data.append(vals)
	return data

def makeAvgVals(vals):
	m  = 0.
	wm = 0.
	G  = 0.
	wG = 0.
	for v in vals:
		m  += v[2]/v[3]
		wm += 1./v[3]
		G  += v[4]/v[5]
		wG += 1./v[5]
	return m/wm, G/wG

def makeResonancePlots(massFileNameBase, outFileNameBase, totalFile = None, LaTeXname = "",  mMin =0.96, mMax = 2.5, fakk = 1000, mLimMin = .700, mLimMax = 1.30, g1LimMin = .003, g1LimMax = .15, g2LimMin = -.10, g2LimMax = .20, isobLaTeX = "", maxErrG = None, maxErrM = None):
	if not "<mode>" in outFileNameBase:
		raise ValueError("No '<mode>' location specifiec in outFileNameBase")
	style       = modernplotting.mpplot.PlotterStyle()
	style.titleRight = LaTeXname
	colors      = [modernplotting.colors.colorScheme.blue,modernplotting.colors.colorScheme.red,modernplotting.colors.colorScheme.green,modernplotting.colors.colorScheme.orange]
	if totalFile:
		totalVals   = parseFile(totalFile, fakk = fakk)
	massFileName = outFileNameBase.replace("<mode>","mass")
	print LaTeXname
	with  modernplotting.toolkit.PdfWriter(massFileName) as pdfOutput:
		plot = style.getPlot1D()
		plot.setXlabel(LaTeX_strings.m3Pi)
		plot.setYlim(mLimMin*fakk,mLimMax*fakk)
		plot.setXlim(mMin, mMax)
		plot.setYlabel(r"$m_{"+isobLaTeX+r"}\,[\text{MeV}/c^2]$")
		for i in range(4):
			if totalFile:
				plot.plot([mMin, mMax], [totalVals[i][2],totalVals[i][2]], color = colors[i], markersize = 0.)
			massFileName = massFileNameBase.replace("<tBin>", str(i))
			tBinVals     = parseFile(massFileName, fakk = fakk)
			x = [tBinVals[j][1] for j in range(len(tBinVals))]
			y = [tBinVals[j][2] for j in range(len(tBinVals))]
			e = [tBinVals[j][3] for j in range(len(tBinVals))]

			if maxErrM:
				nX = []
				nY = []
				nE = []
				for iE,err in enumerate(e):
					if err < maxErrM*fakk:
						nX.append(x[iE])
						nY.append(y[iE])
						nE.append(err)
				x = nX
				y = nY
				e = nE
			plot.axes.errorbar(x,y, yerr = e, color = colors[i], linestyle =  " ")
			plot.plot(x,y, color = colors[i], markeredgecolor = 'none', linestyle =  " ", markersize = 2.)

			fit,c2fit,c2global = makeLineFit(totalVals[i][2],y,e)
			if c2fit < 10.:
				print "tBin",i, "Fitted mass:",fit,"with a chi2/ndf of:",c2fit," Global mass:",totalVals[i][2],"with a chi2/ndf of:",c2global
		plot.finishAndSaveAndClose(pdfOutput)


	g1FileName = outFileNameBase.replace("<mode>","g1")
	with  modernplotting.toolkit.PdfWriter(g1FileName) as pdfOutput:
		plot = style.getPlot1D()
		plot.setXlabel(LaTeX_strings.m3Pi)
		plot.setYlim(g1LimMin, g1LimMax)
		plot.setYlabel(r"$g_\pi\,[(\text{GeV}/c^2)^2]$")
		plot.setXlim(mMin, mMax)
		for i in range(4):
			if totalFile:
				plot.plot([mMin, mMax], [totalVals[i][4],totalVals[i][4]], color = colors[i], markersize = 0.)
			massFileName = massFileNameBase.replace("<tBin>", str(i))
			tBinVals     = parseFile(massFileName, fakk = fakk)
			x = [tBinVals[j][1] for j in range(len(tBinVals))]
			y = [tBinVals[j][4] for j in range(len(tBinVals))]
			e = [tBinVals[j][5] for j in range(len(tBinVals))]
			if maxErrG:
				nX = []
				nY = []
				nE = []
				for iE,err in enumerate(e):
					if err < maxErrG:
						nX.append(x[iE])
						nY.append(y[iE])
						nE.append(err)
				x = nX
				y = nY
				e = nE
			plot.axes.errorbar(x,y, yerr = e, color = colors[i], linestyle =  " ")
			plot.plot(x,y, color = colors[i], markeredgecolor = 'none', linestyle =  " ", markersize = 2.)

			fit,c2fit,c2global = makeLineFit(totalVals[i][4],y,e)
			if c2fit < 10.:
				print "tBin",i, "Fitted g1:",fit,"with a chi2/ndf of:",c2fit," Global g1:",totalVals[i][4],"with a chi2/ndf of:",c2global
		plot.finishAndSaveAndClose(pdfOutput)

	g2FileName = outFileNameBase.replace("<mode>","g2")
	with  modernplotting.toolkit.PdfWriter(g2FileName) as pdfOutput:
		plot = style.getPlot1D()
		plot.setXlabel(LaTeX_strings.m3Pi)
		plot.setYlim(g2LimMin, g2LimMax)
		plot.setYlabel(r"$g_\text{K}\,[(\text{GeV}/c^2)^2]$")
		plot.setXlim(mMin, mMax)
		for i in range(4):
			if totalFile:
				plot.plot([mMin, mMax], [totalVals[i][6],totalVals[i][6]], color = colors[i], markersize = 0.)
			massFileName = massFileNameBase.replace("<tBin>", str(i))
			tBinVals     = parseFile(massFileName, fakk = fakk)
			x = [tBinVals[j][1] for j in range(len(tBinVals))]
			y = [tBinVals[j][6] for j in range(len(tBinVals))]
			e = [tBinVals[j][7] for j in range(len(tBinVals))]
			if maxErrG:
				nX = []
				nY = []
				nE = []
				for iE,err in enumerate(e):
					if err < maxErrG:
						nX.append(x[iE])
						nY.append(y[iE])
						nE.append(err)
				x = nX
				y = nY
				e = nE
			plot.axes.errorbar(x,y, yerr = e, color = colors[i], linestyle =  " ")
			plot.plot(x,y, color = colors[i], markeredgecolor = 'none', linestyle =  " ", markersize = 2.)

			fit,c2fit,c2global = makeLineFit(totalVals[i][6],y,e)
			if c2fit < 10.:
				print "tBin",i, "Fitted mass:",fit,"with a chi2/ndf of:",c2fit," Global mass:",totalVals[i][6],"with a chi2/ndf of:",c2global

		plot.finishAndSaveAndClose(pdfOutput)


	chi2FileName = outFileNameBase.replace("<mode>","chi2")
	with  modernplotting.toolkit.PdfWriter(chi2FileName) as pdfOutput:
		plot = style.getPlot1D()
		plot.setXlabel(LaTeX_strings.m3Pi)
		plot.setYlabel(r"$\chi^2/\text{ndf}$")
		plot.setYlim(0.,20.)
		plot.setXlim(mMin, mMax)
		for i in range(4):
			massFileName = massFileNameBase.replace("<tBin>", str(i))
			tBinVals     = parseFile(massFileName, fakk = fakk)
			x = [tBinVals[j][1] for j in range(len(tBinVals))]
			y = [tBinVals[j][8] for j in range(len(tBinVals))]
#			print "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
#			print x
#			print y
#			print "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
			plot.plot(x,y, color = colors[i], linestyle =  " ", markeredgecolor = 'none' )
		plot.finishAndSaveAndClose(pdfOutput)


def main():
	checkLaTeX()
	waves = [
	("./f0MassesAndWidthForScript/f0MassesAndWidths_0mp_<tBin>.dat","./massWidthPlots/0mp0p0pp_<mode>.pdf","./f0MassesAndWidthForScript/f0MassesAndWidths_0-+0+0++_global.dat",r"$0^{-+}0^+[\pi\pi]_{0^{++}}\pi$S"),
	("./f0MassesAndWidthForScript/f0MassesAndWidths_1pp_<tBin>.dat","./massWidthPlots/1pp0p0pp_<mode>.pdf","./f0MassesAndWidthForScript/f0MassesAndWidths_1++0+0++_global.dat",r"$1^{++}0^+[\pi\pi]_{0^{++}}\pi$P"),
	("./f0MassesAndWidthForScript/f0MassesAndWidths_2mp_<tBin>.dat","./massWidthPlots/2mp0p0pp_<mode>.pdf","./f0MassesAndWidthForScript/f0MassesAndWidths_2-+0+0++_global.dat",r"$2^{-+}0^+[\pi\pi]_{0^{++}}\pi$D")
	]
	for wave in waves:
		makeResonancePlots(*wave, isobLaTeX = r"\text{f}_0", mLimMin = 0.95, mLimMax = 1.00, maxErrG = 0.15, maxErrM = 0.15)
	
if __name__ == "__main__":
	main()
