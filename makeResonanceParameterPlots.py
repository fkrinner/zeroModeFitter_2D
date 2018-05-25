import matplotlib
matplotlib.use("Agg")
import modernplotting.mpplot
import modernplotting.toolkit
import modernplotting.specialPlots as mpsp

import LaTeX_strings
from utils import checkLaTeX

from math import isnan

def parseFile(inFileName, fakk = 1):
	data  = []
	nVals = None
	with open(inFileName, 'r') as inFile:
		for line in inFile.readlines():
			chunks = line.split()
			if nVals in None:
				nVals = len(chunks)
			else:
				if not nVals == len(chunks):
					raise ValueError("Numbers of values do not match")
			if nVals == 0:
				raise ValueError("No values found")
			vals = [int(chunks[0]), float(chunks[1])]
			for i in range(2,6):
				vals.append(float(chunks[i])*fakk)
			if len(chunks) > 6:
				vals.append(float(chunks[6]))
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

def makeLineFit(glob,y,e):
	c2global = 0.
	ndf      = 0
	weight   = 0.
	avg      = 0.
	for i in range(len(y)):
#		if ((glob-y[i])/e[i])**2 > 1000.:
#			print i,((glob-y[i])/e[i])**2, glob, y[i],e[i]
		if e[i] == 0. or isnan(e[i]) or isnan(y[i]):
			continue
		ndf      += 1
		c2global += ((glob-y[i])/e[i])**2
		avg      += y[i]/e[i]
		weight   += 1./e[i]
	avg/=weight
	c2fit = 0.
	for i in range(len(y)):
		if e[i] == 0. or isnan(e[i]) or isnan(y[i]):
			continue
		c2fit += ((y[i] - avg)/e[i])**2
#	print ndf
	return avg, c2fit/(ndf-1), c2global/ndf


def makeResonancePlots(massFileNameBase, outFileNameBase, totalFile = None, LaTeXwaveName = "", mMin = 1.0, mMax = 2.5, fakk = 1000, mLimMin = .700, mLimMax = .850, GLimMin = .050, GLimMax = .250, isobLaTeX = "", maxErrVal = None, PDG = []):
	print outFileNameBase,"<<<<<<<<<<<<<<<<<<<<<<<<<<><<<<<<<<<<<<<<<<<<<<<<<<<<<><<<<<<<<<<<<<<<<<<<<<<>"
	if not "<mode>" in outFileNameBase:
		raise ValueError("No '<mode>' location specifiec in outFileNameBase")
	style       = modernplotting.mpplot.PlotterStyle()
	style.titleRight = LaTeXwaveName
	color1 = modernplotting.colors.colorScheme.blue
	color2 = modernplotting.colors.colorScheme.red
	color3 = modernplotting.colors.colorScheme.green
	color4 = modernplotting.colors.colorScheme.orange
#	color4 = modernplotting.colors.makeColorDarker(modernplotting.colors.colorScheme.orange,.2)

	colors      = [color1,color2,color3,color4]
	if totalFile is not None and not totalFile == "":
		totalVals   = parseFile(totalFile, fakk = fakk)
	massFileName = outFileNameBase.replace("<mode>","mass")
	with  modernplotting.toolkit.PdfWriter(massFileName) as pdfOutput:
		plot = style.getPlot1D()
		plot.setXlabel(LaTeX_strings.m3Pi)
		plot.setYlim(mLimMin*fakk,mLimMax*fakk)
		plot.setXlim(mMin, mMax)
		plot.setYlabel(r"$m_{"+isobLaTeX+r"}\,[\text{MeV}/c^2]$")
		for i in range(4):
			if totalFile is not None and not totalFile == "":
				plot.plot([mMin, mMax], [totalVals[i][2],totalVals[i][2]], color = colors[i], markersize = 0.)
			massFileName = massFileNameBase.replace("<tBin>", str(i))
			tBinVals     = parseFile(massFileName, fakk = fakk)
			x = [tBinVals[j][1] for j in range(len(tBinVals))]
			y = [tBinVals[j][2] for j in range(len(tBinVals))]
			e = [tBinVals[j][3] for j in range(len(tBinVals))]
			if maxErrVal:
				nX = []
				nY = []
				nE = []
				for iE,err in enumerate(e):
					if err < maxErrVal*fakk:
						nX.append(x[iE])
						nY.append(y[iE])
						nE.append(err)
				x = nX
				y = nY
				e = nE
			plot.axes.errorbar(x,y, yerr = e, color = colors[i], linestyle =  " ")
			plot.plot(x,y, color = colors[i], markeredgecolor = 'none', linestyle =  " ", markersize = 2. )
#			if totalFile and not totalFile == "":
#				fit,c2fit,c2global = makeLineFit(totalVals[i][2],y,e)
#				if c2fit < 10.:
#					print "tBin",i, "Fitted mass:",fit,"with a chi2/ndf of:",c2fit," Global mass:",totalVals[i][2],"with a chi2/ndf of:",c2global
		if len(PDG) == 4:
			xLim = plot.axes.get_xlim()
			band = matplotlib.patches.Rectangle((xLim[0], (PDG[0]-PDG[1])*fakk), xLim[1]-xLim[0], 2*PDG[1]*fakk, color = "gray")
			plot.axes.add_patch(band)
		plot.finishAndSaveAndClose(pdfOutput)
	widthFileName = outFileNameBase.replace("<mode>","width")
	print ' - - - - - - - - - - - -'
	with  modernplotting.toolkit.PdfWriter(widthFileName) as pdfOutput:
		plot = style.getPlot1D()
		plot.setXlabel(LaTeX_strings.m3Pi)
		plot.setYlim(GLimMin*fakk,GLimMax*fakk)
		plot.setYlabel(r"$\Gamma_{"+isobLaTeX+r"}\,[\text{MeV}/c^2]$")
		plot.setXlim(mMin, mMax)
		for i in range(4):
			if totalFile is not None and not totalFile == "":
				plot.plot([mMin, mMax], [totalVals[i][4],totalVals[i][4]], color = colors[i], markersize = 0.)
			massFileName = massFileNameBase.replace("<tBin>", str(i))
			tBinVals     = parseFile(massFileName, fakk = fakk)
			x = [tBinVals[j][1] for j in range(len(tBinVals))]
			y = [tBinVals[j][4] for j in range(len(tBinVals))]
			e = [tBinVals[j][5] for j in range(len(tBinVals))]
			if maxErrVal:
				nX = []
				nY = []
				nE = []
				for iE,err in enumerate(e):
					if err < maxErrVal*fakk:
						nX.append(x[iE])
						nY.append(y[iE])
						nE.append(err)
				x = nX
				y = nY
				e = nE

			plot.axes.errorbar(x,y, yerr = e, color = colors[i], linestyle =  " ")
			plot.plot(x,y, color = colors[i], markeredgecolor = 'none', linestyle =  " ", markersize = 2.)
#			if totalFile and not totalFile == "":
#				fit,c2fit,c2global = makeLineFit(totalVals[i][4],y,e)
#				if c2fit < 10.:
#					print "tBin",i,"Fitted width:",fit,"with a chi2/ndf of:",c2fit," Global width:",totalVals[i][4],"with a chi2/ndf of:",c2global
		if len(PDG) == 4:
			xLim = plot.axes.get_xlim()
			band = matplotlib.patches.Rectangle((xLim[0], (PDG[2]-PDG[3])*fakk), xLim[1]-xLim[0], 2*PDG[3]*fakk, color = "gray")
			plot.axes.add_patch(band)

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
			y = [tBinVals[j][6] for j in range(len(tBinVals))]
			plot.plot(x,y, color = colors[i], linestyle =  " ", markeredgecolor = 'none')
		plot.finishAndSaveAndClose(pdfOutput)


def main():
	checkLaTeX()
	waves = [
	("./rhoMassesAndWidthForScript/rhoMassesAndWidths_0mp_<tBin>.dat","./massWidthPlots/0mp0pRho_<mode>.pdf","./rhoMassesAndWidthForScript/rhoMassesAndWidths_0-+0+1--_global.dat",r"$0^{-+}0^+[\pi\pi]_{1^{--}}\pi$P"),
	("./rhoMassesAndWidthForScript/rhoMassesAndWidths_0mp_<tBin>_range1.2.dat","./massWidthPlots/0mp0pRho_<mode>_range1.2.pdf","./rhoMassesAndWidthForScript/rhoMassesAndWidths_0-+0+1--_global_range1.2.dat",r"$0^{-+}0^+[\pi\pi]_{1^{--}}\pi$P"),
	("./rhoMassesAndWidthForScript/rhoMassesAndWidths_1pp_<tBin>.dat","./massWidthPlots/1pp0pRho_<mode>.pdf","./rhoMassesAndWidthForScript/rhoMassesAndWidths_1++0+1--_global.dat",r"$1^{++}0^+[\pi\pi]_{1^{--}}\pi$S"),
	("./rhoMassesAndWidthForScript/rhoMassesAndWidths_1pp_<tBin>_range1.2.dat","./massWidthPlots/1pp0pRho_<mode>_range1.2.pdf","./rhoMassesAndWidthForScript/rhoMassesAndWidths_1++0+1--_global_range1.2.dat",r"$1^{++}0^+[\pi\pi]_{1^{--}}\pi$S"),
#	("./rhoMassesAndWidthForScript/rhoMassesAndWidths_2mp_<tBin>.dat","./massWidthPlots/2mp0pRho_<mode>.pdf","./rhoMassesAndWidthForScript/rhoMassesAndWidths_2-+0+1--_global.dat",r"$2^{-+}0^+[\pi\pi]_{1^{--}}\pi$P/F"),
#	("./rhoMassesAndWidthForScript/rhoMassesAndWidths_2mpP_<tBin>.dat","./massWidthPlots/2mp0pRhoP_<mode>.pdf","./rhoMassesAndWidthForScript/rhoMassesAndWidths_2-+0+1--P_global.dat",r"$2^{-+}0^+[\pi\pi]_{1^{--}}\pi$P"),
#	("./rhoMassesAndWidthForScript/rhoMassesAndWidths_2mpP_<tBin>_range1.2.dat","./massWidthPlots/2mp0pRhoP_<mode>_range1.2.pdf","./rhoMassesAndWidthForScript/rhoMassesAndWidths_2-+0+1--P_global_range1.2.dat",r"$2^{-+}0^+[\pi\pi]_{1^{--}}\pi$P"),
#	("./rhoMassesAndWidthForScript/rhoMassesAndWidths_2mpF_<tBin>.dat","./massWidthPlots/2mp0pRhoF_<mode>.pdf","./rhoMassesAndWidthForScript/rhoMassesAndWidths_2-+0+1--F_global.dat",r"$2^{-+}0^+[\pi\pi]_{1^{--}}\pi$F"),
#	("./rhoMassesAndWidthForScript/rhoMassesAndWidths_2mpF_<tBin>_range1.2.dat","./massWidthPlots/2mp0pRhoF_<mode>_range1.2.pdf","./rhoMassesAndWidthForScript/rhoMassesAndWidths_2-+0+1--F_global_range1.2.dat",r"$2^{-+}0^+[\pi\pi]_{1^{--}}\pi$F"),
#	("./rhoMassesAndWidthForScript/rhoMassesAndWidths_1++1+_<tBin>.dat","./massWidthPlots/1pp1pRho_<mode>.pdf","./rhoMassesAndWidthForScript/rhoMassesAndWidths_1++1+_global.dat",r"$1^{++}1^+[\pi\pi]_{1^{--}}\pi$S"),
#	("./rhoMassesAndWidthForScript/rhoMassesAndWidths_2-+1+_<tBin>.dat","./massWidthPlots/2mp1pRho_<mode>.pdf","./rhoMassesAndWidthForScript/rhoMassesAndWidths_2-+1+_global.dat",r"$2^{-+}1^+[\pi\pi]_{1^{--}}\pi$P"),
#	("./rhoMassesAndWidthForScript/rhoMassesAndWidths_2++1+_<tBin>.dat","./massWidthPlots/2pp1pRho_<mode>.pdf","./rhoMassesAndWidthForScript/rhoMassesAndWidths_2++1+_global.dat",r"$2^{++}1^+[\pi\pi]_{1^{--}}\pi$D"),
#	("./rhoMassesAndWidthForScript/rhoMassesAndWidths_3pp_<tBin>.dat","./massWidthPlots/3pp0pRho_<mode>.pdf","./rhoMassesAndWidthForScript/rhoMassesAndWidths_3pp_global.dat",r"$3^{++}0^+[\pi\pi]_{1^{--}}\pi$D"),
#	("./rhoMassesAndWidthForScript/rhoMassesAndWidths_bigger1ppD_<tBin>.dat","./massWidthPlots/1pp0pRhoD_<mode>.pdf","./rhoMassesAndWidthForScript/rhoMassesAndWidths_bigger1++0+1--D_global.dat",r"$1^{++}0^+[\pi\pi]_{1^{--}}\pi$D"),
	("./rhoMassesAndWidthForScript/rhoMassesAndWidths_exotic_<tBin>.dat","./massWidthPlots/1mp1pRho_<mode>.pdf","./rhoMassesAndWidthForScript/rhoMassesAndWidths_1-+1+1--_global.dat",r"$1^{-+}1^+[\pi\pi]_{1^{--}}\pi$P"),
	("./rhoMassesAndWidthForScript/rhoMassesAndWidths_exotic_<tBin>_range1.2.dat","./massWidthPlots/1mp1pRho_<mode>_range1.2.pdf","./rhoMassesAndWidthForScript/rhoMassesAndWidths_1-+1+1--_global_range1.2.dat",r"$1^{-+}1^+[\pi\pi]_{1^{--}}\pi$P"),
#	("./rhoMassesAndWidthForScript/rhoMassesAndWidths_4++1+1--_<tBin>.dat","./massWidthPlots/4pp1pRho_<mode>.pdf","./rhoMassesAndWidthForScript/rhoMassesAndWidths_4++1+1--_global.dat",r"$4^{++}1^+[\pi\pi]_{1^{--}}\pi$G"),
#	("./rhoMassesAndWidthForScript/rhoMassesAndWidths_4-+0+_<tBin>.dat","./massWidthPlots/4mp0pRho_<mode>.pdf","./rhoMassesAndWidthForScript/rhoMassesAndWidths_4-+0+_global.dat",r"$4^{-+}0^+[\pi\pi]_{1^{--}}\pi$F"),
#	("./rhoMassesAndWidthForScript/rhoMassesAndWidths_6-+0+_<tBin>.dat","./massWidthPlots/6mp0pRho_<mode>.pdf","./rhoMassesAndWidthForScript/rhoMassesAndWidths_6-+0+_global.dat",r"$6^{-+}0^+[\pi\pi]_{1^{--}}\pi$H"),
	("./f2MassesAndWidthForScript/f2MassesAndWidths_2mp_<tBin>.dat","./massWidthPlots/2mpF2_<mode>.pdf","./f2MassesAndWidthForScript/f2MassesAndWidths_2-+0+2++_global.dat",r"$2^{-+}0^+[\pi\pi]_{2^{++}}\pi$S"),
#	("./f2MassesAndWidthForScript/f2MassesAndWidths_2mp_<tBin>_range1.8.dat","./massWidthPlots/2mpF2_<mode>_range1.8.pdf","./f2MassesAndWidthForScript/f2MassesAndWidths_2-+0+2++_global_range1.8.dat",r"$2^{-+}0^+[\pi\pi]_{2^{++}}\pi$S"),
#	("./f2MassesAndWidthForScript/f2MassesAndWidths_bigger1pp_<tBin>.dat","./massWidthPlots/bigger1ppF2_<mode>.pdf","./f2MassesAndWidthForScript/f2MassesAndWidths_bigger1++0+2++_global.dat",r"$1^{++}0^+[\pi\pi]_{2^{++}}\pi$P"),
	("./f2MassesAndWidthForScript/f2MassesAndWidths_bigger2mpD_<tBin>.dat","./massWidthPlots/bigger2mpF2D_<mode>.pdf","./f2MassesAndWidthForScript/f2MassesAndWidths_bigger2-+0+2++D_global.dat",r"$2^{-+}0^+[\pi\pi]_{2^{++}}\pi$D"),
#	("./f2MassesAndWidthForScript/f2MassesAndWidths_bigger2pp_<tBin>.dat","./massWidthPlots/bigger2ppF2_<mode>.pdf","./f2MassesAndWidthForScript/f2MassesAndWidths_bigger2++1+1++_global.dat",r"$2^{++}1^+[\pi\pi]_{2^{++}}\pi$P"),
#	("./f2MassesAndWidthForScript/f2MassesAndWidths_3++0+_<tBin>.dat","./massWidthPlots/bigger3ppF2_<mode>.pdf","./f2MassesAndWidthForScript/f2MassesAndWidths_3++0+_global.dat",r"$3^{++}0^+[\pi\pi]_{2^{++}}\pi$P"),
#	("./f2MassesAndWidthForScript/f2MassesAndWidths_4++1+2++_<tBin>.dat","./massWidthPlots/bigger4ppF2_<mode>.pdf","./f2MassesAndWidthForScript/f2MassesAndWidths_4++1+2++_global.dat",r"$4^{++}1^+[\pi\pi]_{2^{++}}\pi$P"),
#	("./f2MassesAndWidthForScript/f2MassesAndWidths_3++0+_<tBin>.dat","./massWidthPlots/bigger3ppF2_<mode>.pdf","./f2MassesAndWidthForScript/f2MassesAndWidths_2-+1+2++_global.dat",r"$2^{-+}1^+[\pi\pi]_{2^{++}}\pi$P"),
#	("./f2MassesAndWidthForScript/f2MassesAndWidths_2-+1+2++_<tBin>.dat","./massWidthPlots/bigger2mp1pSF2_<mode>.pdf","./f2MassesAndWidthForScript/f2MassesAndWidths_2-+1+2++_global.dat",r"$2^{-+}1^+[\pi\pi]_{2^{++}}\pi$P"),
	("./rho3MassesAndWidthForScript/rho3MassesAndWidths_3pp_<tBin>.dat","./massWidthPlots/3pp0pRho3_<mode>.pdf","./rho3MassesAndWidthForScript/rho3MassesAndWidths_3pp_global.dat",r"$3^{++}0^+[\pi\pi]_{3^{--}}\pi$S"),
	("./rhoPrimeMassesAndWidths/1pp_rhoPrimeMassesAndWidths_<tBin>.dat","./massWidthPlots/1pp0pRhoPrime_<mode>.pdf","",r"$1^{++}0^+[\pi\pi]_{1^{--}}\pi$S"),
	("./f0_1500massesAndWidths/0mp_f0_1500massesAndWidths_<tBin>.dat","./massWidthPlots/0mp0p_f0_1500_<mode>.pdf","./f0_1500massesAndWidths/0mp_f0_1500_massesAndWidths_global.dat",r"$0^{-+}0^+[\pi\pi]_{0^{++}}\pi$S"),
	("./f0_1500massesAndWidths/2mp_f0_1500massesAndWidths_<tBin>.dat","./massWidthPlots/2mp0p_f0_1500_<mode>.pdf","./f0_1500massesAndWidths/2mp_f0_1500_massesAndWidths_global.dat",r"$2^{-+}0^+[\pi\pi]_{0^{++}}\pi$D")
#	("./f0_1500massesAndWidths/0mp_f0_1500massesAndWidths_<tBin>_withConst.dat","./massWidthPlots/0mp0p_f0_1500_<mode>_withConst.pdf","./f0_1500massesAndWidths/0mp_f0_1500_massesAndWidths_global_withConst.dat",r"$0^{-+}0^+[\pi\pi]_{0^{++}}\pi$S"),
#	("./f0_1500massesAndWidths/2mp_f0_1500massesAndWidths_<tBin>_withConst.dat","./massWidthPlots/2mp0p_f0_1500_<mode>_withConst.pdf","./f0_1500massesAndWidths/2mp_f0_1500_massesAndWidths_global_withConst.dat",r"$2^{-+}0^+[\pi\pi]_{0^{++}}\pi$D")
	]


	PDG1500  = [1.504  , 0.006  , 0.109 , 0.007 ]
	PDGprime = [1.465  , 0.025  , 0.4   , 0.06  ]
	PDGrho3  = [1.6888 , 0.0021 , 0.161 , 0.01  ]
	PDGrho   = [0.77525, 0.00026, 0.1491, 0.0008]
	PDGf2    = [1.2757 , 0.0008 , 0.1867, 0.0025]

	for wave in waves:
		mMin1500   = 1.3
		GMin1500   = 0.
		range1500  =  .5

		mMinPrime  = 1.
		GMinPrime  = 0.
		rangePrime = 1.

		mMinRho3   = 1.4
		GMinRho3   = 0.
		rangeRho3  =  .6

		mMinRho    =  .730
		GMinRho    =  .1
		rangeRho   =  .100

		mMinF2     = 1.1
		GMinF2     =  .05
		rangeF2    =  .4

		maxErrVal  = .15

		mMax1500  = mMin1500  + range1500
		mMaxPrime = mMinPrime + rangePrime
		mMaxRho3  = mMinRho3  + rangeRho3
		mMaxRho   = mMinRho   + rangeRho
		mMaxF2    = mMinF2    + rangeF2

		GMax1500  = GMin1500  + range1500
		GMaxPrime = GMinPrime + rangePrime
		GMaxRho3  = GMinRho3  + rangeRho3
		GMaxRho   = GMinRho   + rangeRho
		GMaxF2    = GMinF2    + rangeF2

		if "1500"       in wave[0]:
			mMin1500 = 1.4
			mMax1500 = 1.7
			GMin1500 = 0.
			GMax1500 = 0.3

			makeResonancePlots(*wave, isobLaTeX = r"\mathrm{f}_0(1500)", mLimMin = mMin1500 , mLimMax = mMax1500 , GLimMin = GMin1500 , GLimMax = GMax1500 , maxErrVal = maxErrVal, PDG = PDG1500 )
		elif "rhoPrime" in wave[0]:
			mMinPrime = .5
			mMaxPrime = 2.5
			makeResonancePlots(*wave, isobLaTeX = r"\rho^\prime"       , mLimMin = mMinPrime, mLimMax = mMaxPrime, GLimMin = GMinPrime, GLimMax = GMaxPrime, maxErrVal = maxErrVal, PDG = PDGprime)
		elif "rho3"     in wave[0]:
			mMinRho3 = 1.5
			mMaxRho3 = 1.75
			GMinRho3 = 0.
			GMaxRho3 = 0.35
			makeResonancePlots(*wave, isobLaTeX = r"\rho_3"            , mLimMin = mMinRho3 , mLimMax = mMaxRho3 , GLimMin = GMinRho3 , GLimMax = GMaxRho3 , maxErrVal = maxErrVal, PDG = PDGrho3 )
		elif "rho"      in wave[0]:
			if wave[3] == r"$1^{++}0^+[\pi\pi]_{1^{--}}\pi$S":
				mMinRho = 0.76
				mMaxRho = 0.82
				GMinRho = 0.13
				GMaxRho = 0.25

			makeResonancePlots(*wave, isobLaTeX = r"\rho"              , mLimMin = mMinRho  , mLimMax = mMaxRho  , GLimMin = GMinRho  , GLimMax = GMaxRho  , maxErrVal = maxErrVal, PDG = PDGrho  )
		elif "f2"       in wave[0]:
			if wave[3].startswith(r"$2^{-+}0^+[\pi\pi]_{2^{++}}\pi$"):
				mMinF2 = 1.2
				mMaxF2 = 1.35
				GMinF2 = 0.1
				GMaxF2 = 0.3
			makeResonancePlots(*wave, isobLaTeX = r"\text{f}_2"        , mLimMin = mMinF2   , mLimMax = mMaxF2   , GLimMin = GMinF2   , GLimMax = GMaxF2   , maxErrVal = maxErrVal, PDG = PDGf2   )
		else:
			raise ValueError("Error in wave '" + str(wave) + "'")
	
if __name__ == "__main__":
	main()
