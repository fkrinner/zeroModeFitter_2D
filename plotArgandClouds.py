# plotArgandClouds.py
import matplotlib
matplotlib.use("Agg")
import modernplotting.root
import modernplotting.mpplot
import modernplotting.colors
import modernplotting.toolkit
import modernplotting.specialPlots

import os, sys

def loadData(dataFileName):
	data = []
	with open(dataFileName, 'r') as inFile:
		for line in inFile.readlines():
			dp = [float(v) for v in line.split()]
			if len(dp) > 1:
				data.append(dp)
	return data

def getPointArgand(dataFileName):
	data = loadData(dataFileName)
	re   = []
	im   = []
	for dp in data:
		re.append(dp[3])
		im.append(dp[5])
	return re, im

def getEllipseData(dataFileName, comaFileName):
	re,im = getPointArgand(dataFileName)
	comas = []
	with open(comaFileName, 'r') as inFile:
		for i,line in enumerate(inFile.readlines()):
			if i >= len(re):
				break
			vals = [float(v) for v in line.split()]
			coma = [[vals[0], vals[1]],[vals[1],vals[2]]]
			comas.append(coma)
	return re, im, comas

def getAvg(res, ims):
	nBin = len(res[0])
	nPts = len(res)
	avgsRe = [0.]*nBin
	avgsIm = [0.]*nBin
	for p in range(nPts):
		for b in range(nBin):
			avgsRe[b] += res[p][b]
			avgsIm[b] += ims[p][b]
	for b in range(nBin):
		avgsRe[b] /= nPts
		avgsIm[b] /= nPts
	comas = [[[0.,0.],[0.,0.]] for i in range(nBin)]
	for p in range(nPts):
		for b in range(nBin):
			delRe = res[p][b] - avgsRe[b]
			delIm = ims[p][b] - avgsIm[b]
			comas[b][0][0] += delRe**2
			comas[b][1][0] += delRe*delIm
			comas[b][0][1] += delRe*delIm
			comas[b][1][1] += delIm**2
	for b in range(nBin):
		comas[b][0][0] /= nPts-1
		comas[b][0][1] /= nPts-1
		comas[b][1][0] /= nPts-1
		comas[b][1][1] /= nPts-1
#		comas[b][0][0] **= .5
#		comas[b][0][1] **= .5
#		comas[b][1][0] **= .5
#		comas[b][1][1] **= .5
	return avgsRe, avgsIm, comas

lineArgs = {'marker': None, 'linestyle':'dashed', 'markersize':0, 'linewidth':.4, 'zorder' :0, 'color':'.5'}

def makePointPlot(outFileName, inFileNameBase, comaFileName = None):
	style        = modernplotting.mpplot.PlotterStyle()
	plot         = style.getPlot1D()
	colorizedBin = 10
	colorRe      = []
	colorIm      = []
	res = []
	ims = []
	with modernplotting.toolkit.PdfWriter(outFileName) as pdfOutput:
		nCl   = 1
		while True:
			inFileName = inFileNameBase.replace("<n>", str(nCl))
			if not os.path.isfile(inFileName):
				break
			re, im = getPointArgand(inFileName)
			res.append(re)
			ims.append(im)
			if colorizedBin:
				colorRe.append(re[colorizedBin])
				colorIm.append(im[colorizedBin])
			if nCl == 0:
				pass
#				plot.plot(re,im,linewidth=0.5,zorder=3,color='orange',markersize=2.,linestyle='solid',markeredgecolor='none')
			else:
				pass
#				plot.plot(re,im,linewidth=0.05,zorder=1,color='black',markersize=1.,linestyle='solid',markeredgecolor='none')
			nCl += 1
		if colorizedBin:
			plot.plot(colorRe,colorIm,linewidth=0.,zorder=2,color='red',markersize=2.,linestyle='solid',markeredgecolor='none')
		reAvg, imAvg, comasAvg = getAvg(res, ims)
		style.errorEllipsesEdgeColor = "red"
		style.errorEllipsesFaceColor = "red"
		modernplotting.specialPlots.plotErrorEllipses(plot,reAvg,imAvg,comasAvg,markerDefinitions={'linestyle':'solid','linewidth':.5,'zorder':0,'color':"red",'markersize':2.})
		style.errorEllipsesEdgeColor = "orange"
		style.errorEllipsesFaceColor = "orange"
		if comaFileName:
			re,im,comas=getEllipseData(inFileNameBase.replace("<n>",str(0)),comaFileName)
			modernplotting.specialPlots.plotErrorEllipses(plot,re,im,comas,markerDefinitions={'linestyle':'solid','linewidth':.5,'zorder':0,'color':"orange",'markersize':2.})
		pdfOutput.savefigAndClose()

def main():
	for mode in ['C','D']:
		outFileName    = "./argandDist_"+mode+".pdf"
		inFileNameBase = "./samplingArgands/argand_"+mode+"_<n>.dat"
		comaFileName   = "./correlations_"+mode+".dat"
		makePointPlot(outFileName, inFileNameBase, comaFileName = comaFileName)

if __name__ == "__main__":
	main()
