#!/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/Python_ultra/Python-2.7.10/bin/python
# plotter.py
# Created: 2018-04-06 10:53:06.182634
# Author: Fabian Krinner
import os, sys
from physUtils import chewMandelstam
import ROOT

def iRho(s,m1,m2):
	sth = (m1+m2)**2
	return (sth/s - 1.+0.j)**.5

def main():
	histRe = ROOT.TH1D("hr","hr",1000, -1., 5.)
	histIm = ROOT.TH1D("hi","hi",1000, -1., 5.)
	for i in range(1000):
		s   = histRe.GetXaxis().GetBinCenter(i+1)
		val = chewMandelstam(s+0.00001j, .5 , .5) - chewMandelstam(0.00001j, .5 , .5)
		histRe.SetBinContent(i+1,val.real)
		histIm.SetBinContent(i+1,val.imag)
	histRe.Draw()
	histIm.Draw("SAME")
	raw_input()
	return


	m1 = 0.1
	m2 = 0.2

	c1 = ROOT.TCanvas()
	c2 = ROOT.TCanvas()
	c3 = ROOT.TCanvas()
	c4 = ROOT.TCanvas()

	nX = 1000
	nY = 1000

	CMhistRe = ROOT.TH2D("re_CM","re_CM",nX,-.2,2., nY, -1.,1.)
	CMhistIm = ROOT.TH2D("im_CM","im_CM",nX,-.2,2., nY, -1.,1.)

	PShistRe = ROOT.TH2D("re_PS","re_PS",nX,-.2,2., nY, -1.,1.)
	PShistIm = ROOT.TH2D("im_PS","im_PS",nX,-.2,2., nY, -1.,1.)

	for i in range(nX):
		x = CMhistRe.GetXaxis().GetBinCenter(i+1)
		for j in range(nY):
			y   = CMhistRe.GetYaxis().GetBinCenter(j+1)
			val = chewMandelstam(x+1.j*y, m1, m2)
			CMhistRe.SetBinContent(i+1,j+1,val.real)
			CMhistIm.SetBinContent(i+1,j+1,val.imag)
			val = iRho(x+1.j*y, m1, m2)
			PShistRe.SetBinContent(i+1,j+1,val.real)
			PShistIm.SetBinContent(i+1,j+1,val.imag)
	c1.cd()
	CMhistRe.Draw("colz")
	c2.cd()
	CMhistIm.Draw("colz")
	c3.cd()
	PShistRe.Draw("colz")
	c4.cd()
	PShistIm.Draw("colz")
	raw_input()


if __name__ == "__main__":
	main()
