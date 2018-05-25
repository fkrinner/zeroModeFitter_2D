# makeChi2DiffPLots.py
import os, sys
import ROOT

def getData(inFileName = "chi2comparisons"):
	data = {}
	with open(inFileName) as inFile:
		for line in inFile.readlines():
			chunks = line.split()
			t = int(chunks[0])
			m = int(chunks[1])
			full = float(chunks[3])
			reim = float(chunks[4])
			none = float(chunks[5])
			data[(t,m)] = (full,reim,none)
	return data

def makeHists(data, nT = 4, nM = 50):
	histReim = ROOT.TH2D("reim","reim", nM, 0., 1., nT, 0., 1.)
	histNone = ROOT.TH2D("none","none", nM, 0., 1., nT, 0., 1.)
	for t in range(nT):
		for m in range(nM):
			histReim.SetBinContent(m+1, t+1, data[(t,m)][1]/data[(t,m)][0])
			histNone.SetBinContent(m+1, t+1, data[(t,m)][2]/data[(t,m)][0])
	return histReim, histNone

def main():
	ROOT.gStyle.SetOptStat(0)
	data = getData()
	h1, h2 = makeHists(data)
	h1.Draw("colz")
	raw_input()
	h2.Draw("colz")
	raw_input()

if __name__ == "__main__":
	main()
