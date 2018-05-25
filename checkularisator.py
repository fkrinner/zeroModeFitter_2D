# checkularisator.py
# Created: 2018-03-05 14:54:54.279890
# Author: Fabian Krinner
import os, sys
import ROOT
from rootfabi import root_open

def main():
	histNameBase = "theoTotal_<mode>_1-+1+[pi,pi]1--PiP_0"
	inFileName   = "1mp1p_totals_t0.root"
	with root_open(inFileName, "READ") as inFile:
		intens = inFile.Get(histNameBase.replace("<mode>", "intensity"))
		real   = inFile.Get(histNameBase.replace("<mode>", "real"))
		imag   = inFile.Get(histNameBase.replace("<mode>", "imag"))
		corr   = inFile.Get(histNameBase.replace("<mode>", "real_imag_correlation"))
		for i in range(intens.GetNbinsX()):
			I    = intens.GetBinContent(i+1)
			J    = real.GetBinContent(i+1)**2 + imag.GetBinContent(i+1)**2
			coma = [[real.GetBinError(i+1)**2, corr.GetBinContent(i+1)],[corr.GetBinContent(i+1), imag.GetBinError(i+1)**2]]
			jac  = [2*real.GetBinContent(i+1),2*imag.GetBinContent(i+1)]
			err  = 0.
			for n in range(2):
				for m in range(2):
					err += jac[m]*jac[n]*coma[m][n]
			print err**.5, intens.GetBinError(i+1)
			print I,J,"|||||||||||||||||||||||||}}{||"

if __name__ == "__main__":
	main()
