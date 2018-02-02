# a3fit.py
import os, sys
from numpy import array
import numpy.linalg as la
from cmath import phase

def parseFile(fileName): 
	with open(fileName, 'r') as inFile:
		first = True
		hessCommand = None
		for line in inFile.readlines():
			if first:
				first = False
				cmd = "ampl = "+line
				exec(cmd)
			else:
				if not hessCommand:
					hessCommand = "hess = "
				hessCommand += line
		exec(hessCommand)
	if not len(ampl[0]) == len(hess[0]):
		raise IOError("Dimension between ampl and hess")
	return ampl[0],hess[0]

def main():
	inFileName  = "./3pp_rho_cpls_0.dat"
	ampls, hess = parseFile(inFileName)
	hessInv     = [la.inv(h) for h in hess]
	nn          = len(ampls)
	with open('argand.deleteMe', 'w') as outFile:
		for i in range(len(ampls)):
			outFile.write(str(ampls[i][0]) + ' ' + str((2*hessInv[i][0,0])**.5) + ' ' + str(ampls[i][1]) + ' ' + str((2*hessInv[i][1,1])**.5) + '\n')
	with open('intens.deleteMe', 'w') as outFile:
		for i in range(len(ampls)):
			bin = 50 - nn + i
			m = 0.52 + 0.04 * bin
			outFile.write(str(bin) + ' ' + str(m) + ' ' + str(ampls[i][0]**2+ampls[i][1]**2) + ' ' + str(phase(ampls[i][0]+ 1.j*ampls[i][1])) + '\n')

if __name__ == "__main__":
	main()
