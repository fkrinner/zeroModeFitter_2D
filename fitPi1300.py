# fitPi1300.py
import os, sys
from cmath import phase
def readFile(inFileName):
	dataPoints = []
	with open(inFileName, 'r') as inFile:
		for line in inFile.readlines():
			vals = [float(v) for v in line.replace('[',' ').replace(']', ' ').replace(',',' ').split()]
			dataPoints.append(vals)
	return dataPoints

def main():
	inFileName  = "./0mp_rho_cpls_3.dat"
	outFileName = "pi1300.dat"
	vals        =  readFile(inFileName)
	with open(outFileName, 'w') as outFile:
		for i,v in enumerate(vals):
			outFile.write(str(i)  + ' ' + str(v[0]**2 + v[1]**2)+ ' ' + str(phase(v[0] + 1.j*v[1])) + ' ' +str(v[0]) + ' ' + str(v[1]) + '\n')


if __name__ == "__main__":
	main()
