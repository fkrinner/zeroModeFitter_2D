#!/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/Python_ultra/Python-2.7.10/bin/python
# parsePhases.py
# Created: 2018-06-05 09:32:19.654953
# Author: Fabian Krinner
import os, sys
from cmath import phase

def main():
	inFileName = "./couplingsFromTheGood"
	phases     = []
	t = None
	m = None
	with open(inFileName, 'r') as inFile:
		for line in inFile.readlines():
			if ">>t<<" in line:
				t = int(line.split()[1])
				continue
			if ">>m<<" in line:
				m = int(line.split()[1])
				continue
			val = float(line.split()[0]) + 1.j * float(line.split()[1])
			phases.append(phase(val))

	print phases

if __name__ == "__main__":
	main()
