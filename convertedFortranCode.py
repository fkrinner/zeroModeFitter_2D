#!/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/Python_ultra/Python-2.7.10/bin/python
# convertedFortranCode.py
# Created: 2018-07-18 14:30:20.563736
# Author: Fabian Krinner
import os, sys
from cmath import exp as cdexp
from math import pi
from math import log as dlog
mpi  = 0.13957018e0
mk   = 0.4957e0
mu   = 0.770e0
tCut = 4.e0

def dabs(x):
	return abs(x)

def strr(v):
	return str(v.real)

def stri(v):
	return str(v.imag)

def Ftilde(s,param):
 # # # # # THIS PARAMETER ORDERING IS FIXED NOW... NEVER EVER CHANGE # # # !!!
	M          = param[ 2]
	G          = param[ 3]
	Mp         = param[ 4]
	Gp         = param[ 5]
	M2p        = param[ 6]
	G2p        = param[ 7]
	alp        = param[ 8]
	phip       = param[ 9]
	al2p       = param[10]
	phi2p      = param[11]

	i = 1.j

	part1 = (M**2+s*(alp*cdexp(i*phip)+al2p*cdexp(i*phi2p)))/(M**2-s+kappa(s,M,G, True)*(A(s,mpi,mu)+0.5e0*A(s,mk,mu))-i*M*Gam(s,M,G, True))
	part2 = -s*alp*cdexp(i*phip)/(Mp**2-s+kappa(s,Mp,Gp)*A(s,mpi,mu)-i*Mp*Gam(s,Mp,Gp))
	part3 =  -s*al2p*cdexp(i*phi2p)/(M2p**2-s+kappa(s,M2p,G2p)*A(s,mpi,mu)-i*M2p*Gam(s,M2p,G2p))

	return part1 + part2 + part3

def Gam(s,mP,Gp, isRho = False):
	if isRho:
	      return Gp*s/mP**2*((sigma(s,mpi))**3+ 1.e0/2.e0*(sigma(s,mk))**3)/ ((sigma(mP**2,mpi))**3+1.e0/2.e0*(sigma(mP**2,mk))**3)
	else:
		return Gp*s/mP**2*(sigma(s,mpi))**3/(sigma(mP**2,mpi))**3

def kappa(s,Mp,Gp, isRho = False):
	if isRho:
		return Gp/Mp*s/(pi*(sigma(Mp**2,mpi))**3+ 1.e0/2.e0*(sigma(Mp**2,mk))**3)
	else:
		return Gp/Mp*s/(pi*sigma(Mp**2,mpi)**3)

def sigma(s,mP):
	if not s > 4.e0*mP**2:
		return 0.e0
	else:
		return (1.e0-4.e0*mP**2/s)**.5

def A(s,mP,mu):
	if not s > 4.e0*mP**2:
		return dlog(mP**2/mu**2)+8.e0*mP**2/s-5.e0/3.e0
	else:
		return dlog(mP**2/mu**2)+8.e0*mP**2/s-5.e0/3.e0+(sigma(s,mP))**3*dlog(dabs((sigma(s,mP)+1.e0)/(sigma(s,mP)-1.e0)))

def main():
	s = 2.3
	params = [0.36674e-01, 0.31230e-02, .83386, .19771, 1.4974, .78518, 1.6855, .80109, 0.17254,-0.97591, 0.23374, 2.1982]
	with open("ftilde.dat", 'r') as inFile, open("ftilde_copy.dat", 'w') as outFile:
		for line in inFile.readlines():
			s = float(line.split()[0])
			val = Ftilde(s, params)
			outFile.write(str(s) + ' ' + str(val.real) + ' ' + str(val.imag) + '\n')

if __name__ == "__main__":
	main()
