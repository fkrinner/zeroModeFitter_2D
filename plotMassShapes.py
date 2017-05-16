import parameterizationClasses as pc

def loadAmplFile(fn):
	amps = {}
	with open(fn, 'r') as inin:
		for line in inin.readlines():
			chunks = line.split()
			amps[float(chunks[0])] = float(chunks[1]) + 1.j * float(chunks[2])
	return amps

def compareFiles(fn1, fn2):
	amps1 = loadAmplFile(fn1)
	amps2 = loadAmplFile(fn2)
	with open('compare', 'w') as out:
		for key in amps1:
			if key in amps2:
				rat = amps1[key]/amps2[key]
				out.write(str(key) + ' ' + str(rat.real) + ' ' + str(rat.imag) + '\n')




def main():

	mPi = 0.13957018
	rho = pc.rpwaBreitWigner(mPi, mPi, mPi, 1, 1)
	rho.parameters = [0.7690,0.1509]
	masses = []
	for i in range(2000):
		m = 0.278 + 0.001*(i+10)
		if m + mPi > 1.5:
			break
		masses.append(m)
	amps = rho(masses, externalKinematicVariables = [1.5])
	fn = "./rebuiltRho.dat"
	with open(fn, 'w') as out:
		for i,m in enumerate(masses):
			out.write(str(m)  + ' ' + str(amps[i].real) + ' ' + str(amps[i].imag) + '\n')
	fn2 = "/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/anything/zeroModes/bwAmplitudes_noBF/amp_0mp0pRhoPiP"
	compareFiles(fn, fn2)

if __name__ == "__main__":
	main()
