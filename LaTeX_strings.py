def getProperWaveName(sector):
	retVal = r"$" + sector[0] + r"^{" + sector[1:3] + r"}" + sector[3] + r"^" + sector[4] + r"[\pi\pi]_{" + sector[12] + r"^{" + sector[13:15] + r"}}\pi " + sector[-1] + r"$"
	return retVal

tBins = [r"$t^\prime = 0.100$--$0.141\,($GeV$/c)^2$",r"$t^\prime = 0.141$--$0.194\,($GeV$/c)^2$",r"$t^\prime = 0.194$--$0.326\,($GeV$/c)^2$",r"$t^\prime = 0.326$--$1.000\,($GeV$/c)^2$"]

def getProperDataSet(fileName, tBin):
	if "extensiveFreedIsobarStudies" in fileName:
		if "MC" in fileName:
			return r"Monte Carlo"
		return tBins[tBin]
	return ""

def getMassString(mMin, mMax):
	retVal = r"$m_{3\pi} = " + "{0:.2f}".format(mMin) + r"$--$" + "{0:.2f}".format(mMax) + r"$GeV$/c$"
	return retVal

m3Pi = r"$m_{3\pi}\,[$GeV$/c]$"
m2Pi = r"$m_{\pi^+\pi^-}\,[$GeV$/c]$"

#intens = r"Intensity [A.U.]"
intens = r"Intensity [number of events per bin]"

real = r"$\Re$"
imag = r"$\Im$"
