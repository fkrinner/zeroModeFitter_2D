def getProperWaveName(sector):
	retVal = r"$" + sector[0] + r"^{" + sector[1:3] + r"}" + sector[3] + r"^" + sector[4] + r"[\pi\pi]_{" + sector[12] + r"^{" + sector[13:15] + r"}}\pi " + sector[-1] + r"$"
	return retVal

def getProperDataSet(fileName):
	if "extensiveFreedIsobarStudies" in fileName:
		if "MC" in fileName:
			return r"Monte Carlo"
		return r"\textsc{Compass} 2008 "
	return ""

m3Pi = r"$m_{3\pi}$"
m2Pi = r"$m_{\pi^+\pi^-}\,[$GeV$/c]$"
