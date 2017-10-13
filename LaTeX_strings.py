
def getFreeSubscript(sector):
	if sector == "0++":
#		return r"\text{S}"
		return r"0^{++}"
	if sector == "1--":
#		return r"\text{P}"
		return r"1^{--}"
	return sector[0] + r"^{" + sector[1:3] + r"}"


def getProperWaveName(sector):
	if sector == "Dp[pi,pi]0++PiS":
		return r"S wave"
	if sector == "Dp[pi,pi]1--PiP":
		return r"P wave"
	retVal = r"$" + sector[0] + r"^{" + sector[1:3] + r"}" + sector[3] + r"^" + sector[4] + r"[\pi\pi]_{" + getFreeSubscript(sector[12:15])+ r"}\pi$" + sector[-1] #+ r" wave"
	return retVal

#tBins = [r"$t^\prime = 0.100$--$0.141\,($GeV$/c)^2$",r"$t^\prime = 0.141$--$0.194\,($GeV$/c)^2$",r"$t^\prime = 0.194$--$0.326\,($GeV$/c)^2$",r"$t^\prime = 0.326$--$1.000\,($GeV$/c)^2$"]
tBins = [r"$0.100<t^\prime<0.141\,(\text{GeV}/c)^2$",r"$0.141<t^\prime<0.194\,(\text{GeV}/c)^2$",r"$0.194<t^\prime<0.326\,(\text{GeV}/c)^2$",r"$0.326<t^\prime<1.000\,(\text{GeV}/c)^2$"]

fullTrange = r"$0.100<t^\prime<1.000\,($GeV$/c)^2$"
def getProperDataSet(fileName, tBin):
	if "extensiveFreedIsobarStudies" in fileName:
		if "MC" in fileName:
			return r"Monte Carlo"
		return tBins[tBin]
	return ""

def getMassString(mMin, mMax):
#	retVal = r"$m_{3\pi} = " + "{0:.2f}".format(mMin) + r"$--$" + "{0:.2f}".format(mMax) + r"$GeV$/c$"
	retVal ="$" + "{0:.2f}".format(mMin) + r"<m_{3\pi}<" + "{0:.2f}".format(mMax) + r"\text{GeV}/c^2$"
	return retVal

m3Pi = r"$m_{3\pi}\,[\text{GeV}/c^2]$"
#m2Pi = r"$m_{\pi^-\pi^+}\,[\text{GeV}/c^2]$"
m2Pi = r"$m_{2\pi}\,[\text{GeV}]$"
#m2Pi = r"$m_{\pi^+\pi^-}\,[$GeV$/c]$"

#intens = r"Intensity [A.U.]"
#intens = r"Intensity [Events/(GeV/$c^2$)]"
intens = r"Intensity [events/40\,MeV]"

real = r"$\Re\,[\sqrt{\text{events/40\,MeV}}]$"
imag = r"$\Im\,[\sqrt{\text{events/40\,MeV}}]$"

#real = r"Re$(\mathcal{T}_\text{bin})\,[(\text{Events}/(\text{GeV}/c^2))^{1/2}]$"
#imag = r"Im$(\mathcal{T}_\text{bin})\,[(\text{Events}/(\text{GeV}/c^2))^{1/2}]$"


# self.overrideMassString = r"$m_{3\pi} = 1.87$GeV$c$"
# self.titleRight = r"S wave"
# self.titleRight = r"P wave"
