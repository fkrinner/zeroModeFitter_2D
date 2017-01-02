allWaveListFileName = "./wavelist.allWaveList"
waveNameMap = {}
def loadWaveNameMap():
	global waveNameMap
	waveNameMap = {}
	with open(allWaveListFileName,'r') as inin:
		for i,line in enumerate(inin.readlines()):
			waveNameMap[line.strip()] = i

strL = ["S","P","D","F","G", "H", "7"]
class waveName:
	def __init__(self, waveName):
		self.waveName = waveName
		self.J    = self.getJ()
		self.P    = self.getP()
		self.C    = self.getC()
		self.M    = self.getM()
		self.eps  = self.geteps()
		self.isob = self.getIsob()
		self.Jisob, self.Pisob, self.Cisob = self.getJPCisob()
		self.mMin, self.mMax = self.getMassRange()
		self.L    = self.getL()
		self.binned = self.isBinned()
		self.sector = self.getSector()
		self.oldWaveName = self.getOldWaveName()

	def getJ(self):
		J = int(self.waveName[4])
		return J

	def getP(self):
		if self.waveName[5] == "+":
			return +1
		elif self.waveName[5] == "-":
			return -1
		raise ValueError("Unknown parity: '" +  self.waveName[5] + "'")

	def getC(self):
		return +1

	def getM(self):
		return int(self.waveName[7])

	def geteps(self):
		if self.waveName[8] == '+':
			return +1
		elif self.waveName[8] == "-":
			return -1
		raise ValueError("Unknown reflectivity: '" + self.waveName[8] + "'")

	def getL(self):
		return int(self.waveName[-8])

	def getMassRange(self):
		if not self.isBinned():
			return 0.,float('inf')
		else:
			chunks = self.isob[1:-4].split(',')
			return float(chunks[0]), float(chunks[1])

	def getJPC(self):
		rs = str(self.J) #rs = returnString
		if self.P > 0:
			rs += '+'
		else:
			rs += '-'
		if self.C > 0:
			rs += '+'
		else:
			rs += '-'
		return rs
	
	def getJPCMeps(self):
		rs = self.getJPC()
		rs += str(self.M)
		if self.eps > 0:
			rs += '+'
		else:
			rs += '-'
		return rs

	def getIsob(self):
		if "rho_770_0" in self.waveName:
			return "rho"
		if "f2_1270_0" in self.waveName:
			return "f2"
		if "rho3_1690_0" in self.waveName:
			return "rho3"
		if "sigma0" in self.waveName:
			return "f0_500"
		if "f0_980_0" in self.waveName:
			return "f0_980"
		if "f0_1500_0" in self.waveName:
			return "f0_1500"
		if not self.isBinned():
			raise valueError("Unknwon isobar for wave '"+self.waveName+"'")
		return self.getBinnedIsob()

	def isBinned(self):
		if "binned" in self.waveName:
			return True
		return False

	def getBinnedIsob(self, useMasses = True):
		chunks = self.waveName.split("binned[")[1].split(']')[0].split(',')
		if useMasses:
			return '['+chunks[2]+','+chunks[3]+']'+chunks[1]
		else:
			return '[pi,pi]'+chunks[1]

	def getJPCisob(self):
		if self.isob == "rho" or "1--" in self.isob:
			return 1,-1,-1
		if self.isob == "f2" or "2++" in self.isob:
			return 2,1,1
		if self.isob == "rho3" or "3--" in self.isob:
			return 3,-1,-1
		if self.isob.startswith("f0_") or "0++" in self.isob:
			return 0,1,1
		raise ValueError("Unkknown JPC for isobar: '"+self.isob+"'")

	def getSector(self):
		if not self.isBinned():
			return self.getOldWaveName()
		else:
			return self.getJPCMeps()+self.getBinnedIsob(useMasses = False)+"Pi"+strL[self.L]

	def getOldWaveName(self):
		wn = self.getJPCMeps()
		wn += self.isob
		wn += "Pi"
		wn += strL[self.L]
		return wn

	def getWaveIndex(self):
		if len(waveNameMap) == 0:
			loadWaveNameMap()
		if not self.waveName in waveNameMap:
			raise ValueError("'" + self.waveName + "' not in waveNameMap")
		return waveNameMap[self.waveName]

	def __str__(self):
		return self.oldWaveName
