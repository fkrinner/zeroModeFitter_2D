from utils import countSectors
##
#
# Post process the zero modes
#
##
def removeCertainZeroModes(zeroHists, eigenvalues):
	if len(zeroHists) == 0.:
		return
	if not len(eigenvalues) == len(zeroHists):
		raise ValueError("removeCertainZeroModes(...): Sizes do not match")
#	removeAllButOne2mpModes(zeroHists, eigenvalues)
#	setExplicitelyToZero(zeroHists, ['zero2_0']) # in the MC results, 'zero3_0', 'zero2_0' and 'zero4_0' are 2mp zero-modes

def removeAllButOne2mpModes(zeroHists, eigenvalues):
	"""
	It is knonw, that for the standard four 2mp waves, only one zero-mode 
	exists, while a second small EV in not a zero mode. Therefore remove all
	but one mode and throw an exception if more than four 2mp waves appear
	"""
	zmpModes = []	
	for i_hist, hist in enumerate(zeroHists):
		if countSectors(hist) > 4:
			raise NotImplementedError("More than four 2mp waves found")
		if "2-+0+" in hist.GetTitle():
			zmpModes.append(i_hist)
	nBins = zeroHists[0].GetNbinsX()
	for bin in range(nBins):
		evList = []
		for i_hist in zmpModes:
			ev = eigenvalues[i_hist].GetBinContent(bin+1)
			evList.append((ev, i_hist))
		evList.sort()
		evList.reverse()
		skipped = False
		for pp in evList:
			if pp[0] == 0.:
				continue
			if not skipped:
				isZero = True
				for binY in range(zeroHists[pp[1]].GetNbinsY()):
					if not zeroHists[pp[1]].GetBinContent(bin+1, binY+1) == 0.:
						isZero = False
						break
				if not isZero:
					skipped = True
			else:
				for binY in range(zeroHists[pp[1]].GetNbinsY()):
					zeroHists[pp[1]].SetBinContent(bin+1, binY+1, 0.)
#	for bin in range(nBins):
#		countnonZer = 0
#		for i in zmpModes:
#			isZero = True
#			for bb in range(zeroHists[i].GetNbinsY()):
#				if not zeroHists[i].GetBinContent(bin+1, bb+1) == 0.:
#					isZero = False
#					break
#			if not isZero:
#				countnonZer += 1
#		print bin, countnonZer

def setExplicitelyToZero(zeroHists, histNames):
	for hist in zeroHists:
		if hist.GetName() in histNames:
			print "Removing: '" + hist.GetName() + "'"
			for i in range(hist.GetNbinsX()):
				for j in range(hist.GetNbinsY()):
					hist.SetBinContent(i+1, j+1, 0.)



