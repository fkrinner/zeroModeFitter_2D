from cmath import log, pi

def breakupMomentumSquared(M, m1, m2):
	mSumm = m1 + m2
	mDiff = m1 - m2
	return (M - mSumm) * (M + mSumm) * (M - mDiff) * (M + mDiff) / (4*M**2)

def breakupMomentum(M, m1, m2):
	q2 = breakupMomentumSquared(M, m1, m2)
	if q2 < 0.:
		return 0.
	return q2**.5

def barrierFactorSquared(L, q, Pr):
	"""
	Barrier factor squared, with L as the angular momentum, q as the 
	breakup momentum and Pr as the momentum corresponding to 
	1/<interactionRadius> as taken from ROOTPWA, with the exception, that
	L is NOT multipliet by 2
	"""
	z = q**2/Pr**2
	if L == 0:
		return 1
	if L == 1:
		return (2 * z) / (z + 1)
	if L == 2:
		return (13 * z**2) / (z * (z + 3) + 9)
	if L == 3:
		return (277 * z**3) / (z * (z * (z + 6) + 45) + 225)
	if L == 4:
		return (12746 * z**4) / (z * (z * (z * (z + 10) + 135) + 1575) + 11025)
	if L == 5:
		return (998881 * z**5) / (z * (z * (z * (z * (z + 15) + 315) + 6300) + 99225) + 893025)
	if L == 6:
		return (118394977 * z**6) / (z * (z * (z * (z * (z * (z + 21) + 630) + 18900) + 496125) + 9823275) + 108056025)
	if L == 7:
		return (19727003738 * z**7)/ (z * (z * (z * (z * (z * (z * (z + 28) + 1134) + 47250) + 1819125) + 58939650) + 1404728325L) + 18261468225)
	raise ValueError("Barrier factor not implemetned for L > 7")

def barrierFactor(L, q, Pr):
	bf2 = barrierFactorSquared(L, q, Pr)
	if bf2 < 0.:
		return 0.
	return bf2**.5

def breitWigner(m, m0, G0, L, q, q0, Pr):
	if q0 == 0.:
		return 0.+0.j
	G = G0 * (m0/m) * (q/q0) * barrierFactorSquared(L, q, Pr)/barrierFactorSquared(L,q0, Pr)
	return m0 * G0 / ( m0**2 - m**2 - 1.j*m0*G)

def flatte(M,M0,g1,g2,m1,m2,m2nd1,m2nd2):
	phi1 = 2*breakupMomentum(M,m1,m2)/M
	if phi1 == 0.:
		return 0.+0.j
	phi2 = 2*breakupMomentum(M,m2nd1,m2nd2)/M
	return 1./(M0**2 - M**2 -1.j*(phi1*g1 + phi2*g2))

def threshChewMandelstam(s,sThresh):
	s       = s+0.j
	sThresh = sThresh+0.j
	retVal  = -1./s*(sThresh-s)**.5*(-s)**.5
	retVal *= log(( (sThresh-s)**.5 + (-s)**.5)/(sThresh**.5))
	retVal -= 3./2.
	retVal *= -2./pi
	return retVal

def sameMassChewMandelstam(s,m):
	return threshChewMandelstam(s,4*m**2)

def chewMandelstam(s, m1, m2):
	if m1 == m2:
		return sameMassChewMandelstam(s,m1)
	s  = s +0.j
	m1 = m1+0.j
	m2 = m2+0.j
	retVal  = -1./s*((m1+m2)**2-s)**.5 * ((m1-m2)**2-s)**.5
	retValBad = retVal
	retVal *=  log((((m1+m2)**2-s)**.5 + ((m1-m2)**2-s)**.5)/(2*(m1*m2)**.5))
	retVal += (m1**2-m2**2)/(2*s)*log(m1/m2)
	retVal -= (m1**2+m2**2)/(2*(m1**2-m2**2))*log(m1/m2)
	retVal -= .5
	retVal *= -2./pi
	return retVal

def rho(s,m1,m2):
	"""
	Phase space volume [rho] = 1
	rho(M**2, m1, m2) = q(M,m1,m2)*2/M | q = breakupMomentum
	"""
	return (s**2 + m1**4 + m2**4 - 2*s*m1**2 - 2*s*m2**2 - 2*m1**2*m2**2)**.5/s

def main():
	m1 = 0.139
	m2 = 0.239
	for i in range(100):
		s = 1. + 0.01*i
		print chewMandelstam(s,m1,m2).imag,rho(s,m1,m2),breakupMomentum(s**.5, m1, m2)*2/s**.5

if __name__ == "__main__":
	main()
	
