
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
