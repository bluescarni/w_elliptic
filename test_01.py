def test_01():
	from mpmath import mp, mpc
	import random
	from weierstrass_elliptic import weierstrass_elliptic as we
	random.seed(0)
	mp.dps = 50
	l = [we(random.uniform(-10,20),random.uniform(-10,10)) for _ in range(0,1000)]
	# Extract the tuple of interest.
	lt = [(w.invariants[0],w.invariants[1],w.roots[0],w.roots[1],w.roots[2],w.periods[0],w.periods[1]) for w in l]
	# Fix negative periods and put everything in mpc (apart from invariants).
	lt = [(t[0],t[1],mpc(t[2]),mpc(t[3]),mpc(t[4]),mpc(t[5]),t[6] if t[6].real > 0 else mpc(-t[6].real,t[6].imag)) for t in lt]
	# Convert to suitable string format.
	lt = [(str(t[0]),str(t[1]),str(t[2].real),str(t[2].imag),str(t[3].real),str(t[3].imag),str(t[4].real),str(t[4].imag),str(t[5].real),str(t[5].imag),str(t[6].real),str(t[6].imag)) for t in lt]
	# Build the string.
	return "\n".join([','.join(t) for t in lt])
