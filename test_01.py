# Test for roots, periods, eta.
def test_01():
    from mpmath import mp, mpc
    import random
    from weierstrass_elliptic import weierstrass_elliptic as we
    random.seed(0)
    mp.dps = 50
    l = [we(random.uniform(-10,20),random.uniform(-10,10)) for _ in range(0,1000)]
    # Add a couple of cases wih g2/g3 = 0.
    l.append(we(random.uniform(-10,20),0.))
    l.append(we(random.uniform(-10,20),0.))
    l.append(we(random.uniform(-10,20),0.))
    l.append(we(0.,random.uniform(-10,10)))
    l.append(we(0.,random.uniform(-10,10)))
    l.append(we(0.,random.uniform(-10,10)))
    # Extract the tuple of interest.
    lt = [(w.invariants[0],w.invariants[1],w.roots[0],w.roots[1],w.roots[2],w.periods[0],w.periods[1]) for w in l]
    # Fix negative periods and put everything in mpc (apart from invariants).
    lt = [[t[0],t[1],mpc(t[2]),mpc(t[3]),mpc(t[4]),mpc(t[5]),t[6] if t[6].real > 0 else mpc(-t[6].real,t[6].imag)] for t in lt]
    # Add also eta.
    lt = [t[1] + [t[0].zeta(t[0].periods[0]/2)] for t in zip(l,lt)]
    # Convert to suitable string format.
    lt = [[str(t[0]),str(t[1]),str(t[2].real),str(t[2].imag),str(t[3].real),str(t[3].imag),str(t[4].real),str(t[4].imag),str(t[5].real),str(t[5].imag),
        str(t[6].real),str(t[6].imag),str(t[7].real),str(t[7].imag)] for t in lt]
    # Build the string.
    return "\n".join([','.join(t) for t in lt])

# Test for real P.
def test_02():
    from mpmath import mp, mpc
    import random
    from weierstrass_elliptic import weierstrass_elliptic as we
    random.seed(0)
    mp.dps = 50
    # Generate 20 random objects.
    l = [we(random.uniform(-10,20),random.uniform(-10,10)) for _ in range(0,20)]
    # Add a couple of cases wih g2/g3 = 0.
    l.append(we(random.uniform(-10,20),0.))
    l.append(we(random.uniform(-10,20),0.))
    l.append(we(random.uniform(-10,20),0.))
    l.append(we(0.,random.uniform(-10,10)))
    l.append(we(0.,random.uniform(-10,10)))
    l.append(we(0.,random.uniform(-10,10)))
    retval = []
    for w in l:
        p1,p2 = w.periods
        g2,g3 = w.invariants
        # Generate 100 random values within 10 real periods.
        z_list = [p1.real*random.uniform(-10,10) for _ in range(0,100)]
        values_list = [w.P(_) for _ in z_list]
        res = [(t[0],t[1].real) for t in zip(z_list,values_list)]
        retval.append([str(g2),str(g3)] + [str(item) for sublist in res for item in sublist])
    # Build the string.
    return "\n".join([','.join(t) for t in retval])

# Test for complex P.
def test_03():
    from mpmath import mp, mpc
    import random
    from weierstrass_elliptic import weierstrass_elliptic as we
    random.seed(0)
    mp.dps = 50
    # Generate 20 random objects.
    l = [we(random.uniform(-10,20),random.uniform(-10,10)) for _ in range(0,20)]
    # Add a couple of cases wih g2/g3 = 0.
    l.append(we(random.uniform(-10,20),0.))
    l.append(we(random.uniform(-10,20),0.))
    l.append(we(random.uniform(-10,20),0.))
    l.append(we(0.,random.uniform(-10,10)))
    l.append(we(0.,random.uniform(-10,10)))
    l.append(we(0.,random.uniform(-10,10)))
    retval = []
    for w in l:
        p1,p2 = w.periods
        g2,g3 = w.invariants
        # Generate 100 random values within 10 cells.
        z_list = [p1*random.uniform(-10,10) + p2*random.uniform(-10,10) for _ in range(0,100)]
        values_list = [w.P(_) for _ in z_list]
        res = [(t[0].real,t[0].imag,t[1].real,t[1].imag) for t in zip(z_list,values_list)]
        retval.append([str(g2),str(g3)] + [str(item) for sublist in res for item in sublist])
    # Build the string.
    return "\n".join([','.join(t) for t in retval])

# Test for real P'.
def test_04():
    from mpmath import mp, mpc
    import random
    from weierstrass_elliptic import weierstrass_elliptic as we
    random.seed(0)
    mp.dps = 50
    # Generate 20 random objects.
    l = [we(random.uniform(-10,20),random.uniform(-10,10)) for _ in range(0,20)]
    # Add a couple of cases wih g2/g3 = 0.
    l.append(we(random.uniform(-10,20),0.))
    l.append(we(random.uniform(-10,20),0.))
    l.append(we(random.uniform(-10,20),0.))
    l.append(we(0.,random.uniform(-10,10)))
    l.append(we(0.,random.uniform(-10,10)))
    l.append(we(0.,random.uniform(-10,10)))
    retval = []
    for w in l:
        p1,p2 = w.periods
        g2,g3 = w.invariants
        # Generate 100 random values within 10 real periods.
        z_list = [p1.real*random.uniform(-10,10) for _ in range(0,100)]
        values_list = [w.Pprime(_) for _ in z_list]
        res = [(t[0],t[1].real) for t in zip(z_list,values_list)]
        retval.append([str(g2),str(g3)] + [str(item) for sublist in res for item in sublist])
    # Build the string.
    return "\n".join([','.join(t) for t in retval])

# Test for complex P'.
def test_05():
    from mpmath import mp, mpc
    import random
    from weierstrass_elliptic import weierstrass_elliptic as we
    random.seed(0)
    mp.dps = 50
    # Generate 20 random objects.
    l = [we(random.uniform(-10,20),random.uniform(-10,10)) for _ in range(0,20)]
    # Add a couple of cases wih g2/g3 = 0.
    l.append(we(random.uniform(-10,20),0.))
    l.append(we(random.uniform(-10,20),0.))
    l.append(we(random.uniform(-10,20),0.))
    l.append(we(0.,random.uniform(-10,10)))
    l.append(we(0.,random.uniform(-10,10)))
    l.append(we(0.,random.uniform(-10,10)))
    retval = []
    for w in l:
        p1,p2 = w.periods
        g2,g3 = w.invariants
        # Generate 100 random values within 10 cells.
        z_list = [p1*random.uniform(-10,10) + p2*random.uniform(-10,10) for _ in range(0,100)]
        values_list = [w.Pprime(_) for _ in z_list]
        res = [(t[0].real,t[0].imag,t[1].real,t[1].imag) for t in zip(z_list,values_list)]
        retval.append([str(g2),str(g3)] + [str(item) for sublist in res for item in sublist])
    # Build the string.
    return "\n".join([','.join(t) for t in retval])
