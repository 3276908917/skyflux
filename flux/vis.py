def visibility(ant1, ant2, source, nu=151e6):
    """
    Visibility integrand evaluated for a single source.

    The most glaring waste of compute is separately calculating
    the values for RA and DEC, although it is not clear to me
    how many clock cycles are actually spent thereon.
    """
    I = source.flux_by_frq[nu / 1e6]
    s = np.array([complex(I), 0, 0, 0])

    ra = np.radians(source.ra_angle)
    dec = np.radians(source.dec_angle)
    r = rot.raddec2lm(ra, dec)

    phi = phase_factor(ant1, ant2, r, nu)
    return np.dot(np.dot(stokes.A_matrix(ra, dec, nu), s), phi)

def visibility_integrand(ant1, ant2, nu=151e6):
    """
    Return visibility integrand evaluated for all sources
    in obj_catalog.
    """
    total = np.array([0j, 0j, 0j, 0j]) # 4 x 1. Visibility has a phase,
    for source in obj_catalog:
        total += visibility(ant1, ant2, source, nu)
    return total
