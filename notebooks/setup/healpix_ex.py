#nside = 16, 32, 64, 128...
nside = 128

theta, phi = hp.pix2ang(nside,np.arange(12 * nside * nside))
alt = np.pi - theta
az = phi
J = spline_beam_func(150e6, az, alt)
hp.mollview(J[:,1,0], rot=[0, 90])

alt = np.pi/2 - theta
J = spline_beam_func(150e6, alt, az)

hp.mollview(J[:, 1, 0], rot=[0, 90])
hp.mollview(J[:, 0, 1], rot=[0, 90])
hp.mollview(J[:, 1, 1], rot=[0, 90])
hp.mollview(J[:, 0, 0], rot=[0, 90])

hp.orthview(J[:, 0, 0], rot=[0, 90])
hp.orthview(J[:, 0, 0], rot=[0, 90], half_sky=True)
hp.orthview(np.log10(J[:, 0, 0]), rot=[0, 90], half_sky=True)

hp.orthview(np.abs(J[:, 0, 0]), rot=[0, 90], half_sky=True)
hp.orthview(np.abs(J[:, 1, 0]), rot=[0, 90], half_sky=True)
hp.orthview(np.abs(J[:, 0, 1]), rot=[0, 90], half_sky=True)
hp.orthview(np.abs(J[:, 1, 1]), rot=[0, 90], half_sky=True)

#https://healpix.sourceforge.io/
