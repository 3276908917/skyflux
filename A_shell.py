import skyflux as sf
import matplotlib.pyplot as plt
import numpy as np
import healpy as hp

nside = 128
theta, phi = hp.pix2ang(nside, np.arange(12 * nside * nside))
az = phi
alt = np.pi/2 - theta

J = sf.stokes.create_J(az=az, alt=alt)

A = np.array([sf.stokes.create_A(J=Ji) for Ji in J])



def project_A(A, func):
    def orth(i, j, panel, ttl=None):
        if ttl is None:
            ttl = str(i) + ', ' + str(j)
        hp.orthview(func(A[:, i, j]), rot=[0, 90], half_sky=True, sub=[5, 4, panel], title=ttl)
    
    orth(0, 0, 1, 'I\' <- I')
    orth(0, 1, 2)
    orth(0, 2, 3)
    orth(0, 3, 4)
    
    orth(1, 0, 5)
    orth(1, 1, 6, 'Q\' <- Q')
    orth(1, 2, 7)
    orth(1, 3, 8)
    
    orth(2, 0, 9)
    orth(2, 1, 10)
    orth(2, 2, 11, 'U\' <- U')
    orth(2, 3, 12)
    
    orth(3, 0, 13)
    orth(3, 1, 14)
    orth(3, 2, 15)
    orth(3, 3, 16, 'V\' <- V')
    
# project_A(A); plt.show()
