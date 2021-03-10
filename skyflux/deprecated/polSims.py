#! /usr/bin/env python

""" Begin section: pared-down copy of
    Chuneeta/PolarizedSims/COSMO_constants.py """

import numpy as np
import astropy.cosmology as CS

C = 2.99e8 # SPEED OF LIGHT IN M/S
F21 = 1.42040575177 # FREQUENCY OF 21 CM HYDROGEN LINE IN GHZ
COSMO = CS.FlatLambdaCDM(H0=100.0, Om0=0.27) # H0 = 100 km/s/Mpc

def fq2z(fq):
   """
   Redshift corresponding the specified frequency

   Input(s)
      fq :  [scalar] frequency in Hz
   """
   return F21 / fq - 1
   
def transverse_comoving_distance(z):
   """
   Transverse comoving distance at redshift z
       corresponding to an
       angular separation of 1 radian in Mpc/h.
   
   Input(s)
      z :  [scalar] redshift
   """
   return COSMO.comoving_distance(z).value # Mpc/h 

def comoving_depth(B, z):
   """
   Comoving line-of-sight depth corresponding to
       specified redshift and bandwidth for redshifted
       21 cm line in Mpc/h
  
   Input(s)
      B :    [scalar] Observing bandwith in Hz
      z :    [scalar] redshift
   """
   return (C / 1e3) * B * (1 + z)**2 / \
       F21/COSMO.H0.value/COSMO.efunc(z) # Mpc/h 

def dkprll_deta(z):
   """
   Constant to transform delays to line-of-sight
       wavenumbers corresponding to redshift and
       21 CM HI line in h/Mpc
   
   Input(s)
      z :  [scalar] redshift
   """
   return 2 * np.pi * COSMO.H0.value * F21 * \
       COSMO.efunc(z) / C /  (1 + z)**2 * 1e3

def k_parallel(delays, z):
   """
   Compute line-of-sight wavenumbers corresponding to
       specified delays and redshift for
       redshifted 21 cm line in h/Mpc

   Input(s):
      z : [scalar] redshift
   """
   return dkprll_deta(z) * delays

def k_perp(z):
   """
   Compute transverse wavenumbers corresponding to
       redshifted 21 cm line in h/Mpc

   Input(s)
      z              : [scalar] redshift
   """
   return 2 * np.pi / transverse_comoving_distance(z)

""" End section: Chuneeta/PolarizedSims/COSMO_constants.py """

""" Begin section: pared-down copy of
    Chuneeta/PolarizedSims/genPowerSpectra.py"""

def genWindow(size):
    """
    Implements Blackmann-Harris filter

    size : Size/Lenth of frequency channel to which
        filter is applied; type:int
    """
    window = np.zeros((size),)
    alpha = 0.16 #! magic number
    _w = lambda n: \
        (1-alpha) / 2. - 0.5 * np.cos((2 * np.pi * n)/(size - 1)) \
        + alpha / 2. * np.cos((4 * np.pi * n)/(size - 1))
    window[:] = _w(np.arange(size))
    return window

def f2etas(freqs):
    """
    Evaluates geometric delay (fourier conjugate of frequency)
   
    -freqs: Frequencies in GHz; type:numpy.ndarray
    """
    df = freqs[1] - freqs[0]
    etas = np.fft.fftfreq(freqs.size, df)
    return etas

def delay_transform(data,fqs,convert=None):
    """
    Fourier transforms visibility along frequency axis

    - data: per baseline visibility; type:numpy.ndarray
    - fqs:  slected frequencies in GHz; dtypw:numpy.ndarray
    """
    N = fqs.size
    df = fqs[1] - fqs[0]
    window = genWindow(N) #! this label is never used
    delaySpec = np.fft.ifft(data) * N * df
    return delaySpec 
      
""" End section: Chuneeta/PolarizedSims/genPowerSpectra.py """

""" Remainder of this file: the code unique to this script.
    i.e. my own work. "power plot sketch" """
