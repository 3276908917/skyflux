#! /usr/bin/env python

import numpy as np
import astropy.cosmology as CS

C = 2.99e8 # SPEED OF LIGHT IN M/S
F21 = 1.42040575177 # FREQUENCY OF 21 CM HYDROGEN LINE IN GHZ
COSMO = CS.FlatLambdaCDM(H0=100.0, Om0=0.27)  # Using H0 = 100 km/s/Mpc

def fq2z(fq):
   '''
   Redshift corresponding the specified frequency

   Input(s)
      fq :  [scalar] frequency in Hz
   '''
   z = F21/fq-1
   return z

def z2fq(z):
   '''
   Frequency (in Hz) corresponding the specified redshift

   Input(s)
      z :  [scalar] redshift 
   '''
   fq = F21/(z+1)
   return fq

   
def transverse_comoving_distance(z):
   '''
   Transverse comoving distance at redshift z corresponding to an angular separation of 1 radian in Mpc/h 

   Input(s)
      z :  [scalar] redshift
 
   '''
   Dz =  COSMO.comoving_distance(z).value # Mpc/h 
   return Dz


def comoving_depth(B,z):
   '''
   Comoving line-of-sight depth corresponding to specified  redshift and bandwidth for redshifted
   21 cm line in Mpc/h
  
   Input(s)
      B :    [scalar] Observing bandwith in Hz

      z :    [scalar] redshift
   '''
   deltaD =  (C/1e3) * B * (1+z)**2 /F21/COSMO.H0.value/COSMO.efunc(z) # Mpc/h 
   return deltaD


def dkprll_deta(z):
   '''
   Constant to transform delays to line-of-sight wavenumbers corresponding to redshift and 21 CM HI line
   in h/Mpc
   
   Input(s)
      z :  [scalar] redshift
   '''
   return 2 * np.pi * COSMO.H0.value * F21 * COSMO.efunc(z) / C /  (1+z)**2 * 1e3

def k_parallel(delays, z):
   '''
   Compute line-of-sight wavenumbers corresponding to specified delays and redshift for redshifted 21 cm line in h/Mpc

   Input(s):
      z : [scalar] redshift
   '''
   return dkprll_deta(z) * delays

def k_perp(z):
   '''
   Compute transverse wavenumbers corresponding for redshifted 21 cm line in h/Mpc

   Input(s)
      z              : [scalar] redshift
   '''
   kperp = 2 * np.pi / transverse_comoving_distance(z)
   return kperp

def horizon_limit(z):
   '''
   Compute the horizon limit in Mpc/h at a given redshift for a specified baseline length (from Thyagarajan 2013)

   Input(s):
      baseline_length : [scalar] baseline length in wavelengths

      z : [scalar] redshift
   '''
   hlimit = COSMO.H0.value * COSMO.efunc(z) * transverse_comoving_distance(z) / C / (1+z) * 1e3
   return hlimit

#! /usr/bin/env python

"""
  Generating power spectra using delay-filtering approach
"""

import numpy as np
import optparse
import scipy.constants as CNST1
import os,sys
import aipy as a

C = 2.99e8 # SPEED OF LIGHT IN M/S
F21 = 1420405751.77 # FREQUENCY OF 21 CM HYDROGEN LINE
curtime, zen = None, None

#===============================================================
def genWindow(size):
    """
    Implements Blackmann-Harris filter

    size : Size/Lenth of frequency channel to which filter is applied; type:int
    """
    
    window = np.zeros((size),)
    alpha = 0.16
    _w = lambda n: (1-alpha)/2. - 0.5*np.cos((2*np.pi*n)/(size-1)) + alpha/2.*np.cos((4*np.pi*n)/(size-1))
    window[:] = _w(np.arange(size))
    return window

#===============================================================
def frequency_channels(freqs,f0,B):
    """
    Evaluates frequency boundaries given the center frequency and bandwidth

    freqs: Frequencies in GHz; type:numpy.ndarray
    f0   : Center frequency in GHz; type:float
    B    : Bandwidth in GHz; type:float
    """
    ch0 = np.argmin(np.abs(freqs-f0))
    df = freqs[1]-freqs[0]
    dCH = int(B/df)
    ch1, ch2 =  ch0-dCH/2, ch0+dCH/2 
    return ch1, ch2

#===============================================================
def f2etas(freqs):
    """
    Evaluates geometric delay (fourier conjugate of frequency)
   
    -freqs: Frequencies in GHz; type:numpy.ndarray
    """
    df = freqs[1] - freqs[0]
    etas = np.fft.fftfreq(freqs.size, df)
    return etas

#=============================================================== 
def delay_transform(data,fqs,convert=None):
    """
    Fourier transforms visibility along frequency axis

    - data: per baseline visibility; type:numpy.ndarray
    - fqs:  slected frequencies in GHz; dtypw:numpy.ndarray
    """
    N = fqs.size
    df = fqs[1] - fqs[0]
    window = genWindow(N)
    delaySpec = np.fft.ifft(data) * N * df
    return delaySpec 

#===============================================================
def compute_omegaB(beampath,jd,freqs,freq_wgts,nside=None,healpix=None):
    """
    Evaluates 3D vloume normlaization from Thyagarajan 2016
    
    - beampath : Path where the beam are stored; type:str
    - jd       : Julian date at which the beam is generated
    - freqs    : Frequency at which the beam is generated
    - freq_wgts: Weights applied to the visibility
    - healpix  : Enable healpix normalization
    """
    A0 = np.load(beampath+'/jones-f%.4g_j%.5f.npz'%(freqs[0]*1e3,np.float(t)))['mueller']
    df = freqs[1] - freqs[0] # frequency resolution in Hz
    A = np.zeros((4,A0.shape[-1],len(freqs)),dtype=complex) # initializing A marix
    for ii, f in enumerate(freqs):
      for p in range(4):
         jones = np.load(beampath+'/jones-f%.4g_j%.5f.npz'%(f*1e3,np.float(jd)))['mueller']
         A[p,:,ii] = jones[p,p,:]
    if healpix:
       domega = (4*np.pi)/(12*nside*nside)
    else:
       domega = 1.0/A0.shape[-1] # dividing by the number of sources
    Aw = A*freq_wgts
    omegaB = {}
    omegaB = np.nansum(np.nansum(Aw**2,axis=1),axis=1) * domega 
    OmegaB['xx'] = omegaB[0] 
    OmegaB['xy'] = omegaB[1]
    OmegaB['yx'] = omegaB[2]
    OmegaB['yy'] = omegaB[3]
    return omegaB     
      































"power plot sketch"

import skyflux as sf

import matplotlib.pyplot as plt
import numpy as np

import pickle

def load_sim(fname):
    global k_par
    global k_starter
    global nu_idxs

    sim_file = open(fname, "rb")

    meta = pickle.load(sim_file)

    frq = meta['frequencies'] #! I don't remember where this
       # naming scheme is in the code
    etas = f2etas(frq)
    
    # we are ranging from 50 to 250 MHz
    center_f = np.average(frq)
    z = fq2z(center_f)
    lambda_ = C / center_f

    k_par = k_parallel(etas, z)
    k_starter = k_perp(z) / lambda_ # this will need to be
    # multiplied on a per-baseline basis

    nu_idxs = range(len(frq))

    return {"fs": nu_idxs,
            "kp": k_par,
            "ks": k_starter,
            "sim": meta['picture']}

def static_visual(sim_dict):
    wedge_data = []
    
    ### aliasing ###
    nu_idxs = sim_dict['fs']
    k_par = sim_dict['kp']
    k_starter = sim_dict['ks']
    sim = sim_dict['sim']

    for ant1 in sim.keys():
        for ant2 in sim[ant1].keys():
            """
            Since we only used one LST,
            we do not need to do any averaging
                (averaging is supposed to happen over LSTs,
                    not frequency)
            """
            k_orth = k_starter * sf.ant.baselength(ant1, ant2)
            for nu_idx in nu_idxs:
                # this is a proportionality.
                # The real deal uses the power equation 6
                    # from Nunhokee et al.
                    
                brightness = sim[ant1][ant2][nu_idx]
                    
                power_prop = np.log10(np.vdot(
                    brightness,
                    brightness
                ))
                
                wedge_datum = np.array([
                    k_orth,
                    k_par[nu_idx],
                    float(power_prop)
                ])

                wedge_data.append(wedge_datum)
                
    return np.array(wedge_data)
    
def dynamic_visual(sim):
    wedge_data = []

    for ant1 in sim.keys():
        for ant2 in sim[ant1].keys():
            k_orth = k_starter * sf.ant.baselength(ant1, ant2)
            
            for nu_idx in nu_idxs:
                system = sim[ant1][ant2][nu_idx]
                powers_prop = []
                
                for t_idx in range(len(system) - 1):
                    this_instant = system[t_idx]
                    next_instant = system[t_idx + 1]
                    # this is a proportionality.
                    # The real deal uses the power equation 6
                    # from Nunhokee et al.
                    powers_prop.append(np.vdot(
                        this_instant,
                        next_instant
                    ))

                wedge_datum = np.array([
                    k_orth,
                    k_par[nu_idx],
                    float(
                        np.log10(
                            np.average(
                                np.array(powers_prop)
                            )
                        )
                    )
                ])

                wedge_data.append(wedge_datum)
                
    return np.array(wedge_data)

def open_visual(wedge):
    k_orth = wedge[:, 0]
    k_parr = wedge[:, 1]
    p_p = wedge[:, 2]

    scaled_pow = (p_p - p_p.min()) / p_p.ptp()
    colors = plt.cm.viridis(scaled_pow)

    plt.scatter(k_orth, k_parr, marker='.', c=colors)
    plt.show()

