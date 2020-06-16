import numpy as np
import time
import astropy
import astropy.time

"""
S ~/prop\nu^alpha
S_{200} = S_{150} * (200/150)^\alpha
"""

def collapse_angle(degree, arcminute=0, arcsecond=0):
    return degree + arcminute / 60 + arcsecond / 3600

def collapse_hour(hour, minute=0, second=0):
    return 15 * hour + minute / 4 + second / 240

hera_lat = -collapse_angle(30, 43, 17)
hera_lon = collapse_angle(21, 25, 42)

def get_lst(lon = hera_lon):
    """
    Return current local sidereal time (LST)
    for longitude @lon (default is HERA array).
    LST is returned in radians.
    """
    t = astropy.time.Time(time.time(), format='unix')
    return t.sidereal_time('apparent', longitude=lon).radian

# The change-of-basis matrix between equatorial and galactic coordinate systems
M_eq_to_gal = np.array([
    [-.054876, -.873437, -.483835],
    [.494109, -.444830, .746982],
    [-.867666, -.198076, .455984]
])

def ha_to_gal(ha, dec, lst, radians=False):
    if not radians:
        ha = np.radians(ha)
        dec = np.radians(dec)
        lst = np.radians(lst)
    rct = rectangle(ha, dec)
    ra_dec = np.dot(np.linalg.inv(M_eq_to_ha(lst)), rct)
    gal = np.dot(M_eq_to_gal, ra_dec)
    return new_sphere(gal, radians)

def gal_to_eq(el, be, radians=False):
    if not radians:
        el = np.radians(el)
        be = np.radians(be)
    rct = rectangle(el, be)
    ra_dec = np.dot(np.linalg.inv(M_eq_to_gal), rct)
    return new_sphere(ra_dec, radians)

def eq_to_gal(ra, dec, radians=False):
    '''
    @radians determines the format of BOTH input and output!
    Given a pair of angles @ra and @dec,
    return a pair of angles relating the associated
    galactic longitude (first?) and latitude (second?).
    '''
    if not radians:
        ra = np.radians(ra)
        dec = np.radians(dec)
    eq_vector = rectangle(ra, dec)
    gal_vector = np.dot(M_eq_to_gal, eq_vector)
    return new_sphere(gal_vector, radians)

def M_eq_to_ha(LST):
    '''
    Return the change-of-basis matrix between the equatorial and
    hour angle declination coordinate systems.
    The conversion depends on the @LST, Local Siderial Time
    '''
    s = np.sin(LST)
    c = np.cos(LST)
    return np.array([[c, s, 0], [s, -c, 0], [0, 0, 1]])

def M_ha_to_topo(phi):
    '''
    Return the change-of-basis matrix between the hour angle declination
    and topocentric coordinate systems.
    The conversion depends on the user's current latitude @phi,
        which must be given in radians.
    '''
    s = np.sin(phi)
    c = np.cos(phi)
    return np.array([[-s, 0, c], [0, -1, 0], [c, 0, s]])

def rectangle(a, b):
    '''
    Given a pair of angles (both angles must be in radians),
    return the corresponding 3x1 rectangular vector.
    '''
    return np.array([np.cos(b) * np.cos(a), np.cos(b) * np.sin(a), np.sin(b)])

def gal_to_topo(el, be, jd, lat, lon, radians=False):
    '''
    @radians determines the format of BOTH input and output!
    Given a pair of angles @el and @be (in galactic coordinates),
    return a pair of angles relating the associated
    azimuth and altitude.

    This does not work in the current environment
    we are using ugradio as a crutch.
    '''
    if not radians:
        el = np.radians(el)
        be = np.radians(be)
        lat = np.radians(lat)
    else:
        lon = np.degrees(lon)
    rct = rectangle(l, b)
    ra_dec = np.dot(np.linalg.inv(M_eq_to_gal), rct)
    lst = get_lst(lon)
    hrd = np.dot(np.linalg.inv(M_eq_to_ha(lst)), ra_dec)
    topo = np.dot(M_ha_to_topo(phi), hrd)
    return new_sphere(topo, radians)

def new_sphere(out_arr, radians=False):
    '''
    Given a 3x1 vector,
    return the corresponding pair of angles
    @radians determines whether the angles are given in radians.
    '''
    gp = np.arctan2(out_arr[1], out_arr[0])
    tp = np.arcsin(out_arr[2])
    if not radians:
        return np.degrees(gp), np.degrees(tp)   
    return gp, tp

def ha_to_topo(ha, dec, lat, radians=False):
    '''
    Take a position in hour-angle right ascension / declination
        to local altitude and azimuth.
    This performs NO precession.
    '''
    if not radians:
        ha = np.radians(ha)
        dec = np.radians(dec)
        lat = np.radians(lat)
    rct = rectangle(ha, dec)
    topo = np.dot(M_ha_to_topo(lat), rct)
    return new_sphere(topo, radians)

def ha_to_eq(ha, dec, lat, radians=False):
    '''
    Take a position in hour-angle right-ascension / declination
        to regular right-ascension / declination.
    '''
    if not radians:
        ha = np.radians(ha)
        dec = np.radians(dec)
        lat = np.radians(lat)
    rct = rectangle(ha, dec)
    eq = np.dot(np.linalg.inv(M_eq_to_ha(lat)), rct)
    return new_sphere(eq, radians)

def eq_to_topo(ra, dec, latitude, lst, radians=False):
    '''
    @radians determines the format of BOTH input and output!
    Given a pair of angles @ra and @dec,
    return a pair of angles relating the associated
    azimuth (first) and altitude (second).
    '''
    if not radians:
        ra = np.radians(ra)
        dec = np.radians(dec)
        latitude = np.radians(latitude)
        lst = np.radians(lst)
    eq_vector = rectangle(ra, dec)
    ha_vector = np.dot(M_eq_to_ha(lst), eq_vector)
    topo_vector = np.dot(M_ha_to_topo(latitude), ha_vector)
    return new_sphere(topo_vector, radians)

"""
The following function was written by C. D. Nunhokee,
'genVisibility.py', polarizedSims, Feb 8 2019
https://github.com/Chuneeta/polarizedSims/blob/master/genVisibility.py
"""
# Unfortunately, the parameter order has been switched wrt the citation.
def raddec2lm(ra, dec, ra0=None, dec0=hera_lat): # ra and dec in radians
    """
    Converts ra/dec to direction cosines l/m
    ra0  : reference/phase right ascension; type: float
    dec0 : reference/phase declination; type:float
    ra   : right ascension in radians; type:float
    dec  : declination in radians; type:float
    """
    # See note at the end about default arguments.
    if ra0 is None:
        ra0 = get_lst()

    l = np.cos(dec) * np.sin(ra0 - ra)
    m = -1 * (np.sin(dec) * np.cos(dec0) - \
        np.cos(dec) * np.sin(dec0) * np.cos(ra-ra0))
    return l, m
