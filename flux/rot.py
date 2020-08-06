import numpy as np
import time
import astropy.time
from astropy import units as u

def collapse_angle(degree, arcminute=0, arcsecond=0, radians=False):
    """
    Return a single angle, based on an
    angle in the form of degrees, arcminutes, and arcseconds.

    @radians
        True: return the angle in radians
        False: in degrees
    """
    deg = degree + arcminute / 60 + arcsecond / 3600
    if radians:
        return np.radians(deg)
    return deg

def collapse_hour(hour, minute=0, second=0, radians=False):
    """
    Return a single angle, based on an
    angle in the form of hours, minutes, and seconds of a day.

    @radians
        True: return the angle in radians
        False: in degrees
    """
    deg = 15 * hour + minute / 4 + second / 240
    if radians:
        return np.radians(deg)
    return deg

hera_lat = -collapse_angle(30, 43, 17)
hera_lon = collapse_angle(21, 25, 42)

def get_lst(lon=None, radians=False):
    """
    Return current local sidereal time (LST)
    for longitude
        @lon : float, degrees BY DEFAULT
        (default is HERA array).
    @radians:
        True: LST is returned in radians, @lon is expected in radians.
        False: LST is returned in degrees, @lon is expected in degrees.
    """
    if lon is None:
        if radians:
            lon = np.radians(hera_lon)
        else:
            lon = hera_lon
    t = astropy.time.Time(time.time(), format='unix')
    if radians:
        return t.sidereal_time('apparent', longitude=lon * u.radian).radian
    return t.sidereal_time('apparent', longitude=lon * u.deg).deg

# The change-of-basis matrix between equatorial and galactic coordinate systems
M_eq_to_gal = np.array([
    [-.054876, -.873437, -.483835],
    [.494109, -.444830, .746982],
    [-.867666, -.198076, .455984]
])

def M_eq_to_ha(lst=None, radians=False):
    """
    Return the change-of-basis matrix between the equatorial and
    hour-angle coordinate systems.
    The conversion depends on the
        local sidereal time @lst : float
    @radians: classifies the units of @lst
    """
    if lst is None:
        lst = get_lst(radians=True)
    if not radians:
        lst = np.radians(lst)
    
    s = np.sin(lst)
    c = np.cos(lst)
    return np.array([[c, s, 0], [s, -c, 0], [0, 0, 1]])

def M_ha_to_topo(phi=None, radians=False):
    """
    Return the change-of-basis matrix between the hour-angle
    and topocentric coordinate systems.
    The conversion depends on the
        latitude of the observer @phi : float
    @radians: classifies the units of @phi
    """
    if phi is None:
        phi = np.radians(hera_lat)
    elif not radians:
        phi = np.radians(phi)
    s = np.sin(phi)
    c = np.cos(phi)
    return np.array([[-s, 0, c], [0, -1, 0], [c, 0, s]])

def rectangle(a, b, radians=False):
    """
    Given a pair of angles
        @a : float
        @b : float
    return the corresponding 3x1 rectangular / Cartesian vector.

    @radians: classifies the units of @a and @b
    """
    if not radians:
        a = np.radians(a)
        b = np.radians(b)
    return np.array([np.cos(b) * np.cos(a), np.cos(b) * np.sin(a), np.sin(b)])

def new_sphere(out_arr, radians=False):
    """
    Given a 3x1 rectangular / Cartesian vector @out_arr,
    return the corresponding pair of angles as a tuple.
    @radians determines whether the angles are given in radians.
    """
    gp = np.arctan2(out_arr[1], out_arr[0])
    tp = np.arcsin(out_arr[2])
    if not radians:
        return np.degrees(gp), np.degrees(tp)   
    return gp, tp

def eq_to_gal(ra, dec, radians=False):
    """
    Convert a position in the equatorial format
        (right ascension = @ra, declination = @dec)
    to the galactic-coordinates position
        (el : galactic longitude, be = galactic latitude)

    @radians determines the interpretation of BOTH the input
    and output. By default everything is in degrees.
    """
    if not radians:
        ra = np.radians(ra)
        dec = np.radians(dec)
    eq_vector = rectangle(ra, dec, radians=True)
    gal_vector = np.dot(M_eq_to_gal, eq_vector)
    return new_sphere(gal_vector, radians)

def eq_to_topo(ra, dec,
    lat=None, lon=None, lst=None, radians=False):
    """
    Convert a position in the equatorial format
        (right ascension = @ra, declination = @dec)
    to the topocentric-coordinates position
        (az : local azimuth, alt = local altitude)

    @radians determines the interpretation of BOTH the input
    and output. By default everything is in degrees.
    """
    if lat is None:
        lat = np.radians(hera_lat)
    elif not radians:
        lat = np.radians(lat)
    if lst is None:
        if lon is None:
            lon = np.radians(hera_lon)
        elif not radians:
            lon = np.radians(lon)
        lst = get_lst(lon, radians=True)
    if not radians:
        ra = np.radians(ra)
        dec = np.radians(dec)
        lst = np.radians(lst)
    eq_vector = rectangle(ra, dec, radians=True)
    ha_vector = np.dot(M_eq_to_ha(lst, radians=True), eq_vector)
    topo_vector = np.dot(M_ha_to_topo(lat, radians=True), ha_vector)
    return new_sphere(topo_vector, radians)

def ha_to_eq(ha, dec, lat, radians=False):
    """
    Convert a position in the hour-angle format
        (hour-angle = @ha, declination = @dec)
    to the equatorial-coordinates position
        (ra : right ascension, dec : declination)

    @radians determines the interpretation of BOTH the input
    and output. By default everything is in degrees.
    """
    if not radians:
        ha = np.radians(ha)
        dec = np.radians(dec)
        lat = np.radians(lat)
    rct = rectangle(ha, dec, radians=True)
    eq = np.dot(np.linalg.inv(M_eq_to_ha(lat, radians=True)), rct)
    return new_sphere(eq, radians)

def ha_to_gal(ha, dec, lst, radians=False):
    """
    Convert a position in the format
        (hour-angle = @ha, declination = @dec),
    and considered at the local sidereal time @lst,
    to the galactic-coordinates position
        (el : galactic longitude, be = galactic latitude)

    Be careful not to confuse the hour-angle coordinate system
        (which this function endeavors to convert)
    with the conventional EQUATORIAL-coordinates representation of
    the right ascension ANGLE in the format (HOUR, minute, second).

    @radians determines the interpretation of BOTH the input
    and output. By default everything is in degrees.
    """
    if not radians:
        ha = np.radians(ha)
        dec = np.radians(dec)
        lst = np.radians(lst)
    rct = rectangle(ha, dec, radians=True)
    ra_dec = np.dot(np.linalg.inv(M_eq_to_ha(lst, radians=True)), rct)
    gal = np.dot(M_eq_to_gal, ra_dec)
    return new_sphere(gal, radians)

def ha_to_topo(ha, dec, lat, radians=False):
    """
    Convert a position in the hour-angle format
        (hour-angle = @ha, declination = @dec),
    and for an observer at latitude @lat,
    to the topocentric-coordinates position
        (az: local azimuth, alt: local altitude)

    Be careful not to confuse the hour-angle coordinate system
    (which this function endeavors to convert)
    with the conventional EQUATORIAL-coordinates representation of
    the right ascension ANGLE in the format (HOUR, minute, second).

    @radians determines the interpretation of BOTH the input
    and output. By default everything is in degrees.

    WARNING: the altitude and azimuth may be listed in the wrong order.
    """
    if not radians:
        ha = np.radians(ha)
        dec = np.radians(dec)
        lat = np.radians(lat)
    rct = rectangle(ha, dec, radians=True)
    topo = np.dot(M_ha_to_topo(lat, radians=True), rct)
    return new_sphere(topo, radians)

def gal_to_eq(el, be, radians=False):
    """
    Convert a position in the galactic format
        (galactic longitude = @el, galactic latitude = @be)
    to the equatorial-coordinates position
        (ra : right ascension, dec = declination)

    @radians determines the interpretation of BOTH the input
    and output. By default everything is in degrees.
    """
    if not radians:
        el = np.radians(el)
        be = np.radians(be)
    rct = rectangle(el, be, radians=True)
    ra_dec = np.dot(np.linalg.inv(M_eq_to_gal), rct)
    return new_sphere(ra_dec, radians)

def gal_to_topo(el, be, lat, lon, radians=False):
    """
    Convert a position in the galactic format
        (galactic longitude = @el, galactic latitude = @be),
    and for an observer at latitude @lat and longitude @lon,
    to the topocentric-coordinates position
        (az : local azimuth, alt : local altitude)

    @radians determines the interpretation of BOTH the input
    and output. By default everything is in degrees.

    WARNING: this function has undergone some dependency changes
    without additional accuracy tests. It may not work correctly, currently.
    """
    if not radians:
        el = np.radians(el)
        be = np.radians(be)
        lat = np.radians(lat)
    rct = rectangle(l, b, radians=True)
    ra_dec = np.dot(np.linalg.inv(M_eq_to_gal), rct)
    lst = get_lst(lon, radians)
    hrd = np.dot(np.linalg.inv(M_eq_to_ha(lst, radians=True)), ra_dec)
    topo = np.dot(M_ha_to_topo(phi, radians=True), hrd)
    return new_sphere(topo, radians)

"""
The following is an adaptation of a function originally
written by C. D. Nunhokee,
'genVisibility.py', polarizedSims, Feb 8 2019
https://github.com/Chuneeta/polarizedSims/blob/master/genVisibility.py
"""
# Unfortunately, the parameter order has been switched wrt the citation.
def radec2lm(ra, dec, ra0=None, dec0=np.radians(hera_lat)):
    """
    Converts equatorial coordinates to direction cosines l & m
    ra   : right ascension in radians; type: float
    dec  : declination in radians; type: float
    ra0  : reference/phase right ascension, in radians; type: float
         (default: current LST)
    dec0 : reference/phase declination; type: float
         (default: latitude of the HERA array)
    """
    # See note at the end about default arguments.
    if ra0 is None:
        ra0 = get_lst(np.radians(hera_lon), radians=True)

    l = np.cos(dec) * np.sin(ra0 - ra)
    m = -1 * (np.sin(dec) * np.cos(dec0) - \
        np.cos(dec) * np.sin(dec0) * np.cos(ra - ra0))
    return l, m

"""
We cannot put ra0 = get_lst() in the function header. Why?
Because Python evaluates all function headers once upon first opening the script.
Consequently, the default argument will be constant: subsequent calls of this
function will use the value of get_lst() calculated when the script was opened.
"""
