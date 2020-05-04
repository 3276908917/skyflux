import numpy as np

def collapse_angle(degree, minute, second):
    return degree + minute / 60 + second / 3600

def collapse_hour(hour, minute, second):
    return 15 * hour + minute / 4 + second / 240

lat = -collapse_angle(30, 43, 17)
lon = collapse_angle(21, 25, 42)

# The change-of-basis matrix between equatorial and galactic coordinate systems
M_eq_to_gal = np.array([
    [-.054876, -.873437, -.483835],
    [.494109, -.444830, .746982],
    [-.867666, -.198076, .455984]
])

def gal_to_eq(el, be, lat=ugradio.nch.lat, radians=False):
    if not radians:
        l = np.radians(el)
        b = np.radians(be)
        phi = np.radians(lat)
    else:
        l = el
        b = be
        phi = lat
    rct = rectangle(l, b)
    ra_dec = np.dot(np.linalg.inv(M_eq_to_gal), rct)
    return new_sphere(ra_dec, radians)

def eq_to_gal(ra, dec, latitude, lst, radians=False):
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

def gal_to_topo(el, be, jd, lat=ugradio.nch.lat,
    lon=ugradio.timing.nch.lon, radians=False):
    '''
    @radians determines the format of BOTH input and output!
    Given a pair of angles @el and @be (in galactic coordinates),
    return a pair of angles relating the associated
    azimuth and altitude.
    '''
    if not radians:
        l = np.radians(el)
        b = np.radians(be)
        phi = np.radians(lat)
        # The lst function expects radians,
        # so we do not convert this quantity.
        theta = lon
    else:
        l = el
        b = be
        phi = lat
        theta = np.degrees(lon)
    rct = rectangle(l, b)
    ra_dec = np.dot(np.linalg.inv(M_eq_to_gal), rct)
    lst = ugradio.timing.lst(jd, theta)
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

def ha_to_topo(ha, dec, lat=ugradio.nch.lat, radians=False):
    '''
    Take a position in hour-angle right ascension / declination
        to local altitude and azimuth.
    This performs NO precession.
    '''
    if not radians:
        r = np.radians(ha)
        d = np.radians(dec)
        phi = np.radians(lat)
    else:
        r = ha
        d = dec
        phi = lat
    rct = rectangle(r, d)
    topo = np.dot(M_ha_to_topo(phi), rct)
    return new_sphere(topo, radians)

def ha_to_eq(ha, dec, lat=ugradio.nch.lat, radians=False):
    '''
    Take a position in hour-angle right-ascension / declination
        to regular right-ascension / declination.
    '''
    if not radians:
        r = np.radians(ha)
        d = np.radians(dec)
        phi = np.radians(lat)
    else:
        r = ha
        d = dec
        phi = lat
    rct = rectangle(r, d)
    eq = np.dot(np.linalg.inv(M_eq_to_ha(phi)), rct)
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
