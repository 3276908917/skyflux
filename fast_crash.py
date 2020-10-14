#import skyflux as sf
import matplotlib.pyplot as plt
import numpy as np
import healpy as hp

def stokes_plotter(dv_dict, dv_label):
    plt.clf()
    
    x = dv_dict[dv_dict['label']]
    dv_list = dv_dict[dv_label]

    im = np.imag(dv_list)
    re = np.real(dv_list)
    ab = np.abs(dv_list)
    plt.plot(x, im, label="imaginary")
    plt.plot(x, re, label="real")
    plt.plot(x, ab, label="abs")

    plt.title(dv_label + " over " + dv_dict['label'])
    plt.legend()
    plt.xlabel(dv_dict['label'] + " [rad]")
    plt.ylabel(dv_label + " [Jy]")

    plt.show()

def shell_crasher(create_J):
    
    dec_decf = np.radians(-30.72)
    ra_decf = 0
    list_xx_over_ra = []
    list_yy_over_ra = []
    list_ra = []

    while ra_decf <= 120:
        ra = np.radians(ra_decf)

        # do a manual shell %run to mimic stokes import
        J_decf = create_J(ra=ra, dec=dec_decf, lst=np.pi/3, radians=True)
        list_ra.append(ra)
        list_xx_over_ra.append(J_decf[0][0][0])
        list_yy_over_ra.append(J_decf[0][1][1])
        ra_decf += 10

    dict_ra = {
        'label': 'RA',
        'RA': list_ra,
        'xx': np.array(list_xx_over_ra),
        'yy': np.array(list_yy_over_ra)
    }    
