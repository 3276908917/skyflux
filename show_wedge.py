import numpy as np
import matplotlib.pyplot as plt
import pickle

import skyflux as sf
import skyflux.deprecated.polSims as pol

def wauto_show(fname, static=False):
    """
    Load simulated visibilities from the file named
        @fname
    and visually interpret the results using a
        2D + color wedge plot.
    """
    fcd, fs, ts, ptitle = load_wedge_sim(fname + ".pickle")
    print("Simulation file loaded.\n")
    
    transformed = transform_wedge(fcd, fs, ts)
    print("Fourier transforms applied to simulation.\n")
    
    if static:
        wedge = static_visual(sim_dict)
    else:
        wedge = collect_wedge_points(transformed, fs, ts)
    print("Wedge points collected.\n")
    
    return plot_3D(wedge, ptitle)
    
def load_wedge_sim(fname):
    """
    Fill out a wedge data structure based on
        the contents of the file with name @fname.
    The returned dictionary has the following key-value pairs:
        "fs": the frequencies used in the simulation
        "kp": array of k-parallel, the proxy-variable for frequency
        "ks": a constant coefficient, to be multiplied by
            baseline lengths to yield k-perpendicular,
            the proxy-variable for baseline
        "sim": a three-dimensional array of visibilities
            the first index corresponds to frequency
            the second index corresponds to LST
            the third index corresponds to Stokes parameter
                0: I
                1: Q
                2: U
                3: V
    """
    sim_file = open(fname, "rb")

    meta = pickle.load(sim_file)

    ptitle = meta['title']

    fs = meta['frequencies']
    num_f = len(fs)
    
    ts = meta['times']
    num_t = len(ts)
    
    sim = meta['picture']
    
    fcd = {} # fourier candidate dictionary
    
    for ant1 in sim.keys():
        fcd[ant1] = {}
        for ant2 in sim[ant1].keys():
            # 0: I    1: Q    2: U    3: V
            fourierc = [[], [], [], []]
            
            for ti in range(num_t):
                for parameter in fourierc:
                    parameter.append([])
                
                for ni in range(num_f):
                    v = sim[ant1][ant2][ni][ti]

                    for p_idx in range(len(fourierc)):
                        fourierc[p_idx][ti].append(v[p_idx])
                
                for parameter in fourierc:
                    parameter[ti] = np.array(parameter[ti])
        
            for parameter in fourierc:
                parameter = np.array(parameter)

            fcd[ant1][ant2] = np.array(fourierc)

    return fcd, fs, ts, ptitle    
    
def static_wedge_vis(sim_dict, fs):
    raise NotImplementedError("I have only updated" + \
    " the dynamic wedge routine at this moment.")
    
    """
    Read from the wedge data structure @sim_dict
        (specifically, one using the format
        established in load_wedge_sim)
    and generate 3D points appropriate for a wedge plot.
    i.e., return a list of triples:
        (k_perpendicular, k_parallel, power*)
    
    * currently, the implementation uses values that
    should be proportional to power. The final constants
    have not yet been considered.
    
    This function is distinct from
        dynamic_wedge_vis
    in assuming that the simulation corresponding to
        @sim_dict
    runs over only one value of LST.
    """
    wedge_data = []
    
    etas = f2etas(fs)    

    center_f = np.average(fs)
    z = fq2z(center_f)
    lambda_ = C / center_f

    k_par = k_parallel(etas, z)
    k_starter = k_perp(z) / lambda_ # this will need to be
    # multiplied on a per-baseline basis
    
    ### aliasing ###
    nu_idxs = sim_dict['fs']

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
   
def transform_wedge(original, fs, ts):
    num_f = len(fs)
    num_t = len(ts)
    
    import copy
    
    fourier_dict = copy.deepcopy(original)

    window = pol.genWindow(num_f)

    for ant1 in fourier_dict.keys():
        for ant2 in fourier_dict[ant1].keys():
            fourier = fourier_dict[ant1][ant2]
            for ti in range(num_t):
                for parameter in fourier:
                    parameter[ti] = np.fft.fft(parameter[ti] * window)
                
    return fourier_dict
    
#! There has to be some way to merge this with the function
# above
def collect_wedge_points(fcd, fs, ts):
    """
    Read from the wedge data structure @sim_dict
        (specifically, one using the format
        established in load_wedge_sim)
    and generate 3D points appropriate for a wedge plot.
    i.e., return a list of triples:
        (k_perpendicular, k_parallel, power*)
    
    * currently, the implementation uses values that
    should be proportional to power. The final constants
    have not yet been considered.
    
    This function is distinct from
        static_wedge_vis
    in accepting multiple LST values from
        @sim_dict.
    """
    num_t = len(ts)
    num_f = len(fs)
    
    visual = []
    
    etas = pol.f2etas(fs)    

    center_f = np.average(fs)
    
    z = pol.fq2z(center_f / 1e9)
    lambda_ = pol.C / center_f

    k_parSym = pol.k_parallel(etas, z)
    #!! Does this address the underlying problem?
    #k_par = np.add(k_parSym, -k_parSym.min())
    # the negative and positive value thing is fine;
    # the wedge plot in Nunhokee is just the positive half
    # of what you already have
    
    k_starter = pol.k_perp(z) / lambda_ # this will need to be
    # multiplied on a per-baseline basis
    
    # Power constants
    B = 200e6 # 200 MHz
    D = sf.deprecated.polSims.transverse_comoving_distance(z)
    DeltaD = sf.deprecated.polSims.comoving_depth(B, z)
    kB = 1.380649e-23
    
    # 1 Jy = 1e-26
    square_Jy = (1e-26) ** 2
    
    # This is an approximation for the normalization volume
    Thyagarajan = 1 / 2 / np.pi / B
    
    p_coeff = (lambda_ ** 2 / 2 / kB) ** 2 * \
         D ** 2 * DeltaD / B * square_Jy #* Thyagarajan
    # end

    for ant1 in fcd.keys():
        for ant2 in fcd[ant1].keys():
            k_orth = k_starter * sf.ant.baselength(ant1, ant2)
            
            special = []
                
            for nu_idx in range(num_f):
                powers_prop = []
                
                for t_idx in range(num_t - 1):
                    this_instant = \
                        fcd[ant1][ant2][:, t_idx, nu_idx]
                    next_instant = \
                        fcd[ant1][ant2][:, t_idx + 1, nu_idx]
                    
                    # [I1, Q1, U1, V1] * [I2*, Q2*, U2*, V2*]
                    
                    # todo:
                    # Fix amplitudes and axes
                    # Create separate plots for I, Q, U, V
                    # Create slices after the fashion of
                        # Nunhokee figure 6
                    # do some linear interpolation via imshow
                        # to make the plot more continuous
                        
                    # try UV tools: it's a wrapper for
                        # imshow
                        
                    # make some notes for HERA team
                    # circulation
                    
                    powers_prop.append(np.abs(np.vdot(
                        this_instant,
                        next_instant
                    )))

                avg = p_coeff * np.average(np.array(powers_prop))
                
                wedge_datum = np.array([
                    k_orth,
                    k_par[nu_idx],
                    float(np.log10(avg))
                ])
                
                special.append(np.array([
                    k_par[nu_idx],
                    float(np.log10(avg))
                ]))
                
                """
                dspecvec = np.array([
                parameter[ti][ni] for parameter in fouriered
                ])
            
                norm = np.linalg.norm(dspecvec)

                visual.append(np.array((
                    etas[ni] * 1e9,
                    ts[ti] * 12 / np.pi,
                    np.log10(norm)
                )))
                """

                visual.append(wedge_datum)
            
            # Figure 6-esque investigation    
            #if ant1 == 1 and ant2 == 50:
            #    return np.array(special)

    visual = np.array(visual)
   
    return np.array(visual)
 
def visual_to_image(visual):
    # We need to order these lists for the wedge to cohere
    xu = np.unique(visual[:, 0])
    yu = np.unique(visual[:, 1])
    
    image = np.zeros((len(xu), len(yu)))
    
    for v in visual:
        xi = np.where(xu == v[0])[0][0]
        yi = np.where(yu == v[1])[0][0]
        image[xi][yi] += v[2]
    
    return image 

"""
def visual_to_image(visual):
    import scipy.interpolate
    x = visual[:, 0]
    y = visual[:, 1]
    z = visual[:, 2]

    xi = np.linspace(x.min(), x.max(), 1000)
    yi = np.linspace(y.min(), y.max(), 1000)
    
    xi, yi = np.meshgrid(xi, yi)
    
    zi = scipy.interpolate.griddata((x, y), z,
        (xi, yi), method='cubic')
    #rbf = scipy.interpolate.Rbf(x, y, z, function='linear')
    #zi = rbf(xi, yi)
    
    plt.imshow(zi, vmin=z.min(), vmax=z.max(), origin='lower',
        )
    #plt.scatter(x, y, c=z)
    plt.colorbar()
    plt.show()
"""
   
def plot_3D(visual, title, scaled=False):
    """
    Primitive 3D plotter.
    
    For use with the return value of either
        static_wedge_vis
    or
        dynamic_wedge_vis
        
    Disable @scaled if you are using values such as logarithms
    """
    x = visual[:, 0]
    y = visual[:, 1]
    z = visual[:, 2]

    colors = None
    if (scaled):
        scaled_z = (z - z.min()) / z.ptp()
        colors = plt.cm.viridis(scaled_z)
    else:
        colors = z

    plt.title(title)

    print("Minimum:", z.min())
    print("PTP:", z.ptp())

    image = visual_to_image(visual)

    plt.imshow(image.T, extent=[
        x.min(), x.max(), y.min(), y.max()
    ], interpolation='nearest', aspect='auto')

    plt.colorbar()
    plt.show()

    plt.scatter(x, y, marker='.', c=colors)
    
    plt.colorbar()
    plt.show()
    
    """
    visual_to_image(visual)
    
    
    import matplotlib.tri as tri
    plt.tricontour(x, y, colors, 15, linewidths=0.5, colors='k')
    plt.tricontourf(x, y, z, 15)
    
    plt.show()
    """
    
    return x, y, visual, image
    
