import numpy as np
import matplotlib.pyplot as plt
import pickle

import skyflux as sf
import skyflux.deprecated.polSims as pol

def ppp(pp, key, power=1):
    plt.plot(pp['fs'] / 1e6, pp[key] ** power)
    plt.title(key + " vs Frequency, 125-175 MHz band")
    plt.xlabel("Frequency [MHz]")
    
    xu = None
    
    if key == 'lambda':
        xu = "m"
    elif key == 'D':
        xu = "m?"
    elif key == 'DeltaD':
        xu = "m?"
    elif key == 'kperp':
        xu = "h / Mpc"
        
    plt.ylabel(key + " [" + xu + "]^" + str(power))
    
    print("Vertical range is:",
        (pp[key] ** power).max() - (pp[key] ** power).min(),
        xu + "^" + str(power))
        
    print("Vertical ratio is:",
        (pp[key] ** power).max() / (pp[key] ** power).min())
    
    plt.show()

def power_parameters(fname, ant1, ant2):
    fcd, fs, ts, ptitle = load_wedge_sim(fname + ".pickle")

    print("Simulation file loaded.\n")

    num_f = len(fs)
    B = fs.max() - fs.min()
    
    etas = pol.f2etas(fs)
    #z_avg = pol.fq2z(np.average(fs) / 1e9)
    #k_para = pol.k_parallel(etas, z_avg)
    
    k_para = np.array([
        pol.k_parallel(etas[i], pol.fq2z(fs[i] / 1e9)) \
        for i in range(num_f)
    ])

    #print(k_par)

    """ Power constants, section 1"""
    kB = 1.380649e-26 # this is in mK.
    # To use K, add 3 orders of magnitude.
    # 1 Jy = 1e-20 J / km^2 / s^2
    square_Jy = (1e-20) ** 2
    universal_p_coeff = square_Jy / (2 * kB) ** 2 / B
    """ """

    baselength = sf.ant.baselength(ant1, ant2)
        
    power_pieces = {
        'z': [], # dimensionless
        'lambda': [], # m
        'D': [], # m?
        'DeltaD': [], 
        'kperp': []
    }    
        
    for nu in fs:
        #nu = np.average(fs)

        """ Power constants, section 2 """
        z = pol.fq2z(nu / 1e9)
        power_pieces['z'].append(z)
        
        lambda_ = pol.C / nu
        power_pieces['lambda'].append(lambda_)
        
        D = pol.transverse_comoving_distance(z)
        power_pieces['D'].append(D)
        
        DeltaD = pol.comoving_depth(B, z)
        power_pieces['DeltaD'].append(DeltaD)
        
        k_perp = baselength * pol.k_perp(z) / lambda_
        power_pieces['kperp'].append(k_perp)
        
    for key in power_pieces.keys():
        power_pieces[key] = np.array(power_pieces[key])
        
        
    power_pieces['fs'] = fs
        
    return power_pieces
        
def wauto_show(fname, sp=None, pt_override=None, static=False, Qi=None, special_request=None):
    """
    Load simulated visibilities from the file named
        @fname
    and visually interpret the results using a
        2D + color wedge plot.
        
    @sp : Stokes parameter, if you want to look at just one
        at a time
        0 : I    1 : Q    2 : U    3 : V    
        
    @pt_override is a way to override the title with which
        a simulation came. It is extremely bad form to use this,
        but it can come in handy during lapses of diligence.
        
    @special_requests:
        (ant1, ant2) tuple, in case you want to look at a
        k_orth cross-section
    """
    fcd, fs, ts, ptitle = load_wedge_sim(fname + ".pickle")
    
    if Qi is None:
        # hard coding for four Stokes param.s
        Qi = np.identity(4, np.complex128)
    
    if pt_override is not None:
        ptitle = pt_override
    
    print("Simulation file loaded.\n")
    
    transformed = transform_wedge(fcd, fs, ts)
    print("Fourier transforms applied to simulation.\n")
    
    if static:
        wedge = static_visual(sim_dict)
    else:
        wedge = collect_wedge_points(
            transformed, fs, ts, Qi, sp, special_request
        )
        if special_request is not None:
            return wedge
        #return wedge # just for individual baseline testing
    print("Wedge points collected.\n")
    
    return plot_3D(wedge, fs, ptitle)
    
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
    #print(fs)
    
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
                    parameter[ti] = np.fft.fft(
                        parameter[ti] * window
                    )
                
    return fourier_dict

def collect_wedge_points(fcd, fs, ts, Qi, sp=None,
    special_request=None):
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
        
    @sp : Stokes parameter, if you want to look at just one
        at a time
        0 : I    1 : Q    2 : U    3 : V
        
    @special_request: tuple of
        (baseline1, baseline2,    
    to really complete the Nunhokee parallel,
        you would want to average over a range of
        k_orth values 
        
        
    TODO: last updated 3/30
    # we need the horizon line: white dots
    # Fix amplitudes?
    # Create separate plots for I, Q, U, V
        # Create slices after the fashion of Nunhokee fig 6
        
    # make some notes for HERA team circulation
    """
    num_t = len(ts)
    num_f = len(fs)
    B = fs.max() - fs.min()
    
    visual = []
    
    #print(fs)
    
    etas = pol.f2etas(fs)
    #z_avg = pol.fq2z(np.average(fs) / 1e9)
    #k_para = pol.k_parallel(etas, z_avg)
    
    k_para = np.array([
        pol.k_parallel(etas[i], pol.fq2z(fs[i] / 1e9)) \
        for i in range(num_f)
    ])

    #print(k_par)

    """ Power constants, section 1"""
    kB = 1.380649e-26 # this is in mK.
    # To use K, add 3 orders of magnitude.
    # 1 Jy = 1e-20 J / km^2 / s^2
    square_Jy = (1e-20) ** 2
    universal_p_coeff = square_Jy / (2 * kB) ** 2 / B
    """ """

    for ant1 in fcd.keys():
        for ant2 in fcd[ant1].keys():
            if special_request is not None:
                if ant1 != special_request[0] or \
                    ant2 != special_request[1]:
                    continue
        
            baselength = sf.ant.baselength(ant1, ant2)
            
            special = [[], [], [], []]
                
            for nu_idx in range(num_f):
                nua = np.average(fs)
                nu = fs[nu_idx]

                """ Power constants, section 2 """
                z = pol.fq2z(nu / 1e9)
                lambda_ = pol.C / nu
                lambda_a = pol.C / nua
                D = pol.transverse_comoving_distance(z)
                DeltaD = pol.comoving_depth(B, z)
                # Finally, condense everything into a
                # power coefficient
                p_coeff = universal_p_coeff * \
                    lambda_ ** 4 * D ** 2 * DeltaD
                """ """
                k_perp = baselength * pol.k_perp(z) / lambda_
                
                powers_prop = []
                # store power results by Stokes index:
                special_powers = [[], [], [], []]
                # store normalized Hadamard vectors:
                special_times = []
                
                for t_idx in range(num_t - 1):
                    this_instant = \
                        fcd[ant1][ant2][:, t_idx, nu_idx]
                    next_instant = \
                        fcd[ant1][ant2][:, t_idx + 1, nu_idx]
                    
                    # [I1, Q1, U1, V1] * [I2*, Q2*, U2*, V2*]
                    
                    hadamard = np.multiply(
                        this_instant, next_instant
                    )
                    p = np.dot(Qi, hadamard)
                    
                    # if no Stokes parameter is specified,
                    # we take the absolute value of the
                    # entire Stokes visibility vector
                    if sp is None:
                        sqBr = np.abs(p)
                    else:
                        sqBr = np.abs(p[sp])
                    
                    powers_prop.append(sqBr)
                    
                    # we leave the normalized Hadamard
                    # products for the cross-section
                    # code to handle
                    special_times.append(p)
                    
                if special_request is not None:
                    for vector in np.array(special_times):
                        # si: Stokes index
                        for si in range(len(vector)):
                            # take the magnitude of each
                            # normalized Hadamard index
                            param = np.abs(vector[si])
                            # and count it as a power
                            special_powers[si].append(param)

                avg = p_coeff * np.average(np.array(powers_prop))
                
                #print("Using k_parallel", k_par[nu_idx])
                
                wedge_datum = np.array([
                    k_perp,
                    k_para[nu_idx],
                    #float(avg)
                    float(np.log10(avg))
                ])
                
                if special_request is not None:
                    for si in range(len(special_powers)):
                        special_powers[si] = np.array(
                            special_powers[si])
                    
                    special_powers = np.array(special_powers)
                
                    #print("Using k_parallel", k_par[nu_idx])
                
                    # si: Stokes index
                    for si in range(len(special_powers)):
                        stokes_param = special_powers[si]
                        avg = p_coeff * np.average(stokes_param)
                        
                        #!!! duplicate reference
                        special[si].append(np.array([
                            k_para[nu_idx],
                            #float(avg)
                            float(np.log10(avg))
                        ]))
                    
                visual.append(wedge_datum)
            
            # Figure 6-esque investigation    
            if special_request is not None:
                if ant1 == special_request[0] and \
                    ant2 == special_request[1]:
                    print("Exiting")
                    return np.array(special)

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
   
def plot_3D(visual, fs, title, scaled=False):
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

    print("Minimum:", z.min())
    print("PTP:", z.ptp())

    image = visual_to_image(visual)

    ### horizon-code
    center_f = np.average(fs)
    z = pol.fq2z(center_f / 1e9)
    lambda_ = pol.C / center_f
    k_starter = pol.k_perp(z) / lambda_
    
    horizonx = np.unique(visual[:, 0])
    horizonyp = []
    horizonym = []
    
    etas = pol.f2etas(fs)    
    k_par = pol.k_parallel(etas, z)
    
    for i in range(len(horizonx)):
        k_orthogonal = horizonx[i]
        
        horizonyp.append(pol.horizon_limit(k_orthogonal, z))
        horizonym.append(-pol.horizon_limit(k_orthogonal, z))
        
    plt.scatter(horizonx, horizonyp, marker='.', c='w')
    plt.scatter(horizonx, horizonym, marker='.', c='w')
    ### end horizon-code
    
    # Interpolation choice doesn't really matter;
    # gaussian looks smooth
    plt.imshow(image.T, extent=[
        x.min(), x.max(), y.min(), y.max()
    ], interpolation='gaussian', aspect='auto')
    
    finalize_plot(title)

    plt.scatter(x, y, marker='.', c=colors)
    finalize_plot(title)
    
    return x, y, visual, image
    
def finalize_plot(title):
    """
    Add automatic labels to the plot.
    
    This function chiefly serves as convenience for the
    programmer, due to its almost entire reliance on
    hard-coded strings.
    """
    cbar = plt.colorbar()
    plt.title(title)
    plt.xlabel("$k_\perp$ [$h$ Mpc$^{-1}$]")
    plt.ylabel("$k_\parallel$ [$h$ Mpc$^{-1}$]")
    cbar.set_label("log$_{10}$ [mK$^2$ ($h^{-1}$ Mpc)^3] ?")
    
    plt.show()

# hard-coding for now
# automatically assumes a full sky
    # ie phi \in [0, 2 * np.pi]
    # and theta \in [0, np.pi]
def calculate_Q(
    B=np.arange(50e6, 250e6 + 0.001, 4e6),
    angular_resolution = 250
):
    
    print("Establishing integration parameters")
    dnu = B[1] - B[0]
    window = pol.genWindow(len(B))
    
    list_phi = np.linspace(0, 2 * np.pi, 4 * angular_resolution)
    dphi = list_phi[1] - list_phi[0]
    
    list_theta = np.linspace(0, np.pi / 2, angular_resolution)
    dtheta = list_theta[1] - list_theta[0]
    
    dummy_A = sf.stokes.create_A(
        az = 0, alt = 0, nu=151e6, radians=True
    )
    Q = np.zeros(dummy_A.shape, dtype=np.complex128)

    d3 = dnu * dphi * dtheta
    
    print("Integration paramaters established. Integrating...")

    for phi in list_phi:
        print("phi is", phi)
        for theta in list_theta:
            for nu_idx in range(len(B)):
                next_A = sf.stokes.create_A(
                    az=phi, alt=theta, nu=B[nu_idx], radians=True
                )
                Aw = next_A * window[nu_idx]
                Q += np.multiply(Aw, np.conj(Aw)) * d3
    
    print("Integration complete.")
    
    #! Disgusting hack
    return Q[:, 0, :]

def diag(matrix):
    # Be careful! We are assuming that a shallow copy is fine
    diagd = np.copy(matrix)
    
    for i in range(len(diagd)):
        for j in range(len(diagd[i])):
            if i != j:
                diagd[i][j] = 0
    
    return diagd
                
