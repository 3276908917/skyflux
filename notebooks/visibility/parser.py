import pickle
import numpy as np
import rotations

def test():
    """
    Play around with the return value if you doubt that
    this module is being reloaded properly.
    """
    return 2

J = np.load("../J.npy", fix_imports=False) # Python 2 deserves to die

"""
In principle, it should not matter which flux-by-frequency
you order by. However, in practice (as we can see in
parser_demos.py), some sources are brighter than others
at certain frequencies.

We have several sources that lack spectral indices,
    so it looks like we are moving ahead with plan A
    (use a power-law regression to get a spectral index function)
"""

# The following section is hard-coded to the GLEAMEGCAT format

# all numbers represent MHz quantities
expected_frequencies = [76, 84, 92, 99, 107, 115, 122, 130,
                    143, 151, 158, 166, 174, 181, 189,
                    197, 204, 212, 220, 227]

class GLEAM_entry:
    def __init__(self, line):
        # Might want to redo this line later to exclude universal "GLEAM " prefix
        self.name = line[:line.index("|")]
        line = line[line.index("|") + 1:]
        
        self.ra = line[:line.index("|")]
        self.format_ra()
        line = line[line.index("|") + 1:]

        self.dec = line[:line.index("|")]
        self.format_dec()
        line = line[line.index("|") + 1:]

        self.flux_by_frq = {}

        # we extract and record fluxes according to expected_frequencies
        # at the same time, we convert mJy -> Jy
        for expected_frq in expected_frequencies:
            try:
                self.flux_by_frq[expected_frq] = \
                    float(line[:line.index("|")].strip()) / 1000
            except ValueError:
                print("Missing flux value for:", self.name,
                      "at frequency:", expected_frq, "MHz.")
                self.flux_by_frq[expected_frq] = np.NaN
            line = line[line.index("|") + 1:]

        try:
            self.alpha = float(line[:line.index("|")])
        except ValueError:
            print("Missing spectral index for:", self.name)
            self.alpha = np.NaN

    def format_ra(self):
        remainder = self.ra
        self.ra_hour = float(remainder[:remainder.index(" ")])

        remainder = remainder[remainder.index(" ") + 1:]
        self.ra_minute = float(remainder[:remainder.index(" ")])

        remainder = remainder[remainder.index(" ") + 1:]
        self.ra_second = float(remainder)

        self.ra_angle = rotations.collapse_hour(
            self.ra_hour, self.ra_minute, self.ra_second)

    def format_dec(self):
        remainder = self.dec
        self.dec_degree = float(remainder[:remainder.index(" ")])

        remainder = remainder[remainder.index(" ") + 1:]
        self.dec_arcminute = float(remainder[:remainder.index(" ")])

        remainder = remainder[remainder.index(" ") + 1:]
        self.dec_arcsecond = float(remainder)

        self.dec_angle = rotations.collapse_angle(
            self.dec_degree, self.dec_arcminute, self.dec_arcsecond)

    def __str__(self):
        return "Name: " + self.name + "\nRight ascension: " + str(self.ra_angle) + \
            "\nDeclination: " + str(self.dec_angle) + \
            "\n151 MHz flux: " + str(self.flux_by_frq[151]) + "\n"
    
    # we will probably want a __repr__ function so that we can see
    # ALL fluxes associated with the object.

f = open("gleam_with_alpha.txt", "r")
obj_catalog = []
# For each line in f, the delimiter is |
for line in f:
    obj_catalog.append(GLEAM_entry(line[1:]))
f.close()

# Antenna section

ant_pos = dict(pickle.load(open("ant_dict.pk", "rb")))

def baseline(ant_ID1, ant_ID2):
    """
    Calculate the baseline between antennae # @ant_ID1 and @ant_ID2
    by a simple difference of their coordinates.
    """
    return ant_pos[ant_ID2] - ant_pos[ant_ID1]


S = .5 * np.array([[1, 1, 0, 0,],
                  [0, 0, 1, 1j],
                  [0, 0, 1, -1j],
                  [1, -1, 0, 0]])

def A(source_idx):
    this_J = J[source_idx]
    J_outer = np.kron(J, np.conj(J))
    return np.dot(np.dot(np.linalg.inv(S), J_outer), S)

c = 299792458 # m / s

def phase_factor(ant1, ant2, r, nu=151e6):
    """
    Calculate the phase factor in the direction @r (l, m)
        (we assume that n is of insignificant magnitude)
    and at the frequency @nu
    between two antennae whose ID #s are @ant1 and @ant2.
    When we calculate the baseline (u, v, w), we
        assume that w is of insignificant magnitude.
    """
    b = baseline(ant1, ant2)[0:2] # kill w
    br = np.dot(b, r)
    return np.exp(-2j * np.pi * nu * br / c)

def visibility(ant1, ant2, source_index, nu=151e6):
    """
    Visibility integrand evaluated for a single source.
    """
    source = obj_catalog[source_index]
    
    I = source.flux_by_frq[nu / 1e6]
    s = np.array([complex(I), 0, 0, 0])

    ra = np.radians(source.ra_angle)
    dec = np.radians(source.dec_angle)
    r = rotations.raddec2lm(ra, dec)

    phi = phase_factor(ant1, ant2, r, nu)
    return np.dot(np.dot(A(source_index), s), phi)
        # I do not like this indexing. I only did it to get the scalar,
            # but how can I be sure it is not a sum of elements, for example?

def visibility_integrand(ant1, ant2, nu=151e6):
    total = complex(0) # 4 x 1. Visibility has a phase,
                        # i.e. a phase
    for i in range(len(obj_catalog)):
        total += visibility(ant1, ant2, i, nu)
    return total
    # shouldn't visibility be a real quantity??

    """
        PING RIDHIMA TO GET FREQUENCY BEAM

        for now, assume that 100-200 MHz is the same frequency
            use spectral index to scale the result
            split it at 1 MHz for now.
            (i.e. same Jones matrix 100-200 MHz for now).
        "I will leave the chunking of frequency and time up to you."
        
        We discussed to do ten-minute intervals:
        each ten minutes, we have one data point
        5 hours, 6 data per hour.
        0-8 hours (colder patch)
        create 1D plot.
            great Slack question: what is the best way
            to visually represent the results?

        DO NOT average it yet, return an array of all
        data points.

        time axis (2D), frequency axis (3D)

        write a power law equation which takes a
        spectral index and then scales the result appropriately
            > amplitude of visibility should be decreasing with frequency
            > (since all of these sources are synchrotrons)
        plt.imshow()

        simple test for the computation:
        at zenith (ra = LST), a source should be real
        (phase is zero)

        feel free to use orthview to get the
        Mueller matrix plots. Figure 1 in (Nunhokee et al)
        did not actually use healpy, so feel free to
        tinker.
            "Return projected map -> store in variable
            variable -> imshow"

        read up on healpy.cartview
    """

"""
(I think the output should be a scalar, not a 4x1)

Want graph of times and positions.

A will change according to time. Equation 3 assumes that we have
just one integration time.

The change in A over time represents a drift scan,
like the Gaussians that you recently visualized.
    Time you can do, you have all the information. I am guessing
    I probably only need to call raddec2lm with different values
    of ra0, i.e. moving what would otherwise be the current LST
    over the entire range [0, 2 np.pi)

use constant flux independent of frequency
(this is the same as saying that alpha is equal to one)
then varying frequency from 100-200 MHz
'assume it is the same beam between 100 and 200 MHz'
    beam variation is what we will eventually be handling anyway

15 m or 30 m baseline East-West
    I do not remember what this means.

I need to plot the A matrix, to see the leakage terms.
Plot it in healpix.

(Ask for Ridhima's power spectrum notebook plot for a
single baseline, once the above work has been completed and tested.)
"""

"""
We cannot put ra0 = get_lst() in the function header. Why?
Because Python evaluates all function headers once upon first opening the script.
Consequently, the default argument will be constant: subsequent calls of this
function will use the value of get_lst() calculated when the script was opened.
"""
