from scipy import integrate
from scipy.optimize import minimize, least_squares
import numpy as np
import matplotlib.pyplot as plt
import sncosmo
import os
import pandas as pd
from astropy.coordinates import Distance
import astropy.units as u
import warnings
from bisect import bisect_left, bisect_right
import matplotlib.ticker as tick
import matplotlib.font_manager as font_manager

h = 6.62607015*1e-34 # кг·м2·с−1
c = 299792458   #  м / с
k = 1.380649*1e-23 # Дж/К
Jy = 1e-26 # W /(Hz * m**2) 
stef_bol = 5.670367*1e-8   #Вт * м-2 * К4


# Colours for plots
cols = {'u': 'dodgerblue', 'g': 'g', 'r': 'r', 'i': 'black', 'z': 'k', 'y': '0.5',
        'Y': '0.5', 'U': 'slateblue', 'B': 'b', 'V': 'yellowgreen', 'R': 'crimson', 'G': 'salmon',
        'I': 'chocolate', 'J': 'darkred', 'H': 'orangered', 'K': 'saddlebrown',
        'S': 'mediumorchid', 'D': 'purple', 'A': 'midnightblue',
        'F': 'hotpink', 'N': 'magenta', 'o': 'darkorange', 'c': 'cyan', 
        'UVW1': 'darkorchid', 'UVW2': 'orchid', 'UVM2': 'deeppink',
        "u'": 'dodgerblue', "g'": 'g', "r'": 'r', "i'": 'black', "z'": 'k'}

#Effective wavelengths (in Angs)
wle = {'u': 3560,  'g': 4830, 'r': 6260, 'i': 7670, 'z': 8890, 'y': 9600, 'Y': 9600,
       'U': 3600,  'B': 4380, 'V': 5450, 'R': 6410, 'G': 6730, 'I': 7980, 'J': 12200, 'H': 16300,
       'K': 21900, 'S': 2030, 'D': 2231, 'A': 2634, 'F': 1516, 'N': 2267, 'o': 6790, 'c': 5330}

#Filter widths (in Angs)
width = {'u': 458,  'g': 928, 'r': 812, 'i': 894,  'z': 1183, 'y': 628, 'Y': 628,
         'U': 485,  'B': 831, 'V': 827, 'R': 1389, 'G': 4203, 'I': 899, 'J': 1759, 'H': 2041,
         'K': 2800, 'S': 671, 'D': 446, 'A': 821,  'F': 268,  'N': 732, 'o': 2580, 'c': 2280}

def width_nu(band):
    l = wle[band] # meters
    delta_l = width[band]
    delta_nu = c * delta_l / l**2 * 1e+10
    return delta_nu

######################################################
# bands - ugriz

bands_sdss = ['u','g','r','i','z']
band_sdss = []
for b in bands_sdss:
    band_sdss.append( sncosmo.get_bandpass('sdss'+b) )



#bands - UBVRI

#bands_snls = ['u','b','v','r','i']
# band_snls = []
# for b in bands_snls:
#     band_snls.append( sncosmo.get_bandpass('standard::'+b) )
bands_snls = ['U','B','V','R','I']
bessel_path = os.path.abspath("/home/timofey/saimsu/SLSN_2021/bessel_bands.csv")
ubvri = pd.read_csv(bessel_path)
band_snls = []
for b in bands_snls:
    band_snls.append( ubvri[b] )
    
L = np.array(ubvri['l']) #in Angs
L *= 1e-10 #in m
#L = np.linspace( 3000, 11000, 5000 ) * 1e-10  # м
L_min, L_max = L[0], L[-1]
W = c / L
W_min, W_max = W[-1], W[0]


#bands - swift-uvot
bands_uvot = ['uvm2','uvw1','uvw2']
band_uvot = []
for b in bands_uvot:
    band_uvot.append( sncosmo.get_bandpass('uvot::'+b) )


def hat(l):
    if l >= 2000 and l <= 10000:
        return 1
    else: 
        return 0

UBVRI, uvot = ['U','B','V','R','I'], ['UVM2', 'UVW1', 'UVW2']
ugriz_prime = ["u'", "g'", "r'", "i'", "z'"]
bands = dict( zip( bands_sdss + ugriz_prime + UBVRI + ["hat"] #+ uvot
                  , band_sdss + band_sdss + band_snls + [hat]  ) )#+ band_uvot
# Wl = np.linspace( 3000, 11000, 100000 )
# maximum_bands = [ max(bands[b](Wl)) for b in  bands_sdss + UBVRI]
# maximum_bands = dict( zip( bands_sdss + UBVRI, maximum_bands ) )
######################################################
#VEGA
# vega_path = os.path.abspath("/media/documents/гаиш/SLSN_2021/alpha_lyr_stis_005.ascii")
vega_path = os.path.abspath("/home/timofey/saimsu/SLSN_2021/alpha_lyr_stis_005.ascii")

vega = []
with open(vega_path) as f:
    for line in f:
        vega.append([float(x) for x in line.split()])
vega = np.array(vega)
vega = vega[(vega[:,0]  >= 3000)*(vega[:,0]  <= 11000)] #cut by range wavelength

def Vega(w):
    l = c / w
    l_ang = l * 1e+10 #angstrom
    #nearest = min(vega[:,0], key=lambda x: abs(x-l_ang))
    #ans = vega[(vega[:,0]) == nearest][0,1]
    ind = bisect_right(vega[:,0], l_ang) - 1
    ans = vega[ind, 1] * 1e+7 #f_lambda в си
    ans *= l**2 / c
    return ans
######################################################


def plank_l(l, T):
    f = 2*h*c**2 / l**5 / ( np.exp(h*c/(l*k*T)) - 1 )
    return f

def ABmagn_to_flux(mag):
    return np.power(10, -(48.6 + mag)/2.5)

def flux_to_magn(flux):
    magn = - 2.5 * np.log10(flux)
    #err = 2.5 / np.log(10) / flux * err
    return magn#, err

def magn_to_flux(magn):
    return np.power(10, -0.4 * (magn))


def plank(w, T):
    f = 2*h*w**3 / c**2 / ( np.exp(h*w/(k*T)) - 1 )
    return f

def flux(T):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        f = integrate.quad(plank, W_min, W_max, args=(T,))[0]
    return f

def Band(w, b):
    #возвращает пропускание фильтра b на частоте w
    l = c / w
    l_ang = l * 1e+10 #angstrom
    
    if b[0].islower():
        return bands[b](l_ang) #/ maximum_bands[b]
    else:
        ind = bisect_right(L*1e10, l_ang) - 1
        ans = bands[b][ind]
        return ans


def band_Plank(w, b, T, z):
    w0 = w
    w *= (1+z)
    f = Band(w0,b) * 2*h*w**3 / c**2 / (np.exp( h*w/(k*T) ) - 1) / (1+z)
    
    return f

def band_Plank_prime_T(w, b, T, z):
    w0 = w
    w *= (1+z)
    f = Band(w0,b) * 2*h**2*(w**4)/(c**2*k)*np.exp( h*w/(k*T) ) \
        / (np.exp( h*w/(k*T) ) - 1)**2 / T**2 / (1+z)
    
    return f

# second deriv
def band_Plank_prime_T_2(w, b, T, z):
    w0 = w
    w *= (1+z)
    f = Band(w0,b) * 2*h**2*w**4 / (k*c**2) * np.exp( h*w/(k*T) ) \
        * ( (np.exp( h*w/(k*T) ) - 1)**2 * (-2*T - h*w/k) + 
           2*(np.exp( h*w/(k*T) ) - 1)*h*w/k ) / T**4   \
            / (np.exp( h*w/(k*T) ) - 1)**4 / (1+z)
    return f

def norm_band(w, b):
    if b[0].islower():
        f = 3631*Jy*Band(w,b)
    else:
        f = Vega(w)*Band(w,b)
    return f

const = [1e4, 1e13]
def band_mag(b, x, z):
    #T, R = x
    T, R = x[0]*const[0], x[1]*const[1]
    d = Distance(z=z, unit=u.m).to_value()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        f = -2.5 * np.log10( integrate.quad(band_Plank, W_min, W_max, args=(b,T,z,), 
                                      limit=50)[0] \
        / integrate.quad(norm_band, W_min, W_max, args=(b,), 
                                      limit=50)[0] * R**2 / d**2)

    return f

def band_mag_primes(x, b, z):
    #T, R = x
    T, R = x[0]*const[0], x[1]*const[1]
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        f = -2.5 * integrate.quad(band_Plank_prime_T, W_min, W_max, 
                                                     args=(b,T,z,))[0] \
            / integrate.quad(band_Plank, W_min, W_max, args=(b,T,z))[0]
    return (f/np.log(10)*const[0], -5/R/np.log(10)*const[1])

#second mag primes
def band_mag_primes_2(x, b, z):
    #T, R = x
    T, R = x[0]*const[0], x[1]*const[1]
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        f = -2.5 / np.log(10) * (integrate.quad(band_Plank_prime_T_2, W_min, W_max, 
                                                     args=(b,T,z,))[0] \
            / integrate.quad(band_Plank, W_min, W_max, args=(b,T,z))[0] \
                - (integrate.quad(band_Plank_prime_T, W_min, W_max, 
                                                     args=(b,T,z,))[0] \
            / integrate.quad(band_Plank, W_min, W_max, args=(b,T,z))[0])**2)
                
    return (f*const[0]**2 , 5/R**2/np.log(10)*const[1]**2)


def band_mag_jac(x, data, z):
    T, R = x
    list_b = data.index[1:]
    list_b = list(set(list_b) & set(list(bands.keys())))
    f = ( np.sum( [ -2 * (data[b] - band_mag(b, x, z)) *
                   band_mag_primes(x,b,z)[0]  
               for b in list_b ] ), 
         np.sum( [ -2 * (data[b] - band_mag(b, x, z)) *   #минус перед скобкой из за минуса внутри скобок
                   band_mag_primes(x,b,z)[1]  
               for b in list_b ] ) )
    
    return f

def band_mag_hess(x, data, z):
    T, R = x
    list_b = data.index[1:]
    list_b = list(set(list_b) & set(list(bands.keys())))
    f = ( (np.sum( [ -2 * (data[b] - band_mag(b, x, z)) *
                   band_mag_primes_2(x,b,z)[0] +
                   2 * (band_mag_primes(x,b,z)[0])**2
               for b in list_b ] ),
           
            np.sum( [ 2 * band_mag_primes(x,b,z)[0] * band_mag_primes(x,b,z)[1]
                     for b in list_b] ) ),
         
         ( np.sum( [ 2 * band_mag_primes(x,b,z)[0] * band_mag_primes(x,b,z)[1]
                     for b in list_b] ),
          
          np.sum( [ -2 * (data[b] - band_mag(b, x, z)) *
                   band_mag_primes_2(x,b,z)[1] +
                   2 * (band_mag_primes(x,b,z)[1])**2
               for b in list_b ] ) ) )
    
    return f


def sum_squared_err(x, data, z):
    T, R = x
    list_b = data.index[1:]
    f = 0
    for b in list_b:
        if b in bands:      #есть ли фильтр в нашей базе
            f += (data[b] - band_mag(b, x, z))**2   
    return f



###############################################################################
#plots
def plot_luminosity(slsn, z, data, save=0):
    data_size = len(data["T"])
    L = []
    for i in range(data_size):
        L.append( flux(data["T"][i])*4*np.pi*data["R"][i]**2 * 1e+7 )
    L = np.array(L)

    max_day = data["mjd"][np.argmax(L)]
    mask = (data["mjd"] <= max_day + 150)
    mjd = np.array(data["mjd"][mask]) / (z + 1)
    
    fig, ax = plt.subplots() #figsize=(18, 12),dpi=400
    plt.plot(mjd, np.log10(L[mask]))
    plt.xlabel('mjd')
    plt.ylabel('$lg_{10}L  [erg/s]$')
    plt.title(slsn)
    if save == 1:
        #fig.savefig( fname='/media/documents/гаиш/SLSN_2021/bol_figures/' + 
                    #slsn, bbox_inches="tight")
        fig.savefig( fname='/home/timofey/saimsu/SLSN_2021/fix_x_bol_output/figures/blc/' + 
                    slsn, bbox_inches="tight", format='pdf')

def y_fmt(x, y):
    return '{:2.1e}'.format(x) #.replace('e', 'x10^')

def plot_sub(slsn, z, data, save=0):
    data_size = len(data["T"])
    L = []
    for i in range(data_size):
        L.append( flux(data["T"][i])*4*np.pi*data["R"][i]**2 * 1e+7 )
    L = np.array(L)

    max_day = data["mjd"][np.argmax(L)]
    mask = (data["mjd"] <= max_day + 150)
    mjd = np.array(data["mjd"][mask]) / (z + 1)
    mjd = mjd - mjd[np.argmax(L)]
    
    fig, axs = plt.subplots(3,sharex=True,figsize=(10, 7),
                            gridspec_kw={'height_ratios': [2, 1.2, 1.2]}) #figsize=(18, 12),dpi=400
    #gs = fig.add_gridspec(3, hspace=0)
    #axs = gs.subplots(sharex=True, sharey=True)
    
    formatter =tick.FormatStrFormatter ("%.1f")
    axs[0].yaxis.set_major_formatter (formatter)
    axs[1].yaxis.set_major_formatter(tick.FuncFormatter(y_fmt))
    axs[2].yaxis.set_major_formatter(tick.FuncFormatter(y_fmt))
    
    axs[0].plot(mjd, np.log10(L[mask]), c='royalblue')
    #plt.xlabel('mjd')
    axs[0].set_ylabel('$lg_{10}L$  $[erg/s]$')
    #plt.title(slsn)
    axs[1].plot(mjd, data["R"][mask]*100, c='royalblue')
    axs[1].set_ylabel('$R$  $[cm]$')
    
    axs[2].plot(mjd, data["T"][mask], c='royalblue')
    axs[2].set_ylabel('$T$  $[K]$')
    axs[2].set_xlabel('rest frame')
    plt.subplots_adjust(hspace=0)
    axs[0].set_title(slsn)
    
    axs[1].xaxis.set_ticks_position('top')
    axs[2].xaxis.set_ticks_position('both')
    for ax in axs:
        ax.label_outer()
    
    for ax in axs:
        ax.grid('on', linestyle='--', alpha=0.5)
    if save == 1:
        #fig.savefig( fname='/media/documents/гаиш/SLSN_2021/bol_figures/' + 
                    #slsn, bbox_inches="tight")
        fig.savefig( fname='/home/timofey/saimsu/SLSN_2021/bol_output_fixedz/figures/' + 
                    slsn, bbox_inches="tight", format='pdf')
        plt.close()


def presets_fig(ax):
    font = {'family' : 'Times New Roman',
    #'weight' : 'bold',
    'size'   : 22}
    plt.rcParams['axes.linewidth'] = 1.2
    plt.rc('font', **font)
    plt.rcParams['lines.linewidth'] = 2

    ax.tick_params(axis='both', direction='in', which='major',  length=8, width=2)
    ax.tick_params(axis='both', direction='in', which='minor',  length=5, width=1.5)
    #ax.xaxis.set_major_locator(ticker.MultipleLocator(5))
    ax.xaxis.set_minor_locator(tick.MultipleLocator(5))
    ax.yaxis.set_minor_locator(tick.MultipleLocator(0.1))

