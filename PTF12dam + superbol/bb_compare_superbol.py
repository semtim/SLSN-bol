import os
import pandas as pd
import matplotlib.pyplot as plt
from scipy import integrate
#from snad.load.curves import OSCCurve
import numpy as np
from sklearn.preprocessing import scale
from astropy.coordinates import Distance
import astropy.units as u

import sys
sys.path.append('D:\\гаиш\\SLSN_2021')
import bb

#print('SLSN name: \n')
name = 'PTF12dam' #input()

# print('save fig? [1/0]')
# save = int(input())
# if save:
#     os.mkdir(name + ' BB parameters')
path = os.path.abspath(name + ' BB parameters')

data_path = os.path.abspath(name + '_bb.csv')
data = pd.read_csv(data_path, sep=",")

data_sup_path = os.path.abspath("bol_" + name + "_gri.txt")
data_sup = pd.read_csv(data_sup_path, sep="\s+")

data_multicol_path = os.path.abspath(name + ".csv")
data_multicol = pd.read_csv(data_multicol_path, sep=",")

#raw_data = OSCCurve.from_name(name, bands=data_multicol.columns[1:],
 #                down_args={"baseurl": "https://sne.space/sne/"})

#raw_data = raw_data.filtered(with_upper_limits=False, with_inf_e_flux=False,
#                                                     sort='filtered')

#raw_data = raw_data.binned(bin_width=3, discrete_time=True)

####################################################################
M_sun, m_sun = 4.74, -26.74
L_sun = 3.828*1e33 #sgs
z = 0.1073
d_sun = 1.496 * 1e13 #cm
d = Distance(z=z, unit=u.cm).to_value()
d_pc = Distance(z=z, unit=u.pc).to_value()

data_size = len(data["T"])
L, M = [], []
for i in range(data_size):
    L.append( bb.flux(data["T"][i])*4*np.pi*data["R"][i]**2 * 1e+7 )
    M.append( M_sun - 2.5*np.log10(L[-1]/L_sun) )
L = np.array(L)

max_day = data["mjd"][np.argmax(L)]
mask = (data["mjd"] <= max_day + 150)
data = data[mask]
L = L[mask]

mask_mult = (data_multicol["mjd"] <= max_day + 150)
data_multicol = data_multicol[mask_mult]
data_multicol = data_multicol[(data_multicol["mjd"] >= data['mjd'][0])]
####################################################################
#model AB magnitudes
model_ABmagn = pd.DataFrame(columns=data_multicol.columns)
model_ABmagn["mjd"] = data['mjd']
for b in ['r','g','i','B']:
    temp = []
    for x in zip(data["T"]/bb.const[0], data["R"]/bb.const[1]):
        temp.append( bb.band_mag(b, x, z) )
    model_ABmagn[b] = temp


fig, ax = plt.subplots() #figsize=(18, 12),dpi=400 
for band in ['r','g','i','B']:
    plt.plot(data_multicol["mjd"], data_multicol[band], ls='--', alpha=0.3,
             color=bb.cols[band])
for band in ['r','g','i','B']:
    plt.plot(model_ABmagn["mjd"], model_ABmagn[band], color=bb.cols[band])
ax.invert_yaxis()

b = ['r','g','i','B']
plt.legend([band + '_obs' for band in b] + [band + '_model' for band in b],
            ncol=2, loc=1)
plt.ylabel('$m$')
plt.xlabel('mjd')
plt.xlim(data['mjd'][0] - 10, data['mjd'][len(data['mjd'])-1] + 10)
plt.ylim(20, 16.5)
if save:
    #os.mkdir(name + ' figures')
    path_fig = os.path.abspath(name + ' figures')
    fig.savefig( fname = path_fig + '/' + name + '_ABmagn'  ,
                bbox_inches="tight")
############
#hat band
# ABmagn_hat = []
# for x in zip(data["T"], data["R"]):
#     ABmagn_hat.append( bb.band_mag('hat', x, z) )
#####################################################################

#magn
# bb_magn = m_sun - 2.5 * np.log10(L / L_sun) + 5 * np.log10( d / d_sun)
# fig, ax = plt.subplots(figsize=(18, 12),dpi=400 )
# plt.plot(data["mjd"], bb_magn, c="black")

# for band in data_multicol.columns[1:-1]:
#     plt.plot(data_multicol["mjd"], data_multicol[band], ls='--', alpha=0.3, 
#              c=bb.cols[band])

# #plt.plot(data_multicol["mjd"], sum_magn, c = 'purple')
# plt.xlim(min(raw_data.X[:,1]), max(data_multicol["mjd"]))
# ax.invert_yaxis()

# plt.legend(["bb fit"] + 
#                data_multicol.columns[1:-1].to_list()
#                + ["pseudobolometric"], ncol=2, loc=1)
# plt.title("Apparent magnitudes")
###########################
#absolute magnitudes
M_sup = M_sun - 2.5*np.log10( np.array(data_sup["L_BB"].to_list()) /L_sun)
M_col = []
for band in data_multicol.columns[1:-1]:
    m = np.array( data_multicol[band].tolist() )
    M_col.append( m - 5*np.log10(d_pc/10) )

fig, ax = plt.subplots()#figsize=(18, 12),dpi=400 
plt.plot(data["mjd"], M)
plt.xlabel('mjd')
plt.ylabel('Absolute magnitude')
ax.invert_yaxis()
for Mb, band in zip(M_col, b[:-1]):
    plt.plot(data["mjd"], Mb, ls='--', alpha=0.3, c=bb.cols[band])

plt.errorbar(data_sup["mjd"], M_sup, fmt='o')
plt.legend(["bb fit"] + data_multicol.columns[1:-1].to_list() + 
           ["superbol fit"], ncol=2, loc=1)

if save:
    fig.savefig( fname = path_fig + '\\' + name + '_abs_magn' ,
                bbox_inches="tight")
    
    
###########################
#Luminosity

mjd_obs, L_obs = [-60, -42, -5, 10, 26, 60, 90, 175], [43.25, 43.67, 44.166, 44.14, 44.03,43.69, 43.5, 43]#from DOI: 1310.4446
day_max = data["mjd"][np.argmax(L)]

pseudobol = np.zeros(len(data_multicol))
for band in ['r','g','i']:
    pseudobol += bb.ABmagn_to_flux(data_multicol[band])* \
        bb.width_nu(band) * 4*np.pi*d**2
# for band in data_multicol.columns[1:-1]:
#     pseudobol += bb.magn_to_flux(data_multicol[band])* \
#         bb.width[band] * bb.fref(band) * 4*np.pi*d**2      
    

fig, ax = plt.subplots() #figsize=(18, 12),dpi=400
plt.plot(data["mjd"], np.log10(L))
plt.scatter(mjd_obs + day_max, L_obs, c='r')
plt.plot(data_multicol["mjd"], np.log10(pseudobol))
plt.errorbar(data_sup["mjd"], np.log10(data_sup["L_BB"]),
           data_sup["err_BB"]/data_sup["L_BB"]/np.log(10), fmt='o')
plt.xlabel('mjd')
plt.ylabel('$log_{10}[ L, erg/s ]$')
plt.legend(['BB fit','pseudobolometric (gri sum)', 'M. Nicholl et al. 2013',
             'superbol fit'])
plt.xlim(data['mjd'][0] - 10, data['mjd'][len(data['mjd'])-1] + 10)

if save:
    fig.savefig( fname = path_fig + '/' + name + '_luminosity' ,
                bbox_inches="tight")

#########################

#Radius
fig, ax = plt.subplots(figsize=(18, 12),dpi=400)
#plt.subplot(122)
R = np.array(data["R"]) * 100 #cm
plt.plot(data["mjd"], R)
plt.xlabel("MJD")
plt.ylabel('R, cm')
#plt.plot(data["mjd"][9:], moving_average(10, data["R"].tolist()) )
#plt.legend(['model R', 'moving average (n=10)'])
plt.title('BB parameters')
if save:
    fig.savefig( fname = path + '\\' + name + '_radius_minimize'  ,
                bbox_inches="tight")
#########################

#Temperature

T_obs = np.array([1.325, 1.625, 1.8, 1.35, 0.95, 0.8]) * 1e4
T_obs_err =  np.array([0.4, 0.4, 0.3, 0.2, 0.15, 0.07]) * 1e4
T_days = [-50, -40, -12.5, 10, 55, 92]
fig, ax = plt.subplots(figsize=(18, 12),dpi=400)#figsize=(18, 12),dpi=400
#plt.subplot(121)
plt.plot(data["mjd"][1:], data["T"][1:])
plt.xlabel('MJD')
plt.ylabel('Temperature, K')
plt.title('BB parameters')

#plt.errorbar(T_days, T_obs, yerr=T_obs_err, c='r', fmt='o')
plt.legend(['BB temperature', 'photometry (M. Nicholl et al. 2013)'])
#T moving average
#plt.plot(data["mjd"][9:], moving_average(10, data["T"].tolist()) )
#plt.legend(['model T', 'moving average (n=10)'])
if save:
    fig.savefig( fname = path + '\\' + name + '_temperature'  
                , bbox_inches="tight")
##########################################################
#sum squared errors
fig, ax = plt.subplots()
plt.plot(data["mjd"], data["fun"])
plt.xlabel('mjd')
plt.ylabel('sum squared errors')
if save:
    fig.savefig( fname = path_fig + '/' + name + '_sse'  
                , bbox_inches="tight")
##########################################################
#saving .dat file for superbol
# os.mkdir(name + ' raw magn')
# path_raw = os.path.abspath(name + ' raw magn')
# band_ledge = 0
# for t, band in enumerate(data_multicol.columns[1:]):
#     raw_magn, magn_err = bb.flux_to_magn( raw_data.odict[band].y, 
#                                         raw_data.odict[band].err)
#     mask =(raw_data.X[:,0] == t) * (raw_data.X[:,1] < np.max(data["mjd"]) - 60) \
#         * (raw_data.X[:,1] > np.min(data["mjd"]))
#     mjd = raw_data.X[mask][:,1]
#     band_size = sum((raw_data.X[:,0] == t))
#     raw_magn = raw_magn[mask[band_ledge:band_ledge + band_size]]
#     magn_err = magn_err[mask[band_ledge:band_ledge + band_size]]
#     band_ledge += band_size
#     columns = ['mjd', band, band + '_err']
    
#     save_raw = pd.DataFrame(columns=columns)
#     save_raw["mjd"], save_raw[band] = mjd, raw_magn
#     save_raw[band + '_err'] = magn_err #+ raw_magn / 50
#     save_raw.to_csv( path_raw + '\\' + name + '_'+ band + '.dat', index=False)
    