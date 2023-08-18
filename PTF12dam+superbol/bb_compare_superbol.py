import os
import pandas as pd
import matplotlib.pyplot as plt
from scipy import integrate
#from snad.load.curves import OSCCurve
import numpy as np
from sklearn.preprocessing import scale
from astropy.coordinates import Distance
import astropy.units as u
import matplotlib.ticker as tick
from matplotlib import rc

import sys
sys.path.append('/home/timofey/saimsu/SLSN-bol')
import bb


data_path = os.path.abspath('PTF12dam_bb.csv')
data = pd.read_csv(data_path, sep=",")

griB_UV_JHK = pd.read_csv('griBUVJHK.csv', sep=",")
griB_JHK = pd.read_csv('griB_JHK.csv', sep=",")


data_sup_path = os.path.abspath("bol_PTF12dam_Bgri.txt")
data_sup = pd.read_csv(data_sup_path, sep="\s+")

data_multicol_path = os.path.abspath("PTF12dam_multicol.csv")
data_multicol = pd.read_csv(data_multicol_path, sep=",")

####################################################################
z = 0.1073

#model AB magnitudes and sse
model_ABmagn = pd.DataFrame(columns=data_multicol.columns)
model_ABmagn["mjd"] = data['mjd']
for b in ['r','g','i','B']:
    temp = []
    for x in zip(data["T"]/bb.const[0], data["R"]/bb.const[1]):
        temp.append( bb.band_mag(b, x, z) )
    model_ABmagn[b] = temp


fig, axs = plt.subplots(2, sharex=True, figsize=(10, 14))
for band, band_it in zip(['r','g','i','B'], [r'\textit{r}', r'\textit{g}', r'\textit{i}', r'\textit{B}']):
    axs[0].plot(data_multicol["mjd"], data_multicol[band], ls='--', alpha=0.3,
             color=bb.cols[band], linewidth=2, label=band_it)
for band, band_it in zip(['r','g','i','B'], [r'\textit{r}', r'\textit{g}', r'\textit{i}', r'\textit{B}']):
    axs[0].plot(model_ABmagn["mjd"], model_ABmagn[band], color=bb.cols[band],
                linewidth=2, label=band_it)
axs[0].invert_yaxis()


#get handles and labels
handles, labels = axs[0].get_legend_handles_labels()
#specify order of items in legend
order = [3,7,1,5,0,4,2,6]
#add legend to plot
axs[0].legend([handles[idx] for idx in order],[labels[idx] for idx in order],
            ncol=4, loc='lower right', columnspacing=0.7, labelspacing=0.3,
                       handletextpad=0.35)
axs[0].text(56078.5, 19.5,r'\textbf{observed:}')
axs[0].text(56090, 19.77,r'\textbf{model:}')

axs[0].set_ylabel('Apparent magnitude')
axs[0].set_xlim(data['mjd'][0] - 10, data['mjd'][len(data['mjd'])-1] + 10)
axs[0].set_ylim(20, 16.5)

axs[1].plot(data["mjd"], np.sqrt(data["fun"] / 4), c='blue')
axs[1].set_xlabel('MJD')
axs[1].set_ylabel('RMSE')
plt.subplots_adjust(hspace=0.05)

for ax in axs:
    ax.label_outer()
    bb.presets_fig(ax)

for ax in axs:
    ax.grid('on', linestyle='--', alpha=0.7, linewidth=1)

axs[0].yaxis.set_minor_locator(tick.MultipleLocator(0.25))
axs[0].yaxis.set_major_locator(tick.MultipleLocator(0.5))
yticks = [16.5, 17, 17.5, 18, 18.5, 19, 19.5]  #axs[0].get_yticklabels()
axs[0].set_yticks(yticks)


#axs[1].yaxis.set_major_locator(tick.MultipleLocator(0.001))
#axs[1].yaxis.set_minor_locator(tick.MultipleLocator(0.0005))
axs[1].xaxis.set_major_locator(tick.MultipleLocator(30))


path_fig = os.path.abspath('PTF12dam figures')
fig.savefig( fname = path_fig + '/PTF12dam_BBaccur.pdf'  
                , bbox_inches="tight", format='pdf', dpi=1000)
    
###########################
#Luminosity
#Matt Nicholl et al. 2013
mjd_obs, L_obs = [-60, -42, -5, 10, 26, 60, 90, 175], [43.25, 43.67, 44.166, 44.14, 44.03,43.69, 43.5, 43]#from DOI: 1310.4446
mjd = data["mjd"] / (1 + z)
logL = data['logL'] + 7
day_max = mjd[np.argmax(logL)]

#Vreeswijk et al. 2016
mjd_vreeswijk = np.array([-78.86, -75.24, -71.68, -68.86, -66.28, -25.71, -23.76, -19.69, -17.06, -14.89,
                 -13.22, -10.37, -6.92, -5.87, -4.54, 1.44, 5.89, 8.45, 9.69, 10.39, 14.42, 19.26,
                 29.34, 40.21, 58.20, 72.52, 79.85, 105.98, 134.35, 140.65, 155.05, 168.56, 181.16, 222.56,
                 224.44, 231.57, 239.73, 246.99, 250.37, 257.69, 281.11, 334.29]) + 56096.7/ (1 + z)
L_vreeswijk = [42.49, 42.54, 42.51, 42.58, 43.07, 44.21, 44.24, 44.25, 44.28, 44.29, 44.29, 
               44.31, 44.31, 44.30, 44.30, 44.28, 44.25, 44.23, 44.23, 44.22, 44.19, 44.16,
               44.07, 43.98, 43.86, 43.78, 43.72, 43.55, 43.38, 43.34, 43.24, 43.14, 43.06, 
               42.81, 42.80, 42.76, 42.71, 42.67, 42.69, 42.61, 42.47, 42.12]

err_vreeswijk = [0, 0, 0, 0, 0, 0.07, 0.06, 0.09, 0.1, 0.11, 0.1, 0.09, 0.09, 0.09, 0.08,
                 0.08, 0.08, 0.08, 0.08, 0.08, 0.07, 0.09, 0.08, 0.1, 0.09, 0.1, 0.09, 0.09,
                 0.07, 0.1, 0.11, 0.13, 0.11, 0.11, 0.11, 0.12, 0.12, 0.09, 0.08, 
                 0.11, 0.15, 0.17]

#pseudobol = np.zeros(len(data_multicol))
#for band in ['r','g','i']:
#    pseudobol += bb.ABmagn_to_flux(data_multicol[band])* \
#        bb.width_nu(band) * 4*np.pi*d**2
#pseudobol *= (1 + z) # redshift correction


superbol_err = data_sup["err_BB"]/data_sup["L_BB"]/np.log(10)
mask = superbol_err < 0.47


fig, ax = plt.subplots(figsize=(10, 7)) #figsize=(18, 12),dpi=400
bb.presets_fig(ax)

ax.plot(mjd, logL, c='b', label=r'vector GP (\textit{Bgri}) + BB fit')
ax.fill_between(mjd, logL - data['sigma_logL'], logL + data['sigma_logL'],
                alpha=0.3, color='b')


ax.plot(griB_UV_JHK['mjd'] / (1+z),
        griB_UV_JHK['logL'] + 7,
        c='r', label=r'vector GP (\textit{UV}+\textit{Bgri}+\textit{JHK}) + BB fit')
ax.fill_between(griB_UV_JHK['mjd'] / (1+z),
                griB_UV_JHK['logL'] + 7 - griB_UV_JHK['sigma_logL'],
                griB_UV_JHK['logL'] + 7 + griB_UV_JHK['sigma_logL'],
                alpha=0.3, color='r')

ax.plot(griB_JHK['mjd'] / (1+z),
        griB_JHK['logL'] + 7,
        c='green', label=r'vector GP (\textit{Bgri}+\textit{JHK}) + BB fit')
ax.fill_between(griB_JHK['mjd'] / (1+z),
                griB_JHK['logL'] + 7 - griB_JHK['sigma_logL'],
                griB_JHK['logL'] + 7 + griB_JHK['sigma_logL'],
                alpha=0.3, color='green')



#ax.plot(data_multicol["mjd"]-data_multicol["mjd"][np.argmax(pseudobol)],
#        np.log10(pseudobol), c='limegreen', label="$gri$ sum")

ax.errorbar(data_sup["mjd"], #data_sup["mjd"][np.argmax(data_sup["L_BB"])]
            np.log10(data_sup["L_BB"]),
           superbol_err, fmt='o', ms=8, c='black',
           label='SuperBol', alpha=0.7, capsize=3)

#ax.errorbar(mjd_obs, L_obs, ms=8, c='r', label='Nicholl et al. 2013', marker='v', fmt='o')
ax.errorbar(mjd_vreeswijk, L_vreeswijk,
           err_vreeswijk, fmt='o', ms=8, c='darkorange',
           label='Vreeswijk et al. 2017', alpha=0.7, capsize=3)

ax.set_xlabel(r'Rest-frame time (days)')
ax.set_ylabel('$log_{10}[ L, erg/s ]$')
ax.set_title('PTF12dam')
ax.legend(ncol=1,
          columnspacing=0.4,
          labelspacing=0.2,
          handletextpad=0.35, fontsize=17,
          bbox_to_anchor=(1.01, 1.015),
          loc="upper right",
          )


plt.xlim(mjd[0] - 7, mjd[len(data['mjd']) - 1] + 5)
plt.ylim(43.2, 44.9)
ax.xaxis.set_major_locator(tick.MultipleLocator(30))
ax.grid('on', linestyle='--', alpha=0.7, linewidth=1)
ax.yaxis.set_minor_locator(tick.MultipleLocator(0.125))


fig.savefig( fname = 'PTF12dam figures/PTF12dam_luminosity.pdf' ,
                bbox_inches="tight", format='pdf')
