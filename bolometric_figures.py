import os
import pandas as pd
import matplotlib.pyplot as plt
from scipy import integrate
import numpy as np
from sklearn.preprocessing import scale
from astropy.coordinates import Distance
import astropy.units as u
import bb

temp = os.path.abspath("second_cut.csv")
name = pd.read_csv(temp, sep=",")
name = pd.DataFrame(name["Name"])
name = np.array( name["Name"].tolist() )

temp = os.path.abspath("redshifts.csv")
redshifts = pd.read_csv(temp, sep=",")

for slsn in name:
    try:
        data_path = os.path.abspath('bol_output/data/'
                                    + slsn + ".csv")
        data = pd.read_csv(data_path, sep=",")
        z = redshifts.loc[redshifts['Name'] == slsn]["z"].iloc[0]
        #z = 0.003319  #SN2018aoq
        bb.plot_sub(slsn, z, data, save=1)
    except:
        continue
# fig, axs = plt.subplots(9, 5, figsize=(30, 60)) #figsize=(18, 12),dpi=400
# #plt.subplots_adjust(wspace=0.2, hspace=0.2)
# for n, slsn in enumerate(name):
#     data_path = os.path.abspath(slsn + ".csv")
#     data = pd.read_csv(data_path, sep=",")

#     data_size = len(data["T"])
#     L, M = [], []
#     for i in range(data_size):
#         L.append( bb.flux(data["T"][i])*4*np.pi*data["R"][i]**2 * 1e+7 )
#         M.append( M_sun - 2.5*np.log10(L[-1]/L_sun) )
#     L = np.array(L)
    
#     axs[n // 5 , n % 5].plot(data["mjd"], np.log10(L))
#     axs[n // 5 , n % 5].set_title(slsn)
#     #axs[i].xlabel('mjd')
#     #plt.ylabel('$log_{10}[ L, erg/s ]$')
# fig.savefig( fname='bolometric_fig', bbox_inches="tight")

# for slsn in name:
#     try:
#         data_path = os.path.abspath('D:\\гаиш\\SLSN_2021\\bol_output_magn\\' + slsn + ".csv")
#         data = pd.read_csv(data_path, sep=",")
        
#         # fig, ax = plt.subplots() #figsize=(18, 12),dpi=400
#         # plt.plot(data["mjd"], data["fun"])
#         # plt.xlabel('mjd')
#         # plt.ylabel('sum squared errors')
#         # plt.title(slsn)
#         # fig.savefig( fname='D:\\гаиш\\SLSN_2021\\bol_output_magn\\figures\\sse\\' +
#         #             slsn, bbox_inches="tight")
        

#         # fig, ax = plt.subplots() #figsize=(18, 12),dpi=400
#         # plt.plot(data["mjd"], data["R"])
#         # plt.xlabel('mjd')
#         # plt.ylabel('Radius')
#         # plt.title(slsn)
#         # fig.savefig( fname='D:\\гаиш\\SLSN_2021\\bol_output_magn\\figures\\Radius\\' + 
#         #             slsn, bbox_inches="tight")
        
#         fig, ax = plt.subplots() #figsize=(18, 12),dpi=400
#         plt.plot(data["mjd"], data["T"])
#         plt.xlabel('mjd')
#         plt.ylabel('Temperature')
#         plt.title(slsn)
#         fig.savefig( fname='D:\\гаиш\\SLSN_2021\\bol_output_magn\\figures\\Temperature\\' + 
#                     slsn, bbox_inches="tight")
#     except:
#         continue