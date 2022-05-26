import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize, least_squares
import os
import pandas as pd
import bb

# data в магнитудах       
temp = os.path.abspath("SN2018aoq\\approx_results\\data\\SN2018aoq.csv")#approx_magn_results\\data\\
ap_data = pd.read_csv(temp, sep=",")
z = 0.003319
#z = 0.1073

columns = ['mjd', 'T', 'R', 'fun', 'success']
BB_parameters = pd.DataFrame(columns=columns)
BB_parameters["mjd"] = ap_data["mjd"]
data_size = len(ap_data)
bb_T, bb_R = [], []
fun = []
succ = []
# T, R = fsolve(mse_init, (10000, 1e+13), args=(ap_data.iloc[0], z,))
# bb_T.append(T)
# bb_R.append(R)

pred_params = np.array([1.61096187, 0.27153882])
for d in range(data_size):
    result = minimize(bb.sum_squared_err, pred_params,
                     args=(ap_data.iloc[d], z,),
                     #method = 'CG',
                     #bounds=([3000,80000],[1e+10, 1e+18]),
                     method = 'Newton-CG',
                     jac=bb.band_mag_jac,
                     hess=bb.band_mag_hess
                     #tol=0.1
                     )
    bb_T.append(result.x[0] * bb.const[0])
    bb_R.append(result.x[1] * bb.const[1])
    succ.append(result.success)
    fun.append(result.fun)
    pred_params = (result.x[0], result.x[1])

BB_parameters["T"], BB_parameters["R"] = bb_T, bb_R
BB_parameters["fun"] = fun
BB_parameters["success"] = succ

BB_parameters.to_csv( 'PTF12dam_bolometric_output' + '.csv', index=False)




# fig, ax = plt.subplots(figsize=(18, 12),dpi=400) #figsize=(18, 12),dpi=400
# plt.plot(BB_parameters["mjd"], BB_parameters["T"])
# plt.xlabel('MJD')
# plt.ylabel('Temperature, K')

# L = []
# for i in range(data_size):
#     L.append( bb.flux(bb_T[i])*4*np.pi*bb_R[i]**2 )
    
# fig, ax = plt.subplots(figsize=(18, 12),dpi=400) 
# plt.plot(BB_parameters["mjd"], L)
# plt.xlabel('MJD')
# plt.ylabel('Luminousity')