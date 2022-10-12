import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize, least_squares
import os
import pandas as pd
import bb
from multiprocessing import Pool

temp = os.path.abspath("second_cut.csv")
name = pd.read_csv(temp, sep=",")
name = pd.DataFrame(name["Name"])
name = np.array( name["Name"].tolist() )

temp = os.path.abspath("redshifts.csv")
redshifts = pd.read_csv(temp, sep=",")
#cur_sample = ['SN2011ke', 'SN2007bi', 'SN2013dg', 'PTF10aagc', 'LSQ14bdq', 'PTF12mkp']
mask = (name == 'LSQ14mo')

def compute_bol(slsn):
    #slsn - name of object
    # data в магнитудах
    temp = os.path.abspath('approx_magn/data/' + slsn + ".csv")
    ap_data = pd.read_csv(temp, sep=",")
    z = float( redshifts.loc[redshifts["Name"] == slsn]['z'])

    columns = ['mjd', 'T', 'R', 'fun', 'success']
    BB_parameters = pd.DataFrame(columns=columns)
    BB_parameters["mjd"] = ap_data["mjd"]
    data_size = len(ap_data)
    bb_T, bb_R = [], []
    fun = []
    succ = []

    pred_params = np.array([1.7, 1])
    for d in range(data_size):
        result = minimize(bb.sum_squared_err, pred_params, args=(ap_data.iloc[d], z,),
                     method = 'Newton-CG',
                     jac=bb.band_mag_jac,
                     hess=bb.band_mag_hess
                     )
        bb_T.append(result.x[0] * bb.const[0])
        bb_R.append(result.x[1] * bb.const[1])
        succ.append(result.success)
        fun.append(result.fun)
        pred_params = (result.x[0], result.x[1])

    BB_parameters["T"], BB_parameters["R"] = bb_T, bb_R
    BB_parameters["fun"] = fun
    BB_parameters["success"] = succ

    BB_parameters.to_csv('bol_output/data/' + slsn + '.csv', index=False)

pool = Pool(processes=len(name[mask]))
pool.map(compute_bol, name[mask])
pool.close()
pool.join()
