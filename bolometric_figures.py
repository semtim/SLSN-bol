import os
import pandas as pd
import matplotlib.pyplot as plt
from scipy import integrate
#from snad.load.curves import OSCCurve
import numpy as np
from sklearn.preprocessing import scale
from astropy.coordinates import Distance
import astropy.units as u

#import sys
#sys.path.append('/media/documents/гаиш/SLSN_2021')
import bb

#temp = os.path.abspath("/media/documents/гаиш/SLSN_2021/sample_SLSN.csv")
temp = os.path.abspath("second_cut.csv")
name = pd.read_csv(temp, sep=",")
name = pd.DataFrame(name["Name"])
name = np.array( name["Name"].tolist() )

temp = os.path.abspath("redshifts.csv")
redshifts = pd.read_csv(temp, sep=",")

for slsn in name:
    try:
        data_path = os.path.abspath('bol_output_fixed/data/'
                                    + slsn + ".csv")
        data = pd.read_csv(data_path, sep=",")
        z = redshifts.loc[redshifts['Name'] == slsn]["z"].iloc[0]
        #z = 0.003319  #SN2018aoq
        bb.plot_sub(slsn, z, data, save=1)
    except:
        continue
