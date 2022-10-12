import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import bb

#temp = os.path.abspath("/media/documents/гаиш/SLSN_2021/sample_SLSN.csv")
temp = os.path.abspath("sample_SLSN.csv")
name = pd.read_csv(temp, sep=",")
name = pd.DataFrame(name["Name"])
name = np.array( name["Name"].tolist() )

second_cut = []
for slsn in name:
    data_path = os.path.abspath('approx_results/data/' + slsn + ".csv")
    data = pd.read_csv(data_path, sep=",")
    
    n = len(data.columns[1:])
    cut = False
    for i in range(0, n, 2):
        cur_arg = np.argmin(data[data.columns[1:][i]]) #смотрим по каждому фильтру, где максимум
        #print(data.columns[1:][i])
        if cur_arg < 2 or cur_arg > len(data) - 3: #если максимум в первых 6 или последних 6 днях, то выбрасываем кривую
            cut = True
            break
    if not cut:
        second_cut.append(slsn)

# for slsn in second_cut:
#     data_path = os.path.abspath('approx_results/data/' + slsn + ".csv")
#     data = pd.read_csv(data_path, sep=",")
  
#     min_mjd = data['mjd'][len(data) - 1]
#     n = len(data.columns[1:])
#     for i in range(0, n, 2):
#         tmp = np.array(bb.magn_to_flux(data[data.columns[1:][i]]))
#         if not tmp[tmp<=0].size:
#             print(i)
#             min_magn = bb.flux_to_magn(tmp[tmp<=0][0])
#             min_mjd = min(min_mjd, data['mjd'].loc[data.columns[1:][i] == min_magn])
second_cut = np.array(second_cut)
second_cut = second_cut[(second_cut != 'SN2017bcc') * (second_cut != 'SN2016ezh')] #SN2016ezh - TDE

pd.DataFrame(second_cut, columns=['Name']).to_csv('second_cut' + '.csv',
                   index=False)