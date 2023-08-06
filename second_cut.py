import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import bb

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




second_cut = np.array(second_cut)
second_cut = second_cut[(second_cut != 'SN2017bcc') * (second_cut != 'SN2016ezh')] #SN2016ezh - TDE

pd.DataFrame(second_cut, columns=['Name']).to_csv('second_cut' + '.csv',
                   index=False)