import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from snad.load.curves import OSCCurve
import os


temp = os.path.abspath("SLSN.csv")
name = pd.read_csv(temp, sep=",")
name = pd.DataFrame(name)
name = name["Name"]

sn=[]
# for i in range(len(name)):
#     sn.append( OSCCurve.from_name(name[i], 
#                 down_args={"baseurl": "https://sne.space/sne/"}) )

name = name.drop(59)#недостаточно данных фотометрии
name = name.drop(60) 
name = name.drop(61)
name = name.drop(62)

name = np.array(name.tolist())
mask = (name != 'MLS150612:112519+081418') * \
    (name != 'CSS121015:004244+132827') * (name != 'SSS120810:231802-560926') *\
(name != 'OGLE15sd') * (name != 'SN2009jh') * (name != 'SN2010hy') *\
	(name != 'SN2016ard') * (name != 'SN2016els')
name = name[mask]

for i in range(len(name)):
    try:
        sn.append( OSCCurve.from_json( os.path.join('./sne', name[i] + '.json') ))
        sn[-1] = sn[-1].filtered(with_upper_limits=False, with_inf_e_flux=False,
                                sort='filtered')
    except FileNotFoundError:
        continue

filt_name = []
for i in range(len(sn)):
    if len(sn[i].bands) >= 3:
        filt_name.append(sn[i].name)
filt_name = pd.DataFrame(filt_name, columns=['Name'])

filt_name.to_csv('sample_SLSN.csv', index=False)
