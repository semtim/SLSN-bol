import numpy as np
import pandas as pd
import os

temp = os.path.abspath("second_cut.csv")
name = pd.read_csv(temp, sep=",")
name = pd.DataFrame(name["Name"])
name = np.array( name["Name"].tolist() )

temp = os.path.abspath("table_of_sample_copy.csv")
table = pd.read_csv(temp, sep=",")

col = ['Name', 'redshift', 'ra', 'dec', 'bandpasses']
SLSNI, SLSNII, SLSNR = pd.DataFrame(columns=col), pd.DataFrame(columns=col), pd.DataFrame(columns=col)
for slsn in name:
    if table.loc[table["Name"] == slsn ]["type"].iloc[0] == 'SLSN-I':
        SLSNI = SLSNI.append(table.loc[table["Name"] == slsn ].drop(columns=['type']),
                             ignore_index=True)
    elif table.loc[table["Name"] == slsn ]["type"].iloc[0] == 'SLSN-II':
        SLSNII = SLSNII.append(table.loc[table["Name"] == slsn ].drop(columns=['type']),
                               ignore_index=True)
    elif table.loc[table["Name"] == slsn ]["type"].iloc[0] == 'SLSN-R':
        SLSNR = SLSNR.append(table.loc[table["Name"] == slsn ].drop(columns=['type']),
                             ignore_index=True)
SLSNI["bandpasses"] = '{' +SLSNI["bandpasses"] + '}'
SLSNII["bandpasses"] = '{' +SLSNII["bandpasses"] + '}'
SLSNR["bandpasses"] = '{' +SLSNR["bandpasses"] + '}'
SLSNI.to_csv('tables\\SLSNI.csv', index=False)
SLSNII.to_csv('tables\\SLSNII.csv', index=False)
SLSNR.to_csv('tables\\SLSNR.csv', index=False)