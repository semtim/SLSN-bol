import numpy as np
import pandas as pd
import bb
import os
import matplotlib.pyplot as plt

temp = os.path.abspath("second_cut.csv")
name = pd.read_csv(temp, sep=",")
name = pd.DataFrame(name["Name"])
name = np.array( name["Name"].tolist() )

temp = os.path.abspath("table_of_sample_copy.csv")
table = pd.read_csv(temp, sep=",")


fig1, ax1 = plt.subplots(figsize=(10, 7)) #figsize=(18, 12),dpi=400
ax1.set_xlabel('mjd')
ax1.set_ylabel('$log_{10}[ L, erg/s ]$')
#ax1.set_title('SLSN-I')
    
fig2, ax2 = plt.subplots(figsize=(10, 7)) #figsize=(18, 12),dpi=400
ax2.set_xlabel('mjd')
ax2.set_ylabel('$log_{10}[ L, erg/s ]$')
ax2.set_title('SLSN-II')
    
figr, axr = plt.subplots(figsize=(10, 7)) #figsize=(18, 12),dpi=400
axr.set_xlabel('mjd')
axr.set_ylabel('$log_{10}[ L, erg/s ]$')
axr.set_title('SLSN-R')

mark1 = ['o', 'v', '^', '<', '>', 's',"P", "*",
         'D', '$\cap$', '$\cup$',]

mark2 = ["1", "2", "3", '4','+', 'x', "|", "_",
         '$\otimes$', '$\dagger$', '$wr$', '$\ominus$', '$\oslash$', 
         '$\odot$']

line = ['-', '--', '-.', ':']
c1, c1_l, c1_2, c2, cr = 0, 0, 0, 0, 0
name1, name2, namer = [], [], []

for slsn in name:
    data_path = os.path.abspath('D:\\гаиш\\SLSN_2021\\bol_output_magn\\' + slsn + ".csv")
    data = pd.read_csv(data_path, sep=",")
    data_size = len(data["T"])
    L = []
    for i in range(data_size):
        L.append( bb.flux(data["T"][i])*4*np.pi*data["R"][i]**2 * 1e+7 )
    L = np.array(L)

    
    if table.loc[table['Name'] == slsn]["type"].iloc[0][:7] == 'SLSN-II':
        z = table.loc[table['Name'] == slsn]["redshift"].iloc[0]
        mjd = np.array(data["mjd"]) / (z + 1)
        mjd = mjd - mjd[np.argmax(L)]
        ax2.scatter(mjd, np.log10(L), marker=mark1[c2], s=30)
        c2 += 1
        name2.append(slsn)
    
    elif table.loc[table['Name'] == slsn]["type"].iloc[0][:6] == 'SLSN-I':
        z = table.loc[table['Name'] == slsn]["redshift"].iloc[0]
        mjd = np.array(data["mjd"]) / (z + 1)
        mjd = mjd - mjd[np.argmax(L)]
        if c1 < 11:
            ax1.scatter(mjd, np.log10(L), marker=mark1[c1], s=25,
                        facecolors='none',
                        edgecolors=list(bb.cols.values())[c1])
            c1 += 1
            name1.append(slsn)
        
        elif c1_2 < 10:
            ax1.scatter(mjd, np.log10(L), marker=mark2[c1_2], s=25,
                        c=list(bb.cols.values())[c1+c1_2])
            c1_2 += 1
            name1.append(slsn)
        else:
            ax1.plot(mjd, np.log10(L), linestyle=line[c1_l])
            c1_l += 1
            name1.append(slsn)
            
        # ax1.plot(mjd, np.log10(L), c=list(bb.cols.values())[c1])
        # ax1.text(min(mjd[-1], 140), np.log10(L)[-1] + 0.01, slsn,
        #          c=list(bb.cols.values())[c1])
        # c1 += 1
        
    elif table.loc[table['Name'] == slsn]["type"].iloc[0][:6] == 'SLSN-R':
        z = table.loc[table['Name'] == slsn]["redshift"].iloc[0]
        mjd = np.array(data["mjd"]) / (z + 1)
        mjd = mjd - mjd[np.argmax(L)]
        axr.scatter(mjd, np.log10(L), marker=mark1[cr], s=30)
        cr += 1
        namer.append(slsn)

ax1.set_xlim(-50, 160)
ax1.set_ylim(42.25, 44.8)
ax1.legend(name1, ncol=3, columnspacing=0.8, labelspacing=0.3,
           handletextpad=0.4)
ax2.legend(name2)
axr.legend(namer)

fig1.savefig( fname='D:\\гаиш\\SLSN_2021\\types_figures\\' + 
                    'SLSN-I', bbox_inches="tight")
fig2.savefig( fname='D:\\гаиш\\SLSN_2021\\types_figures\\' + 
                    'SLSN-II', bbox_inches="tight")
figr.savefig( fname='D:\\гаиш\\SLSN_2021\\types_figures\\' + 
                    'SLSN-R', bbox_inches="tight")

plt.close()
plt.close()
plt.close()
#######################################
#PTF12dam SN2007bi
# figcomp, axcomp = plt.subplots(figsize=(18, 12),dpi=400) #figsize=(18, 12),dpi=400
# axcomp.set_xlabel('mjd')
# axcomp.set_ylabel('$log_{10}[ L, erg/s ]$')

# for slsn in ['SN2007bi', 'PTF12dam']:
#     data_path = os.path.abspath('D:\\гаиш\\SLSN_2021\\bol_output_2cut\\' + slsn + ".csv")
#     data = pd.read_csv(data_path, sep=",")
#     data_size = len(data["T"])
#     L = []
#     for i in range(data_size):
#         L.append( bb.flux(data["T"][i])*4*np.pi*data["R"][i]**2 * 1e+7 )
#     L = np.array(L)
    
#     z = table.loc[table['Name'] == slsn]["redshift"].iloc[0]
#     mjd = np.array(data["mjd"]) / (z + 1)
#     mjd = mjd - mjd[np.argmax(L)]
#     axcomp.plot(mjd, np.log10(L))
# plt.legend(['SN2007bi', 'PTF12dam'])

fig1, ax1 = plt.subplots(figsize=(10, 7)) #figsize=(18, 12),dpi=400
ax1.set_xlabel('mjd (rest frame)')
ax1.set_ylabel('$log_{10}[ L, erg/s ]$')
#ax1.set_title('SLSN-I')
c1, c2, c1_l, c1_2 = 0, 0, 0, 0
for slsn in name:
    data_path = os.path.abspath('D:\\гаиш\\SLSN_2021\\bol_output_magn\\' + slsn + ".csv")
    data = pd.read_csv(data_path, sep=",")
    data_size = len(data["T"])
    L = []
    for i in range(data_size):
        L.append( bb.flux(data["T"][i])*4*np.pi*data["R"][i]**2 * 1e+7 )
    L = np.array(L)
    
    z = table.loc[table['Name'] == slsn]["redshift"].iloc[0]
    mjd = np.array(data["mjd"]) / (z + 1)
    mjd = mjd - mjd[np.argmax(L)]
    if c1 < 11:
        ax1.scatter(mjd, np.log10(L), marker=mark1[c1], s=30,
                    facecolors='none',
                    edgecolors=list(bb.cols.values())[c1])
        c1 += 1
        
    elif c1_2 < 14:
        ax1.scatter(mjd, np.log10(L), marker=mark2[c1_2], s=30,
                    c=list(bb.cols.values())[c1+c1_2])
        c1_2 += 1
        
    else:
        ax1.plot(mjd, np.log10(L), linestyle=line[c1_l])
        c1_l += 1
        
ax1.set_xlim(-50, 166)
ax1.set_ylim(42.25, 44.95)
ax1.legend(name, ncol=3, columnspacing=0.7, labelspacing=0.25,
           handletextpad=0.35)

fig1.savefig( fname='D:\\гаиш\\SLSN_2021\\types_figures\\' + 
                    'ALL', bbox_inches="tight")
