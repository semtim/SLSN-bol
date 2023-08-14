import matplotlib.font_manager as font_manager
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.ticker as ticker

date = [2015, 2015.5, 2016, 2016.5, 2017, 2017.5, 2018, 2018.5, 2019,
        2019.5, 2020, 2020.5, 2021, 2021.5, 2022, 2022.5]
# classified SN from TNS
# allSN = [5061, 5095, 5148, 5427, 5788, 6116, 6535, 7117, 8037, 8979, 10182, 11107, 12292, 13473, 14629, 15578]
# slsn = [0, 1, 1, 5, 7, 10, 12, 21, 45, 77, 100, 111, 134, 153, 161, 173]


# spectra observation num > 0
allSN = [95, 104, 124, 419, 826, 1184, 1687, 2350, 3413, 4439, 5734, 6732, 7980, 9239, 10506, 11547]
slsn = [0, 1, 1, 5, 7, 10, 12, 21, 45, 77, 100, 111, 134, 153, 161, 172]

plt.rc('text', usetex=True)
font = {'family' : 'Times New Roman',
    #'weight' : 'bold',
    'size'   : 22}

plt.rc('font', **font)
plt.rcParams['axes.linewidth'] = 1.2

fig, ax = plt.subplots(figsize=(10, 7))
#plt.title( 'Cumulative distribution of SNe' )

ax.step(date, allSN, linewidth=2, c='b')
ax.step(date, slsn, linewidth=2, c='red')
lab = ['All SNe', 'SLSNe']
ax.legend(lab, loc='best', columnspacing=0.7, labelspacing=0.3,
               handletextpad=0.35,)

ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
#ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.5))


ax.tick_params(axis='both', direction='in', which='major',  length=8, width=2)
ax.tick_params(axis='both', direction='in', which='minor',  length=5, width=1.5)

ax2 = ax.twinx()
ax2.step(date, allSN, linewidth=2, visible=False)
ax2.step(date, slsn, linewidth=2, visible=False)
ax2.tick_params(axis='both', direction='in', which='major',  length=8, width=2)
ax2.tick_params(axis='both', direction='in', which='minor',  length=5, width=1.5)
ax2.semilogy()


ax.grid('on', linestyle='--', alpha=0.7, linewidth=1)

ax.set_ylabel('Cumulative number')
ax.set_xlabel('Year')

ax.semilogy()


fig.savefig( fname = 'fig_for_article/cumulative_corrected.pdf',
                    bbox_inches="tight", format='pdf')