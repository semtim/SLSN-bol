import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from multistate_kernel import MultiStateKernel
from snad.load.curves import OSCCurve
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel,Matern,ConstantKernel,RationalQuadratic
import os
import math
import warnings
from sklearn.exceptions import ConvergenceWarning
import time
from scipy import stats
import bb
import matplotlib.ticker as ticker

##########
temp = os.path.abspath("second_cut.csv")
name = pd.read_csv(temp, sep=",")
name = pd.DataFrame(name)

res=3
sn=[[] for i in range(len(name))]


for i in range(len(name)):
    sn[i] = OSCCurve.from_json(os.path.join('./sne', name['Name'][i] + '.json'))
    sn[i] = sn[i].filtered(with_upper_limits=False, with_inf_e_flux=False, sort='filtered')
    sn[i] = sn[i].binned(bin_width=1, discrete_time=True)
    
# i=28 PTF10uhf
# i=22 PTF10aagc
for i in [22]:
    common_bands = []
    for band in sn[i].bands:
        if band in bb.bands:
            common_bands.append(band) #загружаем только фотометрические данные из обычных систем
    if common_bands != sn[i].bands:
        sn[i] = OSCCurve.from_json(os.path.join('./sne', name['Name'][i] + '.json')
                                   , bands=common_bands)
        sn[i] = sn[i].filtered(with_upper_limits=False, with_inf_e_flux=False, sort='filtered')
        sn[i] = sn[i].binned(bin_width=1, discrete_time=True)
        
    if len( sn[i].bands )<=6:
        kern_size = len( sn[i].bands )
    else:
        kern_size = 6
        
     #######################################
    #убираем лишние фильтры для некоторых объектов
    if sn[i].name == 'PTF09cnd':
        kern_size = 5
        
    if sn[i].name == 'SN2011ke':
        kern_size = 5
      
    #######################################
    
    s = 0    
    for u in range(kern_size):
        band_size = len(sn[i].X[(sn[i].X[:,0]==u)])
        if band_size > 2:
            s+= band_size  # берем только полосы, в которых не менее 3 наблюдений
        
    x,y=[],[]
    x=sn[i].X[:s]
    y=sn[i].y[:s]
    err=sn[i].err[:s]
    
    #######################################
    #upper limits
    if sn[i].name == 'iPTF13ajg':
        ind = np.argsort(y)[-1]
        x, y, err = np.delete(x, ind, axis=0), np.delete(y, ind), np.delete(err, ind)
    
    if sn[i].name == 'PTF12dam':
        ind = np.argsort(y)[-2:]
        x, y, err = np.delete(x, ind, axis=0), np.delete(y, ind), np.delete(err, ind)
        kern_size = 4
        
        s = 0    
        for u in range(kern_size):
            band_size = len(sn[i].X[(sn[i].X[:,0]==u)])
            if band_size > 2:
                s+= band_size  # берем только полосы, в которых не менее 3 наблюдений
        
        x,y=[],[]
        x=sn[i].X[:s]
        y=sn[i].y[:s]
        err=sn[i].err[:s]
        
    #######################################
    #cut по дням
    if sn[i].name == 'PTF10aagc':
        rest_frame = (x[:,1] < (55540))
        x = x[rest_frame]
        y = y[rest_frame]
        err = err[rest_frame]
        
    if sn[i].name == 'SN2011ke':
        rest_frame = (x[:,1] < (55780))
        x = x[rest_frame]
        y = y[rest_frame]
        err = err[rest_frame]
        
    if sn[i].name == 'SN2016eay':
        rest_frame = (x[:,1] < (57600))
        x = x[rest_frame]
        y = y[rest_frame]
        err = err[rest_frame]
        
    if sn[i].name == 'SN2015bn':
        rest_frame = (x[:,1] < (57240))
        x = x[rest_frame]
        y = y[rest_frame]
        err = err[rest_frame]
        
    if sn[i].name == 'SN2013dg':
        rest_frame = (x[:,1] < (56550))
        x = x[rest_frame]
        y = y[rest_frame]
        err = err[rest_frame]
        
    if sn[i].name == 'SN2011kg':
        rest_frame = (x[:,1] < (56050))
        x = x[rest_frame]
        y = y[rest_frame]
        err = err[rest_frame]
    #######################################

    #######################################
    
    kern_size = int( x[-1,0] + 1 ) #поскольку фильтровали полосы
    # if kern_size<3:
    #     continue #игнорируем объект, если меньше 3 фильтров
    
    const_matrix = np.eye(kern_size)
    bound_min = np.eye(kern_size) * (-1e1)
    bound_max = np.eye(kern_size) * (1e1)
    for q in range(kern_size):
        const_matrix[q, 0:q] = 1/2
        bound_min[q, 0:q] = -1e1
        bound_max[q, 0:q] = 1e1
    bounds = [bound_min, bound_max]
    
    #size_sc = 3
    #scales = stats.expon.rvs(size=size_sc, loc=1e-5, scale = 10)
    #mk = [MultiStateKernel( [RBF(scales[p],(1e-5,1e3)) for k in range(kern_size)],
    #                       const_matrix, bounds ) for p in range(size_sc)]
    l_edge, r_edge = 0.5, 100
    if sn[i].name == 'SN2013dg':
        r_edge = 30
        
    if sn[i].name == 'LSQ14bdq':
        l_edge, r_edge = 235, 260
    
    # if sn[i].name == 'SN2011ke':
    #     l_edge, r_edge = 1, 5
    
    if sn[i].name == 'SN2007bi':
        l_edge, r_edge = 0.5, 55
    
    list_kern = [RBF(1,(l_edge, r_edge)) for k in range(kern_size-1)]
    list_kern.append(WhiteKernel())
    mk = MultiStateKernel( list_kern, const_matrix, bounds )
    
    
    
    ###########################
    #50 дней до max и 150 после
    max_ind = np.argmax(y)
    day = x[max_ind][1]
    
    if sn[i].name in ['PS1-14bj', 'LSQ14bdq']:
        rest_frame = (x[:,1] > (day - 50)) * (x[:,1] < (day + 250))
    elif sn[i].name == 'iPTF13ajg':
        rest_frame = (x[:,1] > (day - 50)) * (x[:,1] < (day + 100))
    else:
        rest_frame = (x[:,1] > (day - 50)) * (x[:,1] < (day + 150))
    
    x = x[rest_frame]
    y = y[rest_frame]
    err = err[rest_frame]
    err1 = err**2 + (y/10)**2
    ###########################
    #to magn
    magn_arr = sn[i].convert_arrays(x, y, err)
    cur_bands = list(sn[i].bands)[:kern_size]
    y_magn, err_magn = [], []
    for Band in cur_bands:
        y_magn += (- 2.5 * np.log10( magn_arr.odict[Band].y )).tolist()
        err_magn += (2.5 / np.log(10) / magn_arr.odict[Band].y 
                                    * magn_arr.odict[Band].err).tolist()
    # norm = max(y_magn)
    #y_magn, err_magn = np.array(y_magn)/norm, np.array(err_magn)/norm
    y_magn, err_magn = np.array(y_magn), np.array(err_magn)
    err1_magn = err_magn**2 #+ (y_magn/25)**2
    if np.mean(err1_magn) < 0.02:
        err1_magn += 0.026
        
    ###########################

    mask = []
    for j in range(kern_size):
        mask.append( (x[:,0] == j) )
    
###############################################

    gpr = GaussianProcessRegressor(kernel=mk, alpha=err1_magn, 
                                   n_restarts_optimizer=res).fit(x, y_magn) #n_restarts_optimizer=1
    
    b = sn[i].bands[:kern_size]
   
    bins = 3
    if sn[i].name == 'PTF10aagc':
        bins = 1
    n_days = int( ( max(x[:,1]) - min(x[:,1]) ) // bins )
    X = np.linspace(min(x[:,1]), max(x[:,1]), n_days).reshape(-1,1)
    X = np.block([ [u*np.ones_like(X), X] for u in range(kern_size) ] )

    predict, sigma = gpr.predict(X, return_std =True)
    # predict *= norm
    # sigma *= norm
    # y_magn *= norm
    # err_magn *= norm
    #ищем минимальный допустимый индекс
    Y, Sigma, xx = [], [],[]
    min_ind = n_days
    # for u in range(kern_size):
    #     temp =  predict[n_days*u : n_days*u + n_days]
    #     min_ind = min( len(temp[temp > 0]), min_ind)
    # #если есть flux<=0, обрезаем все кривые по этому дню
    for u in range(kern_size):
        temp = predict[n_days*u : n_days*u + n_days][:min_ind]
        Y += temp.tolist()
        Sigma += sigma[n_days*u : n_days*u + n_days][:min_ind].tolist()
        xx += X[n_days*u : n_days*u + n_days][:min_ind].tolist()
        

    #csfont = {'fontname':'Times New Roman'}
    #font = font_manager.FontProperties(family='Times New Roman',
                                   #weight='bold',
                                   #style='normal',
                                   #size=28
     #                              )

    b_legend = [r'\textit{' + B + '}' for B in b]
    fig, ax = plt.subplots(figsize=(10, 7))#figsize=(18, 12),dpi=400
    bb.presets_fig(ax)
    plt.title( name['Name'][i] )
 
    shift = [0, -1, 1]
    shift_lab = ['', ' $-$ 1', ' + 1']
    for u in range(kern_size):
        #shift = (u//2 + u%2)*(-1)**(u+1)
        ax.plot(X[:min_ind,1],
                np.array( Y[min_ind*u : min_ind*u + min_ind]) + shift[u],
                color=bb.cols[b[u]], #label=b[u] + shift_lab[u],
                linewidth=2
                )
        
        ax.errorbar(x[mask[u],1],
                    np.array(y_magn[mask[u]]) + shift[u], err_magn[mask[u]],
                          marker='o', ls='',ms=8, color=bb.cols[b[u]],
                          label=b_legend[u] + shift_lab[u], elinewidth=2
                          )
        
    ax.set_xlabel('MJD')
    ax.set_ylabel('Apparent magnitude')
    ax.set_xlim([X[:min_ind,1][0] - 5, X[:min_ind,1][-1] + 20])
    ax.set_ylim(18.1, 23.8)
    ax.invert_yaxis()
    ax.tick_params(axis='both', direction='in', which='major',  length=8, width=2)
    ax.tick_params(axis='both', direction='in', which='minor',  length=5, width=1.5)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.25))
    ax.grid('on', linestyle='--', alpha=0.7, linewidth=1)

    #get handles and labels
    handles, labels = plt.gca().get_legend_handles_labels()
    #specify order of items in legend
    order = [1,0,2]
    #add legend to plot
    plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],
               loc='upper left', bbox_to_anchor=(0.725, 0.7), labelspacing=3.2) 
   


######################################################
    bins = 3
    if sn[i].name == 'PTF10aagc':
        bins = 1
    n_days = int( ( max(x[:,1]) - min(x[:,1]) ) // bins )
    X = np.linspace(min(x[:,1]), max(x[:,1]), n_days).reshape(-1,1)    

    gpr = []
    predict, sigma = [], []
    for f in range(kern_size):
        kern =  ConstantKernel(1)*RBF(1, (0.5, 100))
        gpr.append( GaussianProcessRegressor(kernel=kern, alpha=err1_magn[mask[f]],
                    n_restarts_optimizer=res).fit(np.reshape(x[mask[f]][:,1], (-1,1)), y_magn[mask[f]])
                   )
   

        pred, sig = gpr[-1].predict(X, return_std =True)
        predict.append(pred)
        sigma.append(sig)
   

    for u in range(kern_size):
        #shift = (u//2 + u%2)*(-1)**(u+1)
        ax.plot(X, predict[u] + shift[u],
                color=bb.cols[b[u]], ls='--', alpha=0.5, linewidth=2)



    
    fig.savefig( fname = 'fig_for_article/' + str(name['Name'][i]) + '_gp1d.pdf',
                    bbox_inches="tight", format='pdf', dpi=1000)

