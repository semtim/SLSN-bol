import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from multistate_kernel import MultiStateKernel
from snad.load.curves import OSCCurve
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel
import os
import bb
import matplotlib.ticker as ticker

##########
#paths for saving approx data

path_fig = os.path.abspath('approx_magn/figures/')
path_data = os.path.abspath('approx_magn/data/')

##########


temp = os.path.abspath("second_cut.csv")
name = pd.read_csv(temp, sep=",")
name = pd.DataFrame(name)

res=3
sn=[[] for i in range(len(name))]

# for i in range(29):
#     sn.append( OSCCurve.from_name(name['Name'][i], 
#                 down_args={"baseurl": "https://sne.space/sne/"}) )


for i in range(len(name)):
    sn[i] = OSCCurve.from_json(os.path.join('./sne', name['Name'][i] + '.json'))
    sn[i] = sn[i].filtered(with_upper_limits=False, with_inf_e_flux=False, sort='filtered')
    sn[i] = sn[i].binned(bin_width=1, discrete_time=True)
    
####################################################################


for i in range(len(name)):
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
    
    const_matrix = np.eye(kern_size)
    bound_min = np.eye(kern_size) * (-1e1)
    bound_max = np.eye(kern_size) * (1e1)
    for q in range(kern_size):
        const_matrix[q, 0:q] = 1/2
        bound_min[q, 0:q] = -1e1
        bound_max[q, 0:q] = 1e1
    bounds = [bound_min, bound_max]
    


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


    y_magn, err_magn = np.array(y_magn), np.array(err_magn)
    err1_magn = err_magn**2 
    if np.mean(err1_magn) < 0.02:
        err1_magn += 0.026
        
    ###########################

    mask = []
    for j in range(kern_size):
        mask.append( (x[:,0] == j) )
    
    ###########################################

    gpr = GaussianProcessRegressor(kernel=mk, alpha=err1_magn, 
                                   n_restarts_optimizer=res).fit(x, y_magn)
    
    b = sn[i].bands[:kern_size]
   
    bins = 3
    if sn[i].name == 'PTF10aagc':
        bins = 1
    n_days = int( ( max(x[:,1]) - min(x[:,1]) ) // bins )
    X = np.linspace(min(x[:,1]), max(x[:,1]), n_days).reshape(-1,1)
    X = np.block([ [u*np.ones_like(X), X] for u in range(kern_size) ] )

    predict, sigma = gpr.predict(X, return_std =True)
    #ищем минимальный допустимый индекс
    Y, Sigma, xx = [], [],[]
    min_ind = n_days
    for u in range(kern_size):
        temp = predict[n_days*u : n_days*u + n_days][:min_ind]
        Y += temp.tolist()
        Sigma += sigma[n_days*u : n_days*u + n_days][:min_ind].tolist()
        xx += X[n_days*u : n_days*u + n_days][:min_ind].tolist()
        
    Y, Sigma = np.array(Y), np.array(Sigma)



    
    #making plots
    font = {'family' : 'Times New Roman',
        'size'   : 22}

    plt.rc('font', **font)
    plt.rcParams['axes.linewidth'] = 1.2
    
    fig, ax = plt.subplots(figsize=(10, 7))
    plt.title( name['Name'][i] )
 
    for u in range(kern_size):
        ax.plot(X[:min_ind,1], Y[min_ind*u : min_ind*u + min_ind],
                color=bb.cols[b[u]], linewidth=2)
        ax.fill_between(X[:min_ind,1],
                        Y[min_ind*u : min_ind*u + min_ind] + Sigma[min_ind*u : min_ind*u + min_ind],
                        Y[min_ind*u : min_ind*u + min_ind] - Sigma[min_ind*u : min_ind*u + min_ind],
                        alpha=0.3, color=bb.cols[b[u]])
        
        plt.errorbar(x[mask[u],1], y_magn[mask[u]], err_magn[mask[u]],
                          marker='o', ls='',ms=8, color=bb.cols[b[u]],
                          label=b[u], elinewidth=2)
        
        
    ax.set_xlabel('MJD')
    ax.set_ylabel('Apparent magnitude')
    ax.set_xlim([X[:min_ind,1][0] - 3, X[:min_ind,1][-1] + 3])
    ax.invert_yaxis()

    ax.tick_params(axis='both', direction='in', which='major',  length=8, width=2)
    ax.tick_params(axis='both', direction='in', which='minor',  length=5, width=1.5)
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
    ax.grid('on', linestyle='--', alpha=0.7, linewidth=1)




    plt.legend(loc='best', columnspacing=0.7, labelspacing=0.3,
               handletextpad=0.35,) #framealpha=0.0,
    
    
    fig.savefig( fname = path_fig + '/' + str(name['Name'][i]) + '.pdf',
                    bbox_inches="tight", format='pdf')

    plt.close()



################################################
#saving data
    
    columns = ['mjd']
    for r in b:
        columns.append(r)
        columns.append('err_' + r)
    
    approx_data = pd.DataFrame(columns=columns)
    approx_data["mjd"] = X[:min_ind, 1]
    
    for p in range(kern_size):
        approx_data[b[p]] = Y[min_ind*p : min_ind*p + min_ind]
        approx_data['err_' + b[p]] = sigma[min_ind*p : min_ind*p + min_ind]
    
    approx_data.to_csv(path_data + '/' + str(name['Name'][i]) + '.csv',
                         index=False)
