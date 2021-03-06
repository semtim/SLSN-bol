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

##########
#paths for saving approx data
path_fig = os.path.abspath('approx_magn_results\\figures\\')
path_data = os.path.abspath('approx_magn_results\\data\\')
##########
temp = os.path.abspath("second_cut.csv")
name = pd.read_csv(temp, sep=",")
name = pd.DataFrame(name)

res=3
sn=[[] for i in range(len(name))]

# for i in range(29):
#     sn.append( OSCCurve.from_name(name['Name'][i], 
#                 down_args={"baseurl": "https://sne.space/sne/"}) )

#bands = ["u","g","r","i","z"]

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
    colour = ['purple', 'green', 'red', 'brown', 'grey', 'black']               
    col = dict( zip( b, colour[:kern_size] ) )
    
   
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
        

    
    fig, ax = plt.subplots()#figsize=(18, 12),dpi=400
    plt.title( name['Name'][i] )
 
    for u in range(kern_size):
        ax.plot(X[:min_ind,1], Y[min_ind*u : min_ind*u + min_ind],
                color=bb.cols[b[u]])
        
        plt.errorbar(x[mask[u],1], y_magn[mask[u]], err_magn[mask[u]],
                          marker='x', ls='',ms=5, color=bb.cols[b[u]])
        
        #plt.errorbar(x[mask[u],1], y[mask[u]], np.sqrt(err1[mask[u]]),
                      # marker='x', ls='',ms=3, color=bb.cols[b[u]])
    ax.set_xlabel('mjd')
    ax.set_ylabel('Apparent magnitude')
    ax.set_xlim([X[:min_ind,1][0] - 3, X[:min_ind,1][-1] + 3])
    ax.invert_yaxis()

    
    t = [ 'error '+ b[l] for l in range(kern_size) ]
    lab = np.array(b).tolist()+t
    plt.legend(lab,framealpha=0.0,loc='upper right',ncol=2)

    fig.savefig( fname = path_fig + '\\' + str(name['Name'][i]),
                    bbox_inches="tight")


################################################
#saving data and figures
    
    columns = ['mjd']
    for r in b:
        columns.append(r)
        columns.append('err_' + r)
    
    approx_data = pd.DataFrame(columns=columns)
    approx_data["mjd"] = X[:min_ind, 1]
    
    # xx, Y, Sigma = np.array(xx), np.array(Y), np.array(Sigma)
    # data = sn[i].convert_arrays(xx, Y, Sigma) #перевод данных в потоки
    # for band in b:
    #     approx_data[band] = - 2.5 * np.log10( data.odict[band].y )
    #     approx_data['err_' + band] = 2.5 / np.log(10) / data.odict[band].y * data.odict[band].err
    for p in range(kern_size):
        approx_data[b[p]] = Y[min_ind*p : min_ind*p + min_ind]
        approx_data['err_' + b[p]] = sigma[min_ind*p : min_ind*p + min_ind]
    
    approx_data.to_csv(path_data + '\\' + str(name['Name'][i]) + '.csv',
                         index=False)


    # x, y, err = np.array(x), np.array(y), np.array(err)
    # data_raw = sn[i].convert_arrays(x, y, err)
    # fig, ax = plt.subplots() #figsize=(18, 12),dpi=400
    # for b in approx_data.columns[1::2]:
    #     plt.plot(approx_data["mjd"], approx_data[b], c=bb.cols[b])
    #     y_raw =  - 2.5 * np.log10( data_raw.odict[b].y )
    #     x_raw = data_raw.odict[b].x
    #     plt.scatter(x_raw, y_raw, c=bb.cols[b])
    # plt.xlabel('mjd')
    # plt.ylabel('Apparent magnitude')
    # plt.title(sn[i].name)
    # plt.legend(approx_data.columns[1::2].tolist())
    # plt.xlim(approx_data["mjd"][0] - 10, 
    #          approx_data["mjd"][len(approx_data["mjd"])-1] + 10)
    # ax.invert_yaxis()
    
    # if not np.isnan(np.max(approx_data.to_numpy()[:,1::2])):
    #     plt.ylim(np.isnan(np.max(approx_data.to_numpy()[:,1::2])) + 0.5, 
    #          np.min(approx_data.to_numpy()[:,1::2]) - 0.5)
    
    # fig.savefig( fname = path_fig + '\\magn\\' + str(name['Name'][i]),
                # bbox_inches="tight")
#########################################################

