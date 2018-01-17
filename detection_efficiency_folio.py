import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import timeit
import sys
import os
from scipy import spatial
from scipy import optimize
from scipy import special

def detprob(m, params):
    '''
    logit function
    params = (m0, k, c)
    '''
    m0, k, c, a = params
#    logit = c/(1+np.exp(k*(m-m0)))
    logit = a + (c-a)/(1+np.exp(k*(m-m0)))
    return logit
def minusLogP(params, mdet, mnon):
    '''
    Takes negative log of above functions to use in optimizing parameters by finding minimum
    params =     
    Returns negative log of pdf
    '''
    pdet = detprob(mdet,params)
    pnon = detprob(mnon,params)    
    result = np.sum(np.log(pdet))
    result += np.sum(np.log(1-pnon))
    return -result

plt.rcParams['font.size']=14
fitspath = 'combined_final.fits'
coaddpath = 'combined_coadds.fits'

data = fits.getdata(fitspath)
data[np.logical_or(np.abs(data['spread_model'])  <= 0.003, data['expnum']==999999)]
coadddata = fits.getdata(coaddpath)
coadddata = coadddata[np.abs(coadddata['spread_model_i'])<=0.003]


results_table = []

print 'setup done'
sys.stdout.flush()
#%% PERFORM NEAREST NEIGHBOR SEARCH 
start = timeit.default_timer()
#f = open('coadd_detection_results.csv', 'w')
#f.write('exposure, band, m50, k, c, coadds found, coadds missed \n')
# build decision tree
treedata = zip(coadddata['ra'], coadddata['dec'])
tree = spatial.cKDTree(treedata, leafsize=16)

# lookup in tree
corners = fits.getdata('corners.fits')
corners = corners[corners['type']=='segmap']
expnums = np.unique(data['expnum'])
for expnum in expnums[:10]:
    print '----->', expnum
    sys.stdout.flush()  
    band = data[data['expnum']==expnum]['band'][0]
    coadds_exp_found = []
    coadds_exp_missed = []
    corners_exp = corners[corners['exposure']==expnum]
    data_exp_coadd = data[data['expnum']==999999]
    data_exp_single = data[data['expnum']==expnum]
    for ccd in range(1,63):
        corners_ccd = corners_exp[corners_exp['ccd']==ccd]
        if len(corners_ccd) > 1:
            raise TypeError('more than 1 ccd identified')
        ra_center = 0.5*(corners_ccd['racmax'] + corners_ccd['racmin'])
        dec_center = 0.5*(corners_ccd['deccmax'] + corners_ccd['deccmin'])
        if not ra_center:
            continue
        else: 
            # identify near neighbors
            ra_center = ra_center[0]
            dec_center = dec_center[0]
            coadd_near_neighbors = 0
            coadd_near_neighbors = tree.query_ball_point([ra_center, dec_center], r=1)
            if not coadd_near_neighbors:
                continue
            coadd_ball = coadddata[coadd_near_neighbors]            
            near_neighborsRA = coadd_ball['ra']
            near_neighborsDEC = coadd_ball['dec']
            # identify the objects on the ccd
            keep1 = np.logical_and(near_neighborsRA >= corners_ccd['racmin'], \
                                   near_neighborsRA <= corners_ccd['racmax'])
            keep2 = np.logical_and(near_neighborsDEC >= corners_ccd['deccmin'], \
                                   near_neighborsDEC <= corners_ccd['deccmax'])
            keep_coadd = np.logical_and(keep1, keep2)
            coadd_ccd = coadd_ball[keep_coadd]
            if coadd_ccd.size == 0:
                continue

            coadd_comp_bool = np.in1d(data_exp_coadd['coadd_object_id'], coadd_ccd['coadd_object_id'])
            coadds_total = data_exp_coadd[coadd_comp_bool]
            matches_bool = np.in1d(coadds_total['match_id'], data_exp_single['match_id'])
            coadds_found = coadds_total[matches_bool]
            coadds_missed = coadds_total[np.logical_not(matches_bool)]
            
#            print "---> identified all objects on ccd"
            coadds_exp_found.append(coadds_found['mag_auto_%s'%band])                          
            coadds_exp_missed.append(coadds_missed['mag_auto_%s'%band])

#            print "---> identified all matches on ccd"
#            sys.stdout.flush()
    coadds_exp_found = np.array(coadds_exp_found)
    coadds_exp_found = np.hstack(coadds_exp_found)   
    coadds_exp_found = coadds_exp_found[coadds_exp_found >= 18.]                      
    coadds_exp_found = coadds_exp_found[coadds_exp_found <= 28.]                      
    coadds_exp_missed = np.array(coadds_exp_missed)
    coadds_exp_missed = np.hstack(coadds_exp_missed)                            
    coadds_exp_missed = coadds_exp_missed[coadds_exp_missed >= 18.]
    coadds_exp_missed = coadds_exp_missed[coadds_exp_missed <= 28.]
    print 'nearest neighbor lookup complete'
    sys.stdout.flush()
    
    # optimize parameters for logit fit
    print 'optimizing...'
    sys.stdout.flush()
    optimized = optimize.minimize(minusLogP, (24, 2, .95, 0), method='Nelder-Mead', \
                                  args=(coadds_exp_found, coadds_exp_missed), tol=1e-4)
    if optimized.success:
        opt_params = optimized.x
        print opt_params
    else:
        print optimized.message
    '''    
    # second round of optimization with better range
    print '          ...2'
    sys.stdout.flush()
#    coadds_exp_found = coadds_exp_found[coadds_exp_found >= (optimized.x[0]-2)]                      
    coadds_exp_found = coadds_exp_found[coadds_exp_found <= (optimized.x[0]+2)]                      
#    coadds_exp_missed = coadds_exp_missed[coadds_exp_missed >= (optimized.x[0]-2)]
    coadds_exp_missed = coadds_exp_missed[coadds_exp_missed <= (optimized.x[0]+2)]

    optimized = optimize.minimize(minusLogP, (24, 2, .95), method='Nelder-Mead', \
                                  args=(coadds_exp_found, coadds_exp_missed), tol=1e-4)
    if optimized.success:
        opt_params = optimized.x
        print opt_params
    else:
        print optimized.message
    '''
    
#    f.write('%d, %s, %.2f, %.3f, %.4f, %d, %d \n'%(expnum,band,optimized.x[0],optimized.x[1],optimized.x[2],\
#                                            len(coadds_exp_found), len(coadds_exp_missed)))
    
    plt.figure(figsize=(13,9))
    bins = np.linspace(18, 28, 20)
    bins_center = ((bins + np.roll(bins, 1))/2)[1:]
    # first bin the data    
    hist_found = np.array(np.histogram(coadds_exp_found, bins=bins)[0], dtype=float)
    hist_missed = np.array(np.histogram(coadds_exp_missed, bins=bins)[0], dtype=float)
    frac = hist_found / (hist_found + hist_missed)
    # plot binned results    
    ax1 = plt.subplot(211)
    plt.scatter(bins_center, frac, s=25, color='b')
    # plot logit pdf
    logit = detprob(np.linspace(18, 28, 500), opt_params)
    ax1.plot(np.linspace(18, 28, 500), logit, color='k', linewidth=2)
    # plot totals for context     
    ax2 = plt.subplot(212, sharex=ax1)
    plt.bar(bins[:-1], (hist_found + hist_missed),color='r', alpha=.6, width=.4)
    ax2.set_yscale('log')
    ax1.set_title('Detection Results for Exposure %d (%s band)' %(expnum, band))
    plt.xlabel('Magnitude')
    ax1.set_ylabel('Detection Frequency')
    ax2.set_ylabel('Number of Coadd Objects')
    ax1.legend()
    ax1.grid()
    ax2.grid()
    plt.xlim(18, 28)
    plt.show()
    
    plt.savefig('DOUBLEOPT-test/%d.png'%expnum)
    plt.close()
    end = timeit.default_timer()
    print 'time: %.1f seconds' %(end - start)
#f.close()