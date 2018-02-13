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
    m0, k, c = params
    logit = c/(1+np.exp(k*(m-m0)))
    return logit
def minusLogP(params, mdet, mnon):
    '''
    Takes negative log of above functions to use in optimizing parameters by finding minimum
    params =     
    Returns negative log of pdf
    '''
    if params[2] > 1.:
        results_collector.append(results_collector[-1] + 1e3)
        return results_collector[-1] + 1e3
    else:
        pdet = detprob(mdet,params)
        pnon = detprob(mnon,params)    
        result = np.sum(np.log(pdet))
        result += np.sum(np.log(1-pnon))
        results_collector.append(-result)
        return -result

plt.rcParams['font.size']=14
zonepath = 'zone29'
fitspath = zonepath + '/' + zonepath + '-combined_final.fits'
data = fits.getdata(fitspath) # everything!
coadd_data = data[data['expnum']==999999] # TEMP
coadds = fits.getdata('zone29/zone29-combined_coadd.fits') # TEMPORARY FOR ZONE29
coadd_stars = coadds[np.abs(coadds['spread_model_i']) <= 0.003] #temporary
coadd_bool = np.in1d(coadd_data['coadd_object_id'], coadd_stars['coadd_object_id']) # TEMP
coadd_table = coadd_data[coadd_bool] # TEMP
#temp##coadd_table = data[data['expnum']==999999]
#tempt##coadd_stars = coadd_table[np.abs(coadd_table['spread_model_i']) <= 0.003] # only stars from main table
coadd_stars = coadd_table
print 'stars identified'
sys.stdout.flush()

#%% PERFORM NEAREST NEIGHBOR SEARCH 
start = timeit.default_timer()
f = open('%s-coadd_detection_results.csv'%zonepath, 'w')
f.write('exposure, band, m50, k, c, coadds found, coadds missed \n')

# build decision tree
treedata = zip(coadd_stars['ra'], coadd_stars['dec'])
tree = spatial.cKDTree(treedata, leafsize=16)

# identify objects
corners = fits.getdata('corners.fits')
corners = corners[corners['type']=='segmap']
expnums = np.unique(data['expnum'])

for expnum in expnums[0:-1]:
    print '----->', expnum
    sys.stdout.flush()  
    band = data[data['expnum']==expnum]['band'][0]
    coadds_exp = []
    coadds_exp_found = []
    coadds_exp_missed = []
    corners_exp = corners[corners['exposure']==expnum]
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
                print 'NO NEAR NEIGHBORS'
                continue
            coadd_ball = coadd_stars[coadd_near_neighbors]            
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
                print 'NO NEAREST NEIGHBORS'
                continue
            coadds_exp.append(coadd_ccd)
    coadds_exp = np.array(coadds_exp)
    if len(coadds_exp) == 0.:
        print "EXPOSURE EMPTY"
        continue
    coadds_exp = np.hstack(coadds_exp)
    for tile in np.unique(coadds_exp['tile']):
        coadds_exp_tile = coadds_exp[coadds_exp['tile']==tile]
        data_exp_single_tile = data_exp_single[data_exp_single['tile']==tile]
        single_bool = np.in1d(coadds_exp_tile['match_id'], data_exp_single_tile['match_id'])
        coadds_found = coadds_exp_tile[single_bool]
        coadds_missed = coadds_exp_tile[np.logical_not(single_bool)]
        coadds_exp_found.append(coadds_found['mag_auto_%s'%band])
        coadds_exp_missed.append(coadds_missed['mag_auto_%s'%band])

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
    results_collector = []
    optimized = optimize.minimize(minusLogP, (22, 5, .98), method='Nelder-Mead', \
                                  args=(coadds_exp_found, coadds_exp_missed), tol=1e-3)
    if optimized.success:
        opt_params = optimized.x
    else:
        print optimized.message

    f.write('%d, %s, %.2f, %.3f, %.4f, %d, %d \n'%(expnum,band,optimized.x[0],optimized.x[1],optimized.x[2],len(coadds_exp_found), len(coadds_exp_missed)))
    print '%d, %s, %.2f, %.3f, %.4f, %d, %d'%(expnum,band,optimized.x[0],optimized.x[1],optimized.x[2],len(coadds_exp_found), len(coadds_exp_missed))
    
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
    plt.savefig(zonepath + '/plots/%d.png'%expnum)
    plt.close()
    
    end = timeit.default_timer()
    print 'time: %.1f seconds' %(end - start)
f.close()
