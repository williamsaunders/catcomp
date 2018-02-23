import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['font.size']=14
import astropy.io.fits as fits
import timeit
import sys
import os
from scipy import spatial
from scipy import optimize
from scipy import special
import argparse

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
        return -result

parser = argparse.ArgumentParser()
parser.add_argument('zonelist', type=str, help='file with list of zones to run')
args = parser.parse_args()

zonelist = np.loadtxt(args.zonelist, dtype=str)
for zone in zonelist:
    fitspath = zone + '/' + zone + '-combined_final.fits'
    data = fits.getdata(fitspath) # everything!
    coadd_table = data[data['expnum']==999999]
    coadd_stars = coadd_table[np.abs(coadd_table['spread_model_i']) <= 0.003] # only stars from main table
    print 'stars identified'
    sys.stdout.flush()

    #%% PERFORM NEAREST NEIGHBOR SEARCH 
    start = timeit.default_timer()
    f = open('%s-coadd_detection_results.csv'%zone, 'w')
    f.write('exposure, band, m50, k, c, coadds found, coadds, missed minusLogP  \n')

    # build decision tree
    treedata = zip(coadd_stars['ra'], coadd_stars['dec'])
    tree = spatial.cKDTree(treedata, leafsize=16)

    # identify objects
    corners = fits.getdata('y4a1.ccdcorners.fits') # THIS NEEDS TO BE WORKED ON 
#    corners = corners[corners['type']=='segmap']
    expnums = np.unique(data['expnum'])

    for expnum in expnums[0:-1]:
        print '----->', expnum
        sys.stdout.flush()  
        band = data[data['expnum']==expnum]['band'][0]
        coadds_exp = []
        coadds_exp_found = []
        coadds_exp_missed = []
        corners_exp = corners[corners['expnum']==expnum]
        data_exp_single = data[data['expnum']==expnum]
        ccd_list = corners_exp['detpos']
        for ccd in ccd_list:
            corners_ccd = corners_exp[corners_exp['detpos']==ccd]
            if len(corners_ccd) > 1:
                raise TypeError('more than 1 ccd identified')
            if len(corners_ccd) == 0:
                continue
            ra_center = corners_ccd[0]['ra'][4]
            dec_center = corners_ccd[0]['dec'][4]
            if not ra_center:
                continue
            else: 
                # identify near neighbors
#                ra_center = ra_center[0]
#                dec_center = dec_center[0]
                coadd_near_neighbors = 0
                coadd_near_neighbors = tree.query_ball_point([ra_center, dec_center], r=1)
                if not coadd_near_neighbors:
#                    print 'NO NEAR NEIGHBORS'
                    continue
                coadd_ball = coadd_stars[coadd_near_neighbors]            
                near_neighborsRA = coadd_ball['ra']
                near_neighborsDEC = coadd_ball['dec']
                # identify the objects on the ccd
                keep1 = np.logical_and(near_neighborsRA >= np.min(corners_ccd['ra']), \
                                       near_neighborsRA <= np.max(corners_ccd['ra']))
                keep2 = np.logical_and(near_neighborsDEC >= np.min(corners_ccd['dec']), \
                                       near_neighborsDEC <= np.max(corners_ccd['dec']))
                keep_coadd = np.logical_and(keep1, keep2)
                coadd_ccd = coadd_ball[keep_coadd]
                if coadd_ccd.size == 0:
#                    print 'NO NEAREST NEIGHBORS'
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
        results_collector = [] # place to store log of likelihood data 
        optimized = optimize.minimize(minusLogP, (22, 5, .98), method='Nelder-Mead', \
                                  args=(coadds_exp_found, coadds_exp_missed), tol=1e-3)
        if optimized.success:
            opt_params = optimized.x
        else:
            print optimized.message

        f.write('%d, %s, %.2f, %.3f, %.4f, %d, %d %.2f \n'%(expnum,band,optimized.x[0],optimized.x[1],optimized.x[2],len(coadds_exp_found), len(coadds_exp_missed), results_collector[-1]))
        end = timeit.default_timer()
        print 'time: %.1f seconds' %(end - start)
    f.close()
    break
