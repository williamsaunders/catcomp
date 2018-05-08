'''
This program calculates the detection efficiency of point-sources in each exposure. 
'''

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
    Takes logit parameters and list of detected and non-detected magnitudes. 
    Returns negative log of pdf
    '''
    if params[2] > 1.:
        # this ensures that c cannot continue to rise above 1 by referencing 
        # the value of the previous trial in the optimizer
        results_collector.append(results_collector[-1] + (params[2] - 1.)*1e5)
        return results_collector[-1] + 1e3
    elif params[2] <= 0.:
        return results_collector[-1] + 1e3
    else:
        pdet = detprob(mdet,params)
        pnon = detprob(mnon,params)    
        result = np.sum(np.log(pdet))
        result += np.sum(np.log(1-pnon))
        return -result

parser = argparse.ArgumentParser()
parser.add_argument('zonepath', type=str, help='path to the directory with all the zones')
parser.add_argument('zonelist', type=str, help='file with list of zones to run')
parser.add_argument('--plots', action='store_true', help='bool, make sigmoid plots')
args = parser.parse_args()

zonelist = np.loadtxt(args.zonelist, dtype=str)
for zone in zonelist:
    z = str(int(zone)).zfill(3)
    print '------------>', 'zone', z
    fitspath = args.zonepath + '/' + 'zone' + z + '/' + 'zone' + z + '-combined_final.fits'
    data = fits.getdata(fitspath) # load everything
    coadd_table = data[data['expnum']==999999]
    coadd_stars = coadd_table[np.abs(coadd_table['spread_model_i']) <= 0.003] # only stars from main table
    print 'stars identified'
    sys.stdout.flush()

    start = timeit.default_timer()
    
    # create variables to store eventually as fits columns to write to a fits table
    exposure_col = []
    band_col = []
    m50_col = []
    k_col = []
    c_col = []
    coadds_found_col = []
    coadds_missed_col = []
    minusLogP_col = []
    
    # build k-d tree

    # the corners file has ra going from -180 to +180,
    # so we have to adjust the coadd stars going into the decision tree to match 
    great180 = coadd_stars['ra'] > 180.
    coadd_stars_ra = coadd_stars['ra']
    coadd_stars_ra[great180] = coadd_stars_ra[great180] - 360.
    treedata = zip(coadd_stars_ra, coadd_stars['dec'])
    tree = spatial.cKDTree(treedata, leafsize=16)

    # identify objects
    corners = fits.getdata('y4a1.ccdcorners.fits') 
    expnums = np.unique(data['expnum'])

    for expnum in expnums[:-1]:
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
                # identify near neighbors with k-d tree search 
                coadd_near_neighbors = 0
                coadd_near_neighbors = tree.query_ball_point([ra_center, dec_center], r=.5)
                if not coadd_near_neighbors:
                    continue
                coadd_ball = coadd_stars[coadd_near_neighbors]            
                near_neighborsRA = coadd_ball['ra']
                near_neighborsDEC = coadd_ball['dec']

                # brute-force the final search of near neighbors to find observations on CCD
                keep1 = np.logical_and(near_neighborsRA >= np.min(corners_ccd['ra']), \
                                       near_neighborsRA <= np.max(corners_ccd['ra']))
                keep2 = np.logical_and(near_neighborsDEC >= np.min(corners_ccd['dec']), \
                                       near_neighborsDEC <= np.max(corners_ccd['dec']))
                keep_coadd = np.logical_and(keep1, keep2)
                coadd_ccd = coadd_ball[keep_coadd]
                if coadd_ccd.size == 0:
                    continue
                else:
                    coadds_exp.append(coadd_ccd)
        coadds_exp = np.array(coadds_exp)
        if len(coadds_exp) == 0.:
            continue
        coadds_exp = np.hstack(coadds_exp)
        
        # have to iterate through each tile separately because MATCH_IDs repeat 
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
        sys.stdout.flush()
    
        # optimize parameters for logit fit
        results_collector = [0] # place to store log of likelihood data 
        optimized = optimize.minimize(minusLogP, (22, 5, .98), method='Powell', \
                                  args=(coadds_exp_found, coadds_exp_missed), tol=1e-3)
        if optimized.success:
            opt_params = optimized.x
        else:
            print optimized.message
            print len(coadds_exp_found), len(coadds_exp_missed)

        exposure_col.append(expnum)
        band_col.append(band)
        m50_col.append(optimized.x[0])
        k_col.append(optimized.x[1])
        c_col.append(optimized.x[2])
        coadds_found_col.append(len(coadds_exp_found))
        coadds_missed_col.append(len(coadds_exp_missed))
        minusLogP_col.append(results_collector[-1])

        # make plots if argument is triggered
        if args.plots:
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
            plt.savefig('zonetest/%d.png'%expnum)
            plt.close()

 
    
        end = timeit.default_timer()
        print 'time: %.1f seconds' %(end - start)
    exposure_col = np.array(exposure_col)
    band_col = np.array(band_col)
    m50_col = np.array(m50_col)
    k_col = np.array(k_col)
    c_col = np.array(c_col)
    coadd_found_col = np.array(coadds_found_col)
    coadds_missed_col = np.array(coadds_missed_col)
    minusLogP_col = np.array(minusLogP_col)
    
    c1 = fits.Column(name='expnum', array=exposure_col, format='D')
    c2 = fits.Column(name='band', array=band_col, format='A')
    c3 = fits.Column(name='m50', array=m50_col, format='F')
    c4 = fits.Column(name='k', array=k_col, format='F')
    c5 = fits.Column(name='c', array=c_col, format='F')
    c6 = fits.Column(name='found', array=coadds_found_col, format='D')
    c7 = fits.Column(name='missed', array=coadds_missed_col, format='D')
    c8 = fits.Column(name='mLogP', array=minusLogP_col, format='F')

    t = fits.BinTableHDU.from_columns([c1, c2, c3, c4, c5, c6, c7, c8])
    t.writeto('zone_efficiencies/zone' + z + '-coadd_detection_results.fits', clobber=True)

