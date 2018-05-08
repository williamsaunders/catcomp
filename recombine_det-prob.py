'''
This program recombines all the exposure m50 tables into one table, choosing the exposure
that has more objects when duplicates arise. 
Write the exposure final statistics into one table to look up during discovery step. 
'''

import astropy.io.fits as fits
import numpy as np
from glob import glob

expnum_col = np.array([])
band_col = np.array([])
m50_col = np.array([])
c_col = np.array([])
k_col = np.array([])
found_col = np.array([])
missed_col = np.array([])
mLogP_col = np.array([])

files = glob('zone_efficiencies/*.fits')
for file1 in files:
    print file1[18:25]
    data = fits.getdata(file1)
    data1 = data[data['band']!='Y'] # remove Y for printing statistics, but keep in for now
    med = np.median(data1['k'])
    sd = np.std(data1['k'])
    data2 = data[data['k']<=(med + 2*sd)]
    print len(data), '--->', len(data1), '--->', len(data2) # see how much the list thins
    del data, data1
    for line in data2:
        if line['expnum'] in expnum_col:
            f = np.where(line['expnum']==expnum_col)[0]
            if len(f) != 1:
                sys.exit('EXPNUM DUPLICATES')
                # something went wrong
            g = expnum_col[f]
            if (found_col[f] + missed_col[f]) < (line['found'] + line['missed']):
                # see if the one already in the list or the new one has more coadd objects
                c_col[f] = line['c']
                k_col[f] = line['k']
                m50_col[f] = line['m50']
                found_col[f] = line['found']
                missed_col[f] = line['missed']
                mLogP_col[f] = line['mLogP']
            continue
        expnum_col = np.append(expnum_col, line['expnum'])
        band_col = np.append(band_col, line['band'])
        m50_col = np.append(m50_col, line['m50'])
        c_col = np.append(c_col, line['c'])
        k_col = np.append(k_col, line['k'])
        found_col = np.append(found_col, line['found'])
        missed_col = np.append(missed_col, line['missed'])
        mLogP_col = np.append(mLogP_col, line['mLogP'])

c1 = fits.Column(name='expnum', array=expnum_col, format='D')
c2 = fits.Column(name='band', array=band_col, format='A')
c3 = fits.Column(name='m50', array=m50_col, format='F')
c4 = fits.Column(name='k', array=k_col, format='F')
c5 = fits.Column(name='c', array=c_col, format='F')
c6 = fits.Column(name='found', array=found_col, format='D')
c7 = fits.Column(name='missed', array=missed_col, format='D')
c8 = fits.Column(name='mLogP', array=mLogP_col, format='F')

t = fits.BinTableHDU.from_columns([c1, c2, c3, c4, c5, c6, c7, c8])
t.writeto('zone_efficiencies/all-coadd_detection_results.fits', clobber=True)
