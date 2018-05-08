'''
This program creates the histogram of m50 values for each exposure.
'''

import numpy as np
import glob
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size' : 18})
import astropy.io.fits as fits
import sys

data = fits.getdata('zone_efficiencies/all-coadd_detection_results.fits')
m50s = {'g':[], 'r':[], 'i':[], 'z':[], 'Y':[]}
bands =  ['g', 'r', 'i', 'z', 'Y']
for band in bands:
    m50band = data[data['band']==band]['m50']
    m50s[band].append(m50band)
sys.stdout.flush()
bins = np.arange(19,26,.05)
plt.figure(figsize=(15,8))
num = 0
for band, co in zip(bands, ['b', 'g', 'y', 'orange', 'r']):
    m50 = np.hstack(m50s[band])
    plt.hist(m50, bins=bins, align='mid', histtype='step', normed=True, label=band, color=co, linewidth=2)
    plt.axvline(np.median(m50), color=co, linestyle='dashed', label='%s median = %.2f'%(band, np.median(m50)), linewidth=2)
    num += len(m50)
plt.legend()
plt.xlim(19,26)
plt.ylim(0,2.6)
plt.title('Detections Histogram (%d exposures)'%num)
plt.xlabel('m50 Magnitude')
plt.ylabel('Frequency (normalized)')
plt.grid()
plt.tight_layout()
plt.savefig('ALLzones_histogram.png', dpi=100)
plt.show()
