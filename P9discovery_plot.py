'''
This program reads discovery results from pickle dump files and makes 2 plots: 
discovery results plot and a histogram of the number of detections in a season.
'''

import numpy as np
import matplotlib
matplotlib.use('TkAgg')
matplotlib.rcParams.update({'font.size' : 14})
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import pickle
import sys
import timeit
from discovery import Totals, Discovery

# import P9 simulation results
p9 = fits.getdata('P9simulation_results/P9results.fits')
fleet_size = 29388.

# identify the total number of objects with at least 1 detection on a DES CCD
num_detect = []
total_num, num_detect = Totals(p9, num_detect=num_detect)
print '------->', total_num

# load results from discovery
discovery = pickle.load(open('discovery.fits'))
discovery_freq = {} 
discovery_omega = {}
for k in discovery.keys():
    discovery_freq[k] = np.array(discovery[k])/fleet_size
    discovery_omega[k] = (np.array(discovery[k])/fleet_size)*41253*(1./8)
pickle.dump(discovery_freq, open('discovery_freq.fits', 'w'))
pickle.dump(discovery_omega, open('discovery_omega.fits', 'w'))
    
# have to do another discovery step to get the count of the number of detections 
discovered = []
count_unique = []
det3, count3, missing = Discovery(p9, 3, magthresh=99., discovered=discovered, count_unique=count_unique)

# collect and plot statistics
fig, ax1 = plt.subplots(figsize=(12,8))
ax1.plot(range(3,8), np.array(discovery['mag99.0'])/fleet_size, label='Bright limit', \
         color='k', linewidth=3, zorder=6)
ax1.plot(range(3,8), np.array(discovery['mag22.5'])/fleet_size, label='Magnitude 22.5', \
         color='b', linewidth=3, zorder=5)
ax1.plot(range(3,8), np.array(discovery['mag23.0'])/fleet_size, label='Magnitude 23.0', \
         color='g', linewidth=3, zorder=4)
ax1.plot(range(3,8), np.array(discovery['mag23.5'])/fleet_size, label='Magnitude 23.5', \
         color='orange', linewidth=3, zorder=3)
ax1.plot(range(3,8), np.array(discovery['mag24.0'])/fleet_size, label='Magnitude 24.0', \
         color='r', linewidth=3, zorder=2)
ax1.plot(range(3,8), np.array(discovery['mag24.5'])/fleet_size, label='Magnitude 24.5', \
         color='m', linewidth=3, zorder=1)

ax1.set_yticks(np.arange(0, 0.6, .05))
ax1.set_xticks([3,4,5,6,7])
ax1.tick_params('y', colors='b')
ax1.xaxis.grid(True)
ax1.yaxis.grid(True, linestyle='dashed')
ax1.set_xlabel('Detection Number Threashold')
ax1.set_ylabel('Frequency Detected on Footprint', color='b')
plt.title('Full Detection Results')
ax2 = ax1.twinx()
ax2.set_ylabel(r'Effective Area ($\Omega_A$) [sq. deg.]', color='r')
Omega = np.arange(0,0.55, .05)*41253*(1./8)
ax2.set_yticks(Omega)
ax2.tick_params('y', colors='r')
plt.grid(linestyle='dashed')
ax1.set_ylim([0,.5])
ax1.legend()
plt.tight_layout()
plt.savefig('discovery_results.png', dpi=100)
plt.close('all')

# histogram plot of the number of detections for everything 
plt.figure(figsize=(12,6))
plt.hist(num_detect, bins=np.arange(1,50,1), align='mid', color='b', label='With 6 Hour Re-Detections', alpha=.5, histtype='step', lw=3)
plt.hist(count3, bins=np.arange(1,50,1), align='mid', color='r', label='Without 6 Hour Re-Detections', alpha=.5, histtype='step', lw=3)
plt.title('Histogram of Planet Nine Detections')
plt.xlabel('Number of Detections')
plt.ylabel('Number')
plt.grid()
plt.legend()
plt.savefig('detection_histogram.png', dpi=100)
plt.close('all')


