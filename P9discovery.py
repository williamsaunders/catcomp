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

# identify the total number of objects with at least 1 detection on a DES CCD
num_detect = []
total_num, num_detect = Totals(p9, num_detect=num_detect)
print total_num, len(num_detect)
print num_detect

# perform calculation of objects with minimum unique (different date) detections >= (3,4,5,6,7)
# and different mag thresholds for each

discovery_dict = {}

discovered = []
count_unique = []
fleet_size = 29388.

for m in np.append([99.], np.arange(23,25,.5)):
    discovery_dict['mag%.1f'%m] = []
    for i in range(3,9):
        print 'det = ', i, 'mag = ', m
        det, count = Discovery(p9, i, magthresh=m, discovered=discovered, 
                               count_unique=count_unique)
        discovery_dict['mag%.1f'%m].append(len(det3))

# save dictionary to file --> plotting in different script
pickle.dump(discovery_dict, open('discovery.fits', 'w'))

