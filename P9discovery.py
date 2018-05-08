'''
This program loads functions from discovery.py and runs the discovery step of the P9 search.
Results are saved as pickle dump files to be read and plotted by P9discovery_plot.py
'''

import numpy as np
import matplotlib
import astropy.io.fits as fits
import pickle
import sys
import timeit
from discovery import Totals, Discovery

# import P9 simulation results
p9 = fits.getdata('P9simulation_results/P9results.fits')

# perform calculation of objects with minimum unique (different date) detections >= (3,4,5,6,7)
# and different mag thresholds for each

discovery_dict = {} # collect discovery statistics
missing_exp_dict = {} # used to collect exposures that don't have point-source stats


fleet_size = 29388.
start = timeit.default_timer()
for m in np.append([99.], np.arange(22.5,25,.5)):
    discovery_dict['mag%.1f'%m] = []
    for i in range(3,8):
        discovered = []  # can probably get rid of this
        count_unique = [] # this too
        print 'det = ', i, '|',  'mag = ', m
        # do the discovery step 
        det, count, missing_exp = Discovery(p9, i, magthresh=m, discovered=discovered, 
                                            count_unique=count_unique)
        discovery_dict['mag%.1f'%m].append(len(det))
        missing_exp_dict['mag%.1f-det%d'%(m,i)] = np.array(missing_exp)
        end = timeit.default_timer()
        print 'time: %.1f s' %(end - start)

# save dictionary to file --> plotting in different script
pickle.dump(discovery_dict, open('discovery.fits', 'w'))
pickle.dump(missing_exp_dict, open('missing_exposures.fits', 'w'))



