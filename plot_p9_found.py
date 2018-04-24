'''
THIS IS THE DISCOVERY STEP OF THE PLANET NINE FLEET. 
DISCOVERY CONDITIONS:
1. At least 4 detections on different dates 
2. Determine probability of detecting based on exposure detection function 
'''

import numpy as np
from numpy import matrix
import matplotlib
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.time import Time
import sys
import timeit
import argparse
import random 

# First start with the bright limit.  >= 4 observations of an object on different dates = detection!

p9 = fits.getdata('P9simulation_results/P9results.fits')
print len(p9)

mcolors = []
for name, hex in matplotlib.colors.cnames.iteritems():
    mcolors.append(name)

colors = ['r','g','b','k','y','orange','c','m']
icount = 0


for i, ob in enumerate(p9):
    print i, '/', len(p9)
    if i == 0:
        co = random.choice(mcolors)
        plt.scatter(ob['ra'], ob['dec'])
    elif ob['ob_num'] == p9[i-1]['ob_num']:
            continue
    else:
        plt.scatter(ob['ra'], ob['dec'])

plt.title('Scatterplot of Detected P9s in Bright Limit')
plt.xlabel('ra [deg]')
plt.ylabel('dec [deg]')
plt.axis('equal')
plt.savefig('P9detectionsALL.png', dpi=100)
plt.show()

