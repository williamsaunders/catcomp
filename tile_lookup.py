### THIS CODE IS A TEMPORARY WAY FOR THE NEAREST NEIGHBOR SEARCH TO LOOK UP WHICH TILES IT NEEDS TO CONSIDER FOR A GIVEN OBJECT.  ###


import numpy as np
import astropy.io.fits as fits
import os

all_filenames = os.listdir('.')
filenames = [a for a in all_filenames if 'final' in a]

f = open('tile_lookup.csv', 'w')
f.write('filename, median_ra, median_dec \n')

for filename in filenames[:-1]:
    print filename
    t = fits.getdata(filename)
    tmedRA = np.median(t['ra'])
    tmedDEC = np.median(t['dec'])
    f.write('%s, %.2f, %.2f \n' %(filename, tmedRA, tmedDEC))
    del t

f.close()
    
