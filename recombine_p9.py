# THIS PROGRAM COMBINES OUTPUTS OF p9population.py INTO ONE FITS TABLE. IT ALSO MAKES BACKUP

import numpy as np
import astropy.io.fits as fits
from astropy.table import Table, hstack
from glob import glob
'''
p9chunks = glob('p9results/9p_results-chunk*.fits')
table_chunks = []
for chunk in p9chunks:
    t = Table.read(chunk, format='fits')
    table_chunks.append(t)
new = hstack(table_chunks)
new.write('P9results.fits')
'''
p9chunks = glob('P9simulation_results/p9_results-chunk*.fits')
print p9chunks

ob_num_col = []
expnum_col = []
CCD_col = []
ra_col = []
dec_col = []
date_col = []
num_col = []
for chunk in p9chunks:
    print chunk
    t = fits.getdata(chunk)
    for i in range(len(t)):
        ob_num_col.append(t[i]['ob_num'])
        expnum_col.append(t[i]['expnum'])
        CCD_col.append(t[i]['CCD'])
        ra_col.append(t[i]['ra'])
        dec_col.append(t[i]['dec'])
        date_col.append(t[i]['date'])
        num_col.append(t[i]['num'])

ob_num_col = np.array(ob_num_col)
expnum_col = np.array(expnum_col)
CCD_col = np.array(CCD_col)
ra_col = np.array(ra_col)
dec_col = np.array(dec_col)
date_col = np.array(date_col)
num_col = np.array(num_col)

c1 = fits.Column(name='ob_num', array=ob_num_col, format='D')
c2 = fits.Column(name='expnum', array=expnum_col, format='D')
c3 = fits.Column(name='ccd', array=CCD_col, format='A3')
c4 = fits.Column(name='date', array=date_col, format='F')
c5 = fits.Column(name='ra', array=ra_col, format='F')
c6 = fits.Column(name='dec', array=dec_col, format='F')
c7 = fits.Column(name='num', array=num_col, format='D')

t = fits.BinTableHDU.from_columns([c1,c2,c3,c4,c5,c6,c7])
t.writeto('P9simulation_results/P9results.fits', clobber=True)
