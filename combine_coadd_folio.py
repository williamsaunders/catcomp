# -*- coding: utf-8 -*-
"""
Created on Fri Nov 03 18:17:44 2017

@author: William

NOTE: THIS FILE RUNS BETTER ON FOLIO
"""

import numpy as np
import astropy.io.fits as fits
import os


files = os.listdir('.')

ra_col = []
dec_col = []
coadd_object_id_col = []
spread_model_col = []
for file1 in files:
    print file1
    if 'coaddinfo' in file1:
        fitspath = file1
        data = fits.getdata(fitspath)
        ra_col.append(data['ra'])
        dec_col.append(data['dec'])
        coadd_object_id_col.append(data['coadd_object_id'])
        spread_model_col.append(data['spread_model_i'])
        del data

ra_col = np.array(ra_col)
dec_col = np.array(dec_col)
coadd_object_id_col = np.array(coadd_object_id_col)
spread_model_col = np.array(spread_model_col)

ra_col = np.hstack(ra_col)
dec_col = np.hstack(dec_col)
coadd_object_id_col = np.hstack(coadd_object_id_col)
spread_model_col = np.hstack(spread_model_col)

c1 = fits.Column(name='ra', array=ra_col, format='F')
c2 = fits.Column(name='dec', array=dec_col, format='F')
c3 = fits.Column(name='coadd_object_id', array=coadd_object_id_col, format='D')
c4 = fits.Column(name='spread_model_i', array=spread_model_col, format='F')

t = fits.BinTableHDU.from_columns([c1, c2, c3, c4])
t.writeto('combined_coadds.fits')

