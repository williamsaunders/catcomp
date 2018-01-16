# -*- coding: utf-8 -*-
"""
Created on Fri Nov 03 18:17:44 2017

@author: William

NOTE: THIS FILE RUNS BETTER ON FOLIO
"""

import numpy as np
import astropy.io.fits as fits
import os



#files = os.listdir('.')
files = ['DES0135-0458_final.fits', 'DES0136-0624_final.fits']

ra_col = []
dec_col = []
expnum_col = []
mag_auto_g_col = []
mag_auto_r_col = []
mag_auto_i_col = []
mag_auto_z_col = []
mag_auto_y_col = []
band_col = []
match_id_col = []
coadd_object_id_col = []
spread_model_col = []
for file1 in files:
    print file1
    if 'final' in file1:
        fitspath = file1
        data = fits.getdata(fitspath)
        ra_col.append(data['ra'])
        dec_col.append(data['dec'])
        expnum_col.append(data['expnum'])
        mag_auto_g_col.append(data['mag_auto_g'])
        mag_auto_r_col.append(data['mag_auto_r'])
        mag_auto_i_col.append(data['mag_auto_i'])
        mag_auto_z_col.append(data['mag_auto_z'])
        mag_auto_y_col.append(data['mag_auto_y'])
        band_col.append(data['band'])
        match_id_col.append(data['match_id'])
        coadd_object_id_col.append(data['coadd_object_id'])
        spread_model_col.append(data['spread_model'])
        del data

ra_col = np.array(ra_col)
dec_col = np.array(dec_col)
expnum_col = np.array(expnum_col)
mag_auto_g_col = np.array(mag_auto_g_col)
mag_auto_r_col = np.array(mag_auto_r_col)
mag_auto_i_col = np.array(mag_auto_i_col)
mag_auto_z_col = np.array(mag_auto_z_col)
mag_auto_y_col = np.array(mag_auto_y_col)
band_col = np.array(band_col)
match_id_col = np.array(match_id_col)
coadd_object_id_col = np.array(coadd_object_id_col)
spread_model_col = np.array(spread_model_col)

ra_col = np.hstack(ra_col)
dec_col = np.hstack(dec_col)
expnum_col = np.hstack(expnum_col)
mag_auto_g_col = np.hstack(mag_auto_g_col)
mag_auto_r_col = np.hstack(mag_auto_r_col)
mag_auto_i_col = np.hstack(mag_auto_i_col)
mag_auto_z_col = np.hstack(mag_auto_z_col)
mag_auto_y_col = np.hstack(mag_auto_y_col)
band_col = np.hstack(band_col)
match_id_col = np.hstack(match_id_col)
coadd_object_id_col = np.hstack(coadd_object_id_col)
spread_model_col = np.hstack(spread_model_col)

c1 = fits.Column(name='ra', array=ra_col, format='F')
c2 = fits.Column(name='dec', array=dec_col, format='F')
c3 = fits.Column(name='expnum', array=expnum_col, format='D')
c4 = fits.Column(name='mag_auto_g', array=mag_auto_g_col, format='F')
c5 = fits.Column(name='mag_auto_r', array=mag_auto_r_col, format='F')
c6 = fits.Column(name='mag_auto_i', array=mag_auto_i_col, format='F')
c7 = fits.Column(name='mag_auto_z', array=mag_auto_z_col, format='F')
c8 = fits.Column(name='mag_auto_Y', array=mag_auto_y_col, format='F')
c9 = fits.Column(name='band', array=band_col, format='A')
c10 = fits.Column(name='match_id', array=match_id_col, format='D')
c11 = fits.Column(name='coadd_object_id', array=coadd_object_id_col, format='D')
c12 = fits.Column(name='spread_model', array=spread_model_col, format='F')

t = fits.BinTableHDU.from_columns([c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12])
a = t.data
#t.writeto('combined_final.fits')

