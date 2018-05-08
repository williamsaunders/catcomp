import os
import astropy.table as tb
import numpy as np
os.chdir('/data3/garyb/tno/matcher')
coadd = []
final = []
y4a1 = []
for i in os.listdir('.'):
      os.chdir(i)
      for j in os.listdir('.'):
          if j.endswith('_y4a1.fits'):
              y4a1.append(j.split('_')[0])
          if j.endswith('_coadd.fits'):
              coadd.append(j.split('_')[0])
          if j.endswith('_final.fits') and not j.startswith('zone'):
              final.append(j.split('_')[0])
      os.chdir('..')

a = np.setdiff1d(y4a1, final)

#This will say which tiles are *not* done. 
#Then:
os.chdir('/home/wsaund/CatalogComparison')
zones = tb.Table.read('zones.fits') #attached

zz = zones[np.in1d(zones['TILENAME'], a)]

np.unique(zz['ZONE']) #is the list of zones that are *not done*!

#Finally

n = np.setdiff1d(range(216), np.unique(zz['ZONE']))

f = open('zones_doneFINAL.txt', 'w')
for nn in n:
    f.write('%d \n'%nn)
f.close()    


