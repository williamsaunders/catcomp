import numpy as np
from numpy import matrix
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy.io import fits
import astropy.units as u
import astropy.constants as c
from mpl_toolkits.mplot3d import Axes3D
from itertools import product, combinations
import math, random
from scipy import spatial
from astropy.time import Time
import sys
import timeit
import matplotlib.path as mpath
import argparse

'''
THIS PROGRAM IS TO BE RUN ONCE TO CREATE 10 FILES EACH OF WHICH HAS 1/10 OF THE TOTAL 
SPHERE POINTS GENERATED.  SPHERE POINTS ARE USED TO CREATE THE FLEET OF PLANET NINES
SPHERE POINTS ARE ONLY THOSE WITHIN THE DES FOOTPRINT
'''

def footprint(ra, dec):
    p = np.loadtxt('round17-poly.txt', comments='#', dtype=float)
    path = mpath.Path(p)
    inside = path.contains_point((ra, dec))
    return inside

lat_lon = []
N = 500000.
Ncount = 0.
while Ncount < N:
    print Ncount
    if Ncount < N:
        r = 1.
        a = (4*np.pi*r**2)/N
        d = np.sqrt(a)
        M_theta = np.round(np.pi/d)
        d_theta = np.pi/M_theta
        d_phi = a/d_theta
        for m in range(int(M_theta)):
            print m, '/', int(M_theta)
            theta = np.pi*(m + .5)/M_theta
            M_phi = np.round(2*np.pi*np.sin(theta)/d_phi)
            if Ncount >= N:
                break
            for n in range(int(M_phi)):
#                print n, '/', int(M_phi)
                if Ncount >= N:
                    break
                phi = 2*np.pi*n/(M_phi)
                lat = np.pi/2 - theta
                lon = phi
                lat = lat*u.rad.to(u.degree)
                lon = lon*u.rad.to(u.degree)
                # Now test if it's on the DES footprint
                if footprint(lon, lat):
                    lat_lon.append([lat, lon])
                Ncount += 1.
print Ncount
print len(lat_lon)

lat_lon = np.array(lat_lon)

print '%.0f sphere points ready'%len(lat_lon)
sys.stdout.flush()

np.random.shuffle(lat_lon)
s0, s1, s2, s3, s4, s5, s6, s7, s8, s9 = np.array_split(lat_lon, 10)

for i, s in enumerate([s0, s1, s2, s3, s4, s5, s6, s7, s8, s9]):
    print i
    fig = plt.figure(figsize=(13,13))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_aspect("equal")
    uu = np.linspace(0, 2 * np.pi, 1000)
    vv = np.linspace(0, np.pi, 1000)
    x = np.outer(np.cos(uu), np.sin(vv))
    y = np.outer(np.sin(uu), np.sin(vv))
    z = np.outer(np.ones(np.size(uu)), np.cos(vv))
    r = 1.
    xx = r*np.cos(s[:,0]*u.degree)*np.cos(s[:,1]*u.degree)
    yy = r*np.cos(s[:,0]*u.degree)*np.sin(s[:,1]*u.degree)
    zz = r*np.sin(s[:,0]*u.degree)
    ax.plot_surface(x, y, z, color='b', linewidth=0, alpha=.3)
    ax.scatter(xx, yy, zz, c='k', s=10, edgecolor='k')
    plt.tight_layout()
    plt.savefig('sphere_points/sphere_points%d.png'%i, dpi=200)
    #plt.show()
    f = open('sphere_points/sphere_points' + str(i) + '.csv', 'w')
    for dec, ra in s:
        f.write('%f, %f \n'%(ra, dec))
    f.close
