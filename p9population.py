import numpy as np
from numpy import matrix
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy.io import fits
import astropy.units as u
import ephem
from mpl_toolkits.mplot3d import Axes3D
from itertools import product, combinations
import math, random

def fibonacci_sphere(samples, randomize=True):
    rnd = 1.
    if randomize:
        rnd = random.random() * samples

    points = []
    offset = 2./samples
    increment = math.pi * (3. - math.sqrt(5.));

    for i in range(samples):
        y = ((i * offset) - 1) + (offset / 2);
        r = math.sqrt(1 - pow(y,2))

        phi = ((i + rnd) % samples) * increment

        x = math.cos(phi) * r
        z = math.sin(phi) * r

        points.append([x,y,z])

    return points

def ecl_to_cart(lat, lon, r):
    x = r*np.cos(lat)*np.cos(lon)
    y = r*np.cos(lat)*np.sin(lon)
    z = r*np.sin(lat)
    return x, y, z

# first create a set of (ecliptic) lat / lon positions evenly spaced
'''
lat = []
lon = []
num_lon_vals = []
max_num_lon = 5.
lat_vals = np.arange(-90,91,1)
for lat in lat_vals:
    if lat == 90 or lat == -90:
        num_lon = 1.
    else:
        num_lon = np.floor(max_num_lon*np.cos(lat*u.degree).value)
    num_lon_vals.append(num_lon)

coords = {'lat' : [], 'lon' : []}

for lat, num_lon in zip(lat_vals, num_lon_vals):
    lons = np.linspace(0, 360., num_lon)
    for lon in lons:
        coords['lat'].append(lat)
        coords['lon'].append(lon)

# plot sphere

for lat, lon in zip(coords['lat'], coords['lon']):
    x, y, z = ecl_to_cart(lat, lon, 1.)
    ax.scatter(x, y, z, color='k', s=50
'''
coords = fibonacci_sphere(samples=3000)
coords = np.array(coords)

fig = plt.figure(figsize=(13,13))
ax = fig.add_subplot(111, projection='3d')
#ax = fig.gca(projection='3d')
ax.set_aspect("equal")

u = np.linspace(0, 2 * np.pi, 1000)
v = np.linspace(0, np.pi, 1000)
x = np.outer(np.cos(u), np.sin(v))
y = np.outer(np.sin(u), np.sin(v))
z = np.outer(np.ones(np.size(u)), np.cos(v))
ax.plot_surface(x, y, z, color='b', linewidth=0, alpha=.3)
'''
u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
x = np.cos(u)*np.sin(v)
y = np.sin(u)*np.sin(v)
z = np.cos(v)
ax.plot_grid(x, y, z, color="r")
'''
ax.scatter(coords[:,0], coords[:,1], coords[:,2], c='k', s=10, edgecolor='k')
plt.tight_layout()
plt.savefig('sphere_points.png', dpi=200)
plt.show()

