'''
This program does the simulation of Planet Nine orbits and detection on DES.
It has a bunch of command line arguments to make it easier to qsub.  
It requires a number of starting sphere points for P9. 
This takes a long time to run ~days, even when divided into 10 chunks, so I ran chunks 
separately and then combined them. 
'''

import numpy as np
from numpy import matrix
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy.io import fits
import astropy.units as u
import astropy.constants as c
import ephem
from mpl_toolkits.mplot3d import Axes3D
from itertools import product, combinations
import math, random
from scipy import spatial
from astropy.time import Time
import sys
import timeit
import matplotlib.path as mpath
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('sphere', type=str, help='file of sphere points to run as Planet Nine "fleet"')
parser.add_argument('--qsub', type=int, help='optional argument to parallellize and qsub the fleet.  Takes an integer = # of separate chunks to parallelize.  Requires running recombine_p9.py after completion.')
parser.add_argument('--chunk', type=int, help='Number chunk.  Number between 0 and qsub-1')
parser.add_argument('--plot', action='store_true', help='make plots of the p9 orbits. SLOW.')
args = parser.parse_args()

def eq_to_ecl_cart(x,y,z):
    i_e = 23.43928
    a = matrix([[1,0,0],[0,np.cos(i_e*u.degree), np.sin(i_e*u.degree)], \
                [0, -np.sin(i_e*u.degree), np.cos(i_e*u.degree)]])
    b = matrix([[x],[y],[z]])
    c = a*b
    return c[0,0], c[1,0], c[2,0] #x, y, z
    
def eq_to_ecl_cart_v(vx,vy,vz):
    # converts equatorial cartesian coodinates to ecliptic cartesian coordinates
    i_e = 23.43928
    a = matrix([[1,0,0],[0,np.cos(i_e*u.degree), -np.sin(i_e*u.degree)], \
                [0, np.sin(i_e*u.degree), np.cos(i_e*u.degree)]])
    b = matrix([[vx],[vy],[vz]])
    c = a*b
    return c[0,0], c[1,0], c[2,0] #x, y, z   


def ecl_to_cart(lat, lon, r):
    # converts ecliptic spherical coorindates to ecliptic cartesian coodinates
    x = r*np.cos(lat*u.degree)*np.cos(lon*u.degree)
    y = r*np.cos(lat*u.degree)*np.sin(lon*u.degree)
    z = r*np.sin(lat*u.degree)
    xx, yy, zz = eq_to_ecl_cart(x.value, y.value, z.value)
    return xx, yy, zz

def v_xyz(lat, lon, r, PA):
    # converts ecliptic spherical coordinates into ecliptic cartesian velocity vectors 
    v = np.sqrt(c.G*c.M_sun/((r*u.au).to(u.m)))
    v = v.to(u.au/u.year)
    v = v.value
    v_vec = v*matrix([[np.cos(lon*u.degree), np.sin(lon*u.degree), 0], \
                      [-np.sin(lon*u.degree), np.cos(lon*u.degree), 0], \
                      [0, 0, 1]])* \
              matrix([[np.cos(lat*u.degree), 0, -np.sin(lat*u.degree)], \
                      [0, 1, 0], \
                      [np.sin(lat*u.degree), 0, np.cos(lat*u.degree)]])* \
              matrix([[0], [-np.sin(PA*u.degree)], [np.cos(PA*u.degree)]])
    vx, vy, vz = eq_to_ecl_cart_v(v_vec[0,0], v_vec[1,0], v_vec[2,0])
    return np.array([vx, vy, vz])
    
def distance(ra1, ra2, dec1, dec2):
    # calculate distance on spherical geometry 
    dist = np.sqrt((ra2 - ra1)**2*np.cos(dec1*u.degree) + (dec2 - dec1)**2)
    return dist.value

def footprint(ra, dec):
    # determine if point is within the footprint 
    p = np.loadtxt('round17-poly.txt', comments='#', dtype=float)
    path = mpath.Path(p)
    inside = path.contains_point((ra, dec))
    return inside
    
def jd_to_date(jd):
    time = Time(jd, format='jd')
    return time.iso

def mjd_to_date(mjd):
    time = Time(mjd, format='mjd')
    return time.iso

def date_to_jd(date):
    time = Time(date, format='iso')
    return time.jd

def date_to_mjd(date):
    time = Time(date, format='iso')
    return time.mjd

def astro_topo(body):
    # we want astronomical, topological coordinates, but PyEphem doesn't provide them.
    # this calculates them given a body 
    body_ra = body.a_ra + (body.ra - body.g_ra)
    body_dec = body.a_dec + (body.dec - body.g_dec)
    body_c = SkyCoord(ra=body_ra*u.rad, dec=body_dec*u.rad, frame='icrs')
    return body_c

def coord_to_hex(coord):
    # don't think you need this....
    ra_h_decimal = coord.ra.hour
    ra_h = np.floor(ra_h_decimal)
    ra_m_decimal = (ra_h_decimal - ra_h)*60.
    ra_m = np.floor(ra_m_decimal)
    ra_s = (ra_m_decimal - ra_m)*60.
    ra_string = '%d:%d:%f' %(ra_h, ra_m, ra_s)

    dec_d_decimal = coord.dec.degree
    dec_d = np.floor(dec_d_decimal)
    dec_m_decimal = (dec_d_decimal - dec_d)*60.
    dec_m = np.floor(dec_m_decimal)
    dec_s = (dec_m_decimal - dec_m)*60.
    dec_string = '%d:%d:%f' %(dec_d, dec_m, dec_s)

    return ra_string, dec_string
              
# import corners file
loc = 'y4a1-noY.ccdcorners.fits'
corners = fits.getdata(loc)
print 'setup complete'
sys.stdout.flush()

# import the appropriate sphere points
lon, lat = np.loadtxt(args.sphere, delimiter=',', unpack=True)

print '%.0f sphere points ready'%len(lon)
sys.stdout.flush()

# DEFINE ORBITAL PROPERTIES FOR EACH SET OF COORDINATES
r = 400
ob_num + (args.chunk+1)*100000

ob_num_col = []
expnum_col = []
CCD_col = []
date_col = []
ra_col = []
dec_col = []
num_col = []

if args.qsub:
    lat_chunk = np.array_split(lat, args.qsub)
    lon_chunk = np.array_split(lon, args.qsub)
    lat = lat_chunk[args.chunk]
    lon = lon_chunk[args.chunk]

for lat, lon in zip(lat, lon):
    for PA in np.linspace(0,300,6):
        start = timeit.default_timer()
        print 'starting position: (RA, Dec, PA)', lon, lat, PA
        x, y, z = ecl_to_cart(lat, lon, r)
        v_vec = v_xyz(lat, lon, r, PA)
        r_vec = np.array([x,y,z])
        
        # a
        a = np.linalg.norm(r_vec)
        
        # e 
        e = 0
        
        # M
        M = 0
        epoch = '2015/01/01 12:00:00'
        
        # i
        L_vec = np.cross(r_vec, v_vec)
        i = np.arccos(L_vec[2]/np.linalg.norm(L_vec))*u.rad.to(u.degree)
        
        # Omega
        Omega = np.arctan2(L_vec[1], L_vec[0])*u.rad.to(u.degree) + 90.
        
        # omega
        Omega_hat = np.array([np.cos(Omega*u.degree), np.sin(Omega*u.degree), 0])
        r_hat = r_vec/np.linalg.norm(r_vec)
        if z >= 0.:
            omega = np.arccos(np.dot(Omega_hat, r_hat))*u.rad.to(u.degree)
        else: 
            omega = -np.arccos(np.dot(Omega_hat, r_hat))*u.rad.to(u.degree)
    
    
        # USE ORBIT/CCD LOCATING CODE TO IDENTIFY IF IT APPEARED ON DES
        
        # initialize object and add orbital elements
        o = ephem.EllipticalBody()
        
        o._e = e
        o._a = a 
        o._Om = Omega
        o._inc = i
        o._om = omega
        o._M = M 
        o._epoch_M = epoch
        o._epoch = epoch
        
        # observer object for CTIO that is same for all observations
        ob = ephem.Observer()
        ob.lat = -30.169117
        ob.lon = 289.194100
        ob.elevation = 2389
        ob.pressure = 0
        print 'observer set'
        sys.stdout.flush()
        
        # DETERMINE THE CCDs THAT THIS OBJECT LANDED ON 
        
        # Now divide up all the data into months
        remainder = np.remainder(len(corners['mjd_mid']), np.floor(len(corners['mjd_mid'])/40))
        breaks = np.zeros(41)
        for i, b in enumerate(breaks):
            breaks[i] = i*np.floor(len(corners['mjd_mid'])/40.)
        breaks[-1] = breaks[-1] + remainder
        breaks = np.array(breaks, dtype=int)
        breaks = breaks[1:-1]
        month_breaks = np.split(corners['mjd_mid'], breaks)
    
        # For each observation, determine which month we care about
        overlaps = []
        if args.plot:
            fig, ax = plt.subplots(figsize=(16,10))
        for m in month_breaks:
            # pick out month from corners.fits
            month = corners[corners['mjd_mid'] >= np.min(m)]
            month = month[month['mjd_mid'] <= np.max(m)]
        
            # find the middle date and get location of object at that time
            mid = (np.max(m) - np.min(m))/2. + np.min(m)
            mid_day = mjd_to_date(mid)
            ob.date = mid_day
            o.compute(ob)
            o_c = astro_topo(o)
            ra_mid = o_c.ra.value
            dec_mid = o_c.dec.value
        
            # to determine radius of tree search, see how far the object moves in the month
            month_start = np.min(m)
            ob.date = mjd_to_date(month_start)
            o.compute(ob)
            o_c_start = astro_topo(o)
        
            month_end = np.max(m)
            ob.date = mjd_to_date(month_end)
            o.compute(ob)
            o_c_end = astro_topo(o)
        
            sep = o_c_start.separation(o_c_end)
        
            # build tree to search for near neighbors around mid position
            treedata = zip(month['ra'][:,4]*np.cos(month['dec'][:,4]), month['dec'][:,4])
            tree = spatial.cKDTree(treedata)
            #near = tree.query_ball_point((ra_mid*np.cos(dec_mid), dec_mid),r=sep.degree/2.+ .17)

            # the above near 'near' uses the "correct" radius search for how far a P9
            # might travel, but it didn't find all of the CCDs for some reason.
            # the below 'near' uses a bit of a larger radius to be safe

            near = tree.query_ball_point((ra_mid*np.cos(dec_mid), dec_mid),r=sep.degree/2.+ .3)
            if near == []:
                continue
            else:
                month_near = month[near]
                del month
        
            # now brute force final search through candidates
            for i, ccd in enumerate(month_near):
                mjd = ccd['mjd_mid']
                date = mjd_to_date(mjd)
                ob.date = date
                o.compute(date)
                o_c = astro_topo(o)
                ra = o_c.ra.value
                dec = o_c.dec.value
        
                if ((np.min(ccd['ra']) <= ra) and (np.max(ccd['ra']) >= ra) and \
                (np.min(ccd['dec']) <= dec) and (np.max(ccd['dec']) >= dec)):
                    # condition for observation 
                    if args.plot:
                        # this makes the plots of orbits, CCDs, and detections 
                        ax.add_patch(matplotlib.patches.Rectangle((np.min(ccd['ra']), \
                                     np.min(ccd['dec'])), 
                                     np.max(ccd['ra']) - np.min(ccd['ra']), \
                                     np.max(ccd['dec']) - np.min(ccd['dec']), alpha=0.05))
                        plt.scatter(ra, dec, marker='o', c='r', s=50, edgecolors='r')
                    sys.stdout.flush()
                    overlaps.append((ccd, ra, dec))
        print '# overlapping CCDs: ', len(overlaps)
        if args.plot:
            plt.axis('equal')
            plt.title('At Starting Epoch: RA=%.2f Dec=%.2f PA=%.0f' %(lon, lat, PA))
            plt.xlabel('ra [deg]')
            plt.ylabel('dec [deg]')
        # plot the whole orbit across the area we care about
        times = []
        if args.plot:
            print 'plotting whole orbit'
            sys.stdout.flush()
            for row in overlaps:
                times.append(row['mjd_mid'])
            if not times:
                continue
            start_time = np.min(times)
            end_time = np.max(times)
            all_ra = []
            all_dec = []
            for time in np.linspace(start_time-((end_time-start_time)*.1), end_time+((end_time-start_time)*.1), 1000):
                ob.date = mjd_to_date(time)
                o.compute(ob)
                o_c = astro_topo(o)
                all_ra.append(o_c.ra.value)
                all_dec.append(o_c.dec.value)
                plt.scatter(o_c.ra.value, o_c.dec.value, marker='o', c='k', edgecolors='k', s=3)
            plt.plot(all_ra, all_dec, c='k', linewidth=2)
            plt.xlim(np.min(all_ra) - .5, np.max(all_ra) + .5)    
            plt.ylim(np.min(all_dec) - .5, np.max(all_dec) + .5)    
            plt.axis('equal')
            plt.show()
            fig.savefig('P9/RA=%.2f,Dec=%.2f,PA=%.0f.png' %(lon, lat, PA), dpi=500, bbox_inches='tight')
            plt.close()
        
        for obs in overlaps:
            ob_num_col.append(ob_num)
            expnum_col.append(obs[0]['expnum'])
            CCD_col.append(obs[0]['detpos'])
            ra_col.append(obs[1])
            dec_col.append(obs[2])
            date_col.append(obs[0]['mjd_mid'])
            num_col.append(len(overlaps))
        end = timeit.default_timer()
        ob_num += 1
        print 'time %.1f seconds' %(end-start)
           
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

if args.qsub:
    name = 'p9_results-chunk' + str(args.chunk) + '.fits'
else:
    name = 'p9_resuls.fits'
t = fits.BinTableHDU.from_columns([c1,c2,c3,c4,c5,c6,c7])
t.writeto(name, clobber=True)
