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
import matplotlib.path as mpath

def eq_to_ecl_cart(x,y,z):
    i_e = 23.43928
    a = matrix([[1,0,0],[0,np.cos(i_e*u.degree), np.sin(i_e*u.degree)], \
                [0, -np.sin(i_e*u.degree), np.cos(i_e*u.degree)]])
    b = matrix([[x],[y],[z]])
    c = a*b
    return c[0,0], c[1,0], c[2,0] #x, y, z
    
def eq_to_ecl_cart_v(vx,vy,vz):
    i_e = 23.43928
    a = matrix([[1,0,0],[0,np.cos(i_e*u.degree), -np.sin(i_e*u.degree)], \
                [0, np.sin(i_e*u.degree), np.cos(i_e*u.degree)]])
    b = matrix([[vx],[vy],[vz]])
    c = a*b
    return c[0,0], c[1,0], c[2,0] #x, y, z   


def ecl_to_cart(lat, lon, r):
    x = r*np.cos(lat*u.degree)*np.cos(lon*u.degree)
    y = r*np.cos(lat*u.degree)*np.sin(lon*u.degree)
    z = r*np.sin(lat*u.degree)
    xx, yy, zz = eq_to_ecl_cart(x.value, y.value, z.value)
#    print (x.value,y.value,z.value)
#    print (xx,yy,zz)
    return xx, yy, zz

def v_xyz(lat, lon, r, PA):
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
#    print (v_vec[0,0], v_vec[1,0], v_vec[2,0])
#    print (vx, vy, vz)
    return np.array([vx, vy, vz])
    
def distance(ra1, ra2, dec1, dec2):
    dist = np.sqrt((ra2 - ra1)**2*np.cos(dec1*u.degree) + (dec2 - dec1)**2)
    return dist.value

def footprint(ra, dec):
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
    body_ra = body.a_ra + (body.ra - body.g_ra)
    body_dec = body.a_dec + (body.dec - body.g_dec)
    body_c = SkyCoord(ra=body_ra*u.rad, dec=body_dec*u.rad, frame='icrs')
    return body_c

def coord_to_hex(coord):
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
loc = 'y4a1.ccdcorners.fits'
corners = fits.getdata(loc)
print 'setup complete'
sys.stdout.flush()

# create set of x,y,z, coordinates of evenly spaced points on a sphere

coords = ([],[],[])
lat_lon = ([],[])
#N = 2e6
N = 4990.
Ncount = 0.
while Ncount < N:
    if Ncount < N:
        r = 1.
        a = (4*np.pi*r**2)/N
        d = np.sqrt(a)
        M_theta = np.round(np.pi/d)
        d_theta = np.pi/M_theta
        d_phi = a/d_theta
        for m in range(int(M_theta)):
            theta = np.pi*(m + .5)/M_theta
            M_phi = np.round(2*np.pi*np.sin(theta)/d_phi)
            if Ncount >= N:
                break
            for n in range(int(M_phi)):
                if Ncount >= N:
                    break
                phi = 2*np.pi*n/(M_phi)
                lat = np.pi/2 - theta
                lon = phi
                lat = lat*u.rad.to(u.degree)
                lon = lon*u.rad.to(u.degree)
                # Now test if it's on the DES footprint
                if footprint(lon, lat):
                    lat_lon[0].append(lat)
                    lat_lon[1].append(lon)
#                    x, y, z = ecl_to_cart(lat, lon, 400)
#                    coords[0].append(x)
#                    coords[1].append(y)
#                    coords[2].append(z)
                Ncount += 1.

#coords = np.array(coords)
lat_lon = np.array(lat_lon)

print 'sphere points ready'
sys.stdout.flush()

'''
fig = plt.figure(figsize=(13,13))
ax = fig.add_subplot(111, projection='3d')
ax.set_aspect("equal")

u = np.linspace(0, 2 * np.pi, 1000)
v = np.linspace(0, np.pi, 1000)
x = np.outer(np.cos(u), np.sin(v))
y = np.outer(np.sin(u), np.sin(v))
z = np.outer(np.ones(np.size(u)), np.cos(v))
ax.plot_surface(x, y, z, color='b', linewidth=0, alpha=.3)
ax.scatter(coords[0,:], coords[1,:], coords[2,:], c='k', s=10, edgecolor='k')
plt.tight_layout()
plt.savefig('sphere_points.png', dpi=200)
#plt.show()
'''

#%% DIVIDING FOR QUICKER RUNNING DURING TESTING PHASES

# DEFINE ORBITAL PROPERTIES FOR EACH SET OF COORDINATES
r = 400

# use the first coordinate set as a test object

distances = []
for lat, lon in zip(lat_lon[0,:], lat_lon[1,:]):
    print 'starting position:', lon, lat
    x, y, z = ecl_to_cart(lat, lon, r)
    v_vec = v_xyz(lat, lon, r, 0)
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
    if lat >= 0.:
        omega = np.arccos(np.dot(Omega_hat, r_hat))*u.rad.to(u.degree)
    else: 
        omega = -np.arccos(np.dot(Omega_hat, r_hat))*u.rad.to(u.degree)


    # USE ORBIT/CCD LOCATING CODE TO IDENTIFY IF IT APPEARED ON DES
    
    #print 'elements set'
    #sys.stdout.flush()

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
    
    '''
    # perform test to see if everything above is working correctly
    ob.date = epoch
    o.compute(ob)
    o_c = astro_topo(o)
    ra = o_c.ra.value
    dec = o_c.dec.value
    print '------------------------------------------------------------------------------------'
    print 'starting :', '%.2f, %.2f'%(lon, lat)
    print 'ephem    :', '%.2f, %.2f'%(ra, dec)
    dist = distance(ra, lon, dec, lat)
    print 'dist     :', '%.2f'%(dist)
    print '------------------------------------------------------------------------------------'
    distances.append(dist)
    '''
    
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
    fig, ax = plt.subplots(figsize=(16,10))
    for m in month_breaks:
        # pick out month from corners.fits
        month_temp = corners[corners['mjd_mid'] >= np.min(m)]
        month = month_temp[month_temp['mjd_mid'] <= np.max(m)]
    
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
#        print 'search_radius :', sep.degree/2. + .17
    
        # build tree to search for near neighbors around mid position
        treedata = zip(month['ra'][:,4]*np.cos(month['dec'][:,4]), month['dec'][:,4])
        tree = spatial.cKDTree(treedata)
    #    near = tree.query_ball_point((ra_mid*np.cos(dec_mid), dec_mid),r=sep.degree/2.+ .17)
        near = tree.query_ball_point((ra_mid*np.cos(dec_mid), dec_mid),r=2)
    #    near = tree.query_ball_point((ra_mid*np.cos(dec_mid), dec_mid),r=8)
        if near == []:
            print 'NO NEAR NEIGHBORS'
            sys.stdout.flush()
            continue
        else:
            month_near = month[near]
    
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
                print '-----------> FOUND CCD'
                ax.add_patch(matplotlib.patches.Rectangle((np.min(ccd['ra']), np.min(ccd['dec'])), \
                                                          np.max(ccd['ra']) - np.min(ccd['ra']), \
                                                          np.max(ccd['dec']) - np.min(ccd['dec']), alpha=0.1))
                plt.scatter(ra, dec, marker='o', c='r', s=50, edgecolors='r')
                sys.stdout.flush()
                overlaps.append(ccd)
        if overlaps == []:
            print 'NO OVERLAP (BUT NEAR NEIGHBORS)'
            sys.stdout.flush()
    print '# overlapping CCDs: ', len(overlaps)
    plt.axis('equal')
    plt.title('Finding 2012 VR113')
    plt.xlabel('ra [deg]')
    plt.ylabel('dec [deg]')
    # plot the whole orbit across the area we care about
    times = []
    print 'plotting whole orbit'
    sys.stdout.flush()
    for row in overlaps:
        times.append(row['mjd_mid'])
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
    plt.show()
    fig.savefig('2-deg-2012_VR113-CCD.png', dpi=100)
    plt.close()
    '''
    with open('p9overlap_test.txt', 'w') as f:
        f.write('expnum, detpos, ra_mid, dec_mid, mjd_mid \n')
        for row in overlaps:
            f.write('%f, %s, %f, %f, %s \n' %(row['expnum'], row['detpos'], row['ra'][4], row['dec'][4], mjd_to_date(row['mjd_mid'])))
    '''
    break






