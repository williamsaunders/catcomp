'''
William Saunders 

### THIS PROGRAM IDENTIFIES CCDs WHERE POTENTIAL FAKE TNOs WOULD HAVE LANDED ON THE DES FOOTPRINT ###
'''
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

import astropy.io.fits as fits
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import spatial
import ephem
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.time import Time
import sys

print 'setup complete'
sys.stdout.flush()

# import corners file
loc = 'y4a1.ccdcorners.fits'
corners = fits.getdata(loc)
print 'imported corners.fits'
sys.stdout.flush()

#%% Build and test an EllipticalBody for 1992 QB1

# initialize object and add orbital elements
qb = ephem.EllipticalBody()
'''
#2014 QU441
qb._e = .1061652351675586
qb._a =  41.20233278657857
qb._Om = 113.2497090323359
qb._inc =  19.33917312630782
qb._om = 299.6258241396662
#qb._H = 7.1
#qb._G = .150
qb._M = 341.3857129245001
qb._epoch_M = '2015/05/27 00:00:00'
qb._epoch = '1999/12/31 12:00:00'
'''
#2012 VR113
qb._e = 0.1691997
qb._a =  47.4449441
qb._Om = 121.01017
qb._inc = 19.28652
qb._om = 219.23707
qb._M = 41.62698
qb._epoch_M = '2017/09/04 00:00:00'
qb._epoch = '1999/12/31 12:00:00'


# observer object for CTIO that is same for all observations
ob = ephem.Observer()
ob.lat = -30.169117
ob.lon = 289.194100
ob.elevation = 2389
ob.pressure = 0

print 'observer set'
sys.stdout.flush()

#%% Now divide up all the data into months
remainder = np.remainder(len(corners['mjd_mid']), np.floor(len(corners['mjd_mid'])/40))
breaks = np.zeros(41)
for i, b in enumerate(breaks):
    breaks[i] = i*np.floor(len(corners['mjd_mid'])/40.)
breaks[-1] = breaks[-1] + remainder
breaks = np.array(breaks, dtype=int)
breaks = breaks[1:-1]
month_breaks = np.split(corners['mjd_mid'], breaks)

#%% For each observation, determine which month we care about
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
    qb.compute(ob)
    qb_c = astro_topo(qb)
    ra_mid = qb_c.ra.value
    dec_mid = qb_c.dec.value  

    # to determine radius of tree search, see how far the object moves in the month
    month_start = np.min(m)
    ob.date = mjd_to_date(month_start)
    qb.compute(ob)
    qb_c_start = astro_topo(qb)

    month_end = np.max(m)
    ob.date = mjd_to_date(month_end)
    qb.compute(ob)
    qb_c_end = astro_topo(qb)

    sep = qb_c_start.separation(qb_c_end)
    print 'search_radius :', sep.degree

    # build tree to search for near neighbors around mid position
    treedata = zip(month['ra'][:,4]*np.cos(month['dec'][:,4]), month['dec'][:,4])
    tree = spatial.cKDTree(treedata)
    near = tree.query_ball_point((ra_mid*np.cos(dec_mid), dec_mid),r=sep.degree*1.1)
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
        qb.compute(date)
        qb_c = astro_topo(qb)
        ra = qb_c.ra.value
        dec = qb_c.dec.value  

        if ((np.min(ccd['ra']) <= ra) and (np.max(ccd['ra']) >= ra) and \
        (np.min(ccd['dec']) <= dec) and (np.max(ccd['dec']) >= dec)):
            print '-----------> FOUND CCD'
            ax.add_patch(matplotlib.patches.Rectangle((np.min(ccd['ra']), np.min(ccd['dec'])), np.max(ccd['ra']) - np.min(ccd['ra']),  np.max(ccd['dec']) - np.min(ccd['dec']), alpha=0.1))
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
    qb.compute(ob)
    qb_c = astro_topo(qb)
    all_ra.append(qb_c.ra.value)
    all_dec.append(qb_c.dec.value)
#    plt.scatter(qb_c.ra.value, qb_c.dec.value, marker='o', c='k', edgecolors='k', s=3)
plt.plot(all_ra, all_dec, c='k', linewidth=2)
plt.show()
fig.savefig('TEST2012_VR113-CCD.png', dpi=100)
plt.close() 
print overlaps


