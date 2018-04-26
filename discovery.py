'''
THIS IS THE DISCOVERY STEP OF THE PLANET NINE FLEET. 
DISCOVERY CONDITIONS:
1. At least 4 detections in the same season on different dates 
2. Determine probability of detecting based on exposure detection function <-- not yet!
'''
# this code defines a function that should be imported and run by a separate driver

    # import P9 simulation results
    p9 = fits.getdata('P9simulation_results/P9results.fits')
    p9det = p9[p9['num']>=4]


def Totals(p9):
    print 'doing totals'
    total_objects = 0
    for i, ob in enumerate(p9):
        if i % 10000 == 0:
            print i, '/', len(p9)
            if i == 0:
                total_objects += 1
            elif ob['ob_num'] == p9[i-1]['ob_num']:
                continue
            else:
                total_objects += 1
    return total_objectsx

def Discovery(p9, detthresh, magthresh=99):
    import numpy as np
    from astropy.io import fits
    from astropy.time import Time
    import sys
    import timeit
    import random 
    def mjd_to_date(mjd):
        time = Time(mjd, format='mjd')
        return time.iso
    def date_to_mjd(date):
        time = Time(date, format='iso')
        return time.mjd
    def season(mjd):
        # identify season of detections
        s1 = date_to_mjd('2013-01-01 12:00:00')
        s2 = date_to_mjd('2014-01-01 12:00:00')
        s3 = date_to_mjd('2015-01-01 12:00:00')
        s4 = date_to_mjd('2016-01-01 12:00:00')
        s5 = date_to_mjd('2017-01-01 12:00:00')
        if mjd > s1 and mjd < s2:
            return 's1'
        elif mjd > s2 and mjd < s3:
            return 's2'
        elif mjd > s3 and mjd < s4:
            return 's3'
        elif mjd > s4 and mjd < s5:
            return 's4'
        elif mjd > s5:
            return 's5???'
    
    # determine the total number of objects detected once on DES CCDs

    print 'doing discovery'
    count = 0
    discoverd = []
    for i, ob in enumerate(p9det):
        if i % 10000 == 0:
            print i, '/', len(p9det)
        if i == 0:
            collect = []
            dates = []
            dates.append(ob['date'])
            collect.append(ob)
        elif ob['ob_num'] == p9det[i-1]['ob_num']:
            dates.append(ob['date'])
            collect.append(ob)
        else:
            # first remove duplicate dates (can't count 2 detections within 6 hours)
            dates = np.array(dates)
            dates_rounded = np.round(dates*4.)/4.
            unique_dates, count_dates = np.unique(dates_rounded, return_counts=True)
            # second use only unique dates to determine unique seasons
            seasons = []
            for date in unique_dates:
                seasons.append(season(date))
            seasons = np.array(seasons, dtype=str)
            unique_seasons, count_seasons = np.unique(seasons, return_counts=True)
            detect = np.greater_equal(count_seasons, np.full(len(count_seasons), \
                                                             detthresh))
            if np.any(detect):
                discovered.append(collect)
                count +=1
            dates = []
            collect = []
            dates.append(ob['date'])
            collect.append(ob)
    return discovered


