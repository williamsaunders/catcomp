'''
This program is two functions to be used bo P9discovery.py to do the discovery step of 
the Planet Nine search
'''

def Totals(p9, num_detect=[]):
    '''
    Determine the total number of objects observed. 
    Takes the p9 fits table.
    Returns the number of total objects and an the object numbers
    (this won't be necessary when the P9 orbit simulation is run again, 
    because the object IDs will be properly unique.  They weren't on my 1 run).
    '''
    import numpy as np
    print 'doing totals'
    total_objects = 0
    for i, ob in enumerate(p9):
        if i % 10000 == 0:
            print i, '/', len(p9)
        if i == 0:
            total_objects += 1
            num_detect.append(ob['num'])
        elif ob['ob_num'] == p9[i-1]['ob_num']:
            continue
        else:
            num_detect.append(ob['num'])
            total_objects += 1
    return total_objects, np.array(num_detect)

def Discovery(p9, detthresh, magthresh=99, discovered=[], count_unique=[]):
    '''
    The actual discovery step for P9.
    Takes the P9 fits table, a number of observations required for detection, 
    the magnitude threashold of Planet Nine, and two blank lists (not sure why Python
    required these.  Threw errors and I didn't have time to fix so I left them)
    '''
    import numpy as np
    from astropy.io import fits
    from astropy.time import Time
    import sys
    import timeit
    import random 
    
    # load the point-source detection data from all exposures
    detstats = fits.getdata('zone_efficiencies/all-coadd_detection_results.fits')

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

    def MagDetection(mag, expnum):
        '''
        This function determines the likelihood of observing a P9 of different magnitudes.  
        It takes magnitude and the exposure number and returns
        a boolean of whether the observation is "realized"
        '''
        expstats = detstats[detstats['expnum']==expnum]
        if len(expstats) > 1 :
            sys.exit('MORE THAN ONE EXPSURE')
            # if this happens, than the exposure point-source data was combined wrong
        if len(expstats) < 1:
            # there were few enough to just return True, but in the future won't need this.
            return (True, expnum)
        
        # determine if magnitude is too low to see for one random draw
        detprob = expstats['c']/(1+np.exp(expstats['k']*(mag-expstats['m50'])))
        if detprob >= np.random.random():
            return (True, 0)
        else:
            return (False, 0)
    
    # determine the total number of objects detected once on DES CCDs
    missing_exp = []
    count = 0
    for i, ob in enumerate(p9):
        realize = 0 
        del realize
        if i % 10000 == 0:
            print i, '/', len(p9)
        if i == 0:
            # because of the way I wrote my P9 simulation output, each object is recorded
            # the number of times it was observed, so I had to write this clunky 
            # way of identifying the number of detections 
            collect = []
            dates = []
            # for each object, see if the observation is realized at that magthresh, expnum
            # NOTE: magnitude of 99 constitutes the bright limit
            if magthresh != 99.: 
                realize = MagDetection(magthresh, ob['expnum'])
            if magthresh == 99.:
                realize = (True, 0)
            if realize[0]:
                dates.append(ob['date'])
                collect.append(ob)
            if realize[1] != 0.:
                missing_exp.append(realize[1])
        elif ob['ob_num'] == p9[i-1]['ob_num']:
            # if the next object in the list is the same as the previous
            if magthresh != 99: 
                realize = MagDetection(magthresh, ob['expnum'])
            if magthresh == 99:
                realize = (True, 0)
            if realize[0]:
                dates.append(ob['date'])
                collect.append(ob)
            if realize[1] != 0.:
                missing_exp.append(realize[1])

        else:
            # if the next object in the last is different --> new object
            # only those who passed the magnitude test make it this far

            # first remove duplicate dates (can't count 2 detections within 6 hours)
            dates = np.array(dates)
            dates_rounded = np.round(dates*4.)/4.
            unique_dates, count_dates = np.unique(dates_rounded, return_counts=True)
            count_unique.append(len(unique_dates))

            # second determine if there are enough detections in one season to discover
            seasons = []
            for date in unique_dates:
                seasons.append(season(date))
            seasons = np.array(seasons, dtype=str)
            unique_seasons, count_seasons = np.unique(seasons, return_counts=True)
            detect = np.greater_equal(count_seasons, np.full(len(count_seasons), \
                                                             detthresh))
            if np.any(detect): # doesn't matter which season
                discovered.append(collect)
                count +=1
            
            dates = []
            collect = []
            if magthresh != 99: 
                realize = MagDetection(magthresh, ob['expnum'])
            if magthresh == 99:
                realize = (True, 0)
            if realize[0]:
                dates.append(ob['date'])
                collect.append(ob)
            if realize[1] != 0.:
                missing_exp.append(realize[1])
    return discovered, count_unique, missing_exp


