def Totals(p9, num_detect=[]):
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
    print 'detthreash = %d, magthresh = %d' %(detthresh, magthresh)
    import numpy as np
    from astropy.io import fits
    from astropy.time import Time
    import sys
    import timeit
    import random 
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
        # This is the big one, the function that determines the likelihood of observing a 
        # P9 of different magnitudes.  It takes magnitude and the exposure number and returns
        # a boolean of whether the observation is "realized" 
                
        word code:
        look up the detection statistics for that exposure 
        determine at that magnitude what the likelihood of observing is 
        do a draw from a random number generator to see if that observation is realized
        return realized = True or False

        expstats = detstats[detstats['expnum']==expnum]
        if len(expstats) > 1 :
            sys.exit('MORE THAN ONE EXPSURE')
        if len(expstats) < 1:
            print 'EXPSURE NOT FOUND'
            return False
        
        detprob = expstats['c']/(1+np.exp(expstats['k']*(mag-expstats['m50'])))
        print detprob
        if detprob >= np.random.random():
            return True
        else:
            return False
        
    
    # determine the total number of objects detected once on DES CCDs

    print 'doing discovery'
    count = 0
    for i, ob in enumerate(p9):
        if i % 10000 == 0:
            print i, '/', len(p9)
        if i == 0:
            collect = []
            dates = []
            # for each object, see if the observation is realized at that magthresh, expnum
            if magthresh != 99: 
                realize = MagDetection(magthresh, ob['expnum'])
                if realize:
                    dates.append(ob['date'])
                    collect.append(ob)
                    print 'reazlied'
                else:
                    print 'not realized'
        elif ob['ob_num'] == p9[i-1]['ob_num']:
            if magthresh != 99: 
                realize = MagDetection(magthresh, ob['expnum'])
                if realize:
                    dates.append(ob['date'])
                    collect.append(ob)
                    print 'reazlied'
                else:
                    print 'not realized'
        else:
            # first remove duplicate dates (can't count 2 detections within 6 hours)
            # NOTE: by this point only objects considering are those realized at magthresh
            dates = np.array(dates)
            dates_rounded = np.round(dates*4.)/4.
            unique_dates, count_dates = np.unique(dates_rounded, return_counts=True)
            count_unique.append(len(unique_dates))
            # second use only realized observations to determine unique seasons
            seasons = []
            for date in unique_dates:
                seasons.append(season(date))
            seasons = np.array(seasons, dtype=str)
            unique_seasons, count_seasons = np.unique(seasons, return_counts=True)
            detect = np.greater_equal(count_seasons, np.full(len(count_seasons), \
                                                             detthresh))
            if np.any(detect):
                discovered.append(collect)
                # NOTE: DUMPING THE collect OBJECT INTO discovered ISN'T QUITE RIGHT BECAUSE
                # SOME OF THOSE OBSERVATION WERE NOT REALIZED EITHER FOR MAGNITUDE, 6-HOUR,
                # OR DIFFERENT SEASON REASONS. THE POINT IS TO LOOK AT THE NUMBER OF TIMES 
                # ANY OBJECTS GETS A collect DUMPED INTO discovered. IF ALL THE OBSERVATIONS
                # ARE DUMPED, IT MEANS AT LEAST DETTHRESH NUMBER OF THEM FIT THE CRITERIA
                count +=1

            dates = []
            collect = []
            if magthresh != 99: 
                realize = MagDetection(magthresh, ob['expnum'])
                if realize:
                    dates.append(ob['date'])
                    collect.append(ob)
                    print 'reazlied'
                else:
                    print 'not realized'
    return discovered, count_unique


