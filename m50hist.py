import numpy as np
import glob
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size' : 18})
import astropy.io.fits as fits
import sys

#files = glob.glob('zone_efficiencies/*.fits')
data = fits.getdata('zone_efficiencies/all-coadd_detection_results.fits')
m50s = {'g':[], 'r':[], 'i':[], 'z':[], 'Y':[]}
bands =  ['g', 'r', 'i', 'z', 'Y']
for band in bands:
    m50band = data[data['band']==band]['m50']
    m50s[band].append(m50band)
sys.stdout.flush()
bins = np.arange(19,26,.05)
plt.figure(figsize=(15,8))
num = 0
for band, co in zip(bands, ['b', 'g', 'y', 'orange', 'r']):
    m50 = np.hstack(m50s[band])
    plt.hist(m50, bins=bins, align='mid', histtype='step', normed=True, label=band, color=co, linewidth=2)
    plt.axvline(np.median(m50), color=co, linestyle='dashed', label='%s median = %.2f'%(band, np.median(m50)), linewidth=2)
    num += len(m50)
plt.legend()
plt.xlim(19,26)
plt.ylim(0,2.6)
plt.title('Detections Histogram (%d exposures)'%num)
plt.xlabel('m50 Magnitude')
plt.ylabel('Frequency (normalized)')
plt.grid()
plt.tight_layout()
plt.savefig('ALLzones_histogram.png', dpi=100)
plt.show()
'''



data = np.genfromtxt('zone_efficiencies/zone001-coadd_detection_results.csv', delimiter=',', dtype=None)
data = np.array(data)

band = data[1]
m50 = data[2]
bands =  ['g', 'r', 'i', 'z', 'Y']

for i, band in zip(range(5), ['g', 'r', 'i', 'z', 'Y']):
    print band
    m50s = []
    for row in data:
        if band in row[1]:
            m50s.append(row[2])
    m50s = np.array(m50s, dtype='float')
    ii = np.zeros(len(m50s))
    ii.fill(i)
    avg = np.median(m50s)
    print avg
    plt.scatter(ii, m50s, alpha=.3)
    plt.plot(i, avg, 'r*', markersize=20, alpha=1)
plt.xticks(range(5), ['g', 'r', 'i', 'z', 'Y'])
plt.title('Comparison of Depth of Bands in Zone 001')
plt.ylabel('M50')
plt.xlabel('band')
plt.grid()
#plt.savefig('zone_efficiencies/zone001-band_efficiency.png')
plt.show()


# Now make the histogram within each band
print 'making histogram'

for i, band, co in zip(range(5), ['g', 'r', 'i', 'z', 'Y'], ['r', 'orange', 'y', 'g', 'b']):
    print band
    m50s = []
    for row in data:
        if band in row[1]:
            m50s.append(row[2])
    m50s = np.array(m50s, dtype='float')
    plt.hist(m50s, bins=np.arange(19,26,.1), label=band, fill=True, color=co, alpha=.5)
    ii = np.zeros(len(m50s))
    ii.fill(i)
    avg = np.median(m50s)
    print avg
#    plt.scatter(ii, m50s, alpha=.3)
#    plt.plot(i, avg, 'r*', markersize=20, alpha=1)
plt.legend()
plt.show()


'''


        
