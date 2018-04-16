import numpy as np
import glob
import matplotlib.pyplot as plt

files = glob.glob('zone_efficiencies/*')
m50s = {'g':[], 'r':[], 'i':[], 'z':[], 'Y':[]}
bands =  ['g', 'r', 'i', 'z', 'Y']
for file in files:
    data = np.genfromtxt(file, delimiter=',', dtype=None, skip_header=1)
    for band in bands:
        for d in data:
            if band in d[1]:
                m50s[band].append(d[2])
bins = np.arange(19,26,.1)
plt.figure(figsize=(15,8))
num = 0
for band, co in zip(bands, ['b', 'g', 'y', 'orange', 'r']):
    plt.hist(m50s[band], bins=bins, align='mid', histtype='step', normed=True, label=band, color=co, linewidth=2)
    plt.axvline(np.median(m50s[band]), color=co, linestyle='dashed', label='%s median'%band, linewidth=2)
    num += len(m50s[band])
plt.legend()
plt.xlim(19,26)
plt.ylim(0,2.6)
plt.title('10 Zones Detections Histogram (n = %d)'%num)
plt.xlabel('m50 mag')
plt.ylabel('frequency (normalized)')
plt.grid()
plt.savefig('10zones_histogram.png', dpi=100)
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


        
