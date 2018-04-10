import numpy as np
import matplotlib.pyplot as plt

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
        
