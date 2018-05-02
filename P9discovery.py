import numpy as np
import matplotlib
matplotlib.use('TkAgg')
matplotlib.rcParams.update({'font.size' : 14})
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import sys
import timeit
from discovery import Totals, Discovery

# import P9 simulation results
p9 = fits.getdata('P9simulation_results/P9results.fits')

# identify the total number of objects with at least 1 detection on a DES CCD
num_detect = []
total_num, num_detect = Totals(p9, num_detect=num_detect)
print total_num, len(num_detect)
print num_detect

# perform calculation of objects with minimum unique (different date) detections >= (3,4,5,6,7)
# and different mag thresholds for each

discovery_dict = {}

discovered = []
count_unique = []
det3, count3 = Discovery(p9, 3, discovered=discovered, count_unique=count_unique)
print len(det3), len(count3)

discovered = []
count_unique = []
det4, count4 = Discovery(p9, 4, discovered=discovered, count_unique=count_unique)
print len(det4), len(count4)

discovered = []
count_unique = []
det5, count5 = Discovery(p9, 5, discovered=discovered, count_unique=count_unique)
print len(det5), len(count5)

discovered = []
count_unique = []
det6, count6 = Discovery(p9, 6, discovered=discovered, count_unique=count_unique)
print len(det6), len(count6)

discovered = []
count_unique = []
det7, count7 = Discovery(p9, 7, discovered=discovered, count_unique=count_unique)
print len(det7), len(count7)

discovered = []
count_unique = []
det8, count8 = Discovery(p9, 8, discovered=discovered, count_unique=count_unique)
print len(det8), len(count8)


# collect and plot statistics
found = np.array([len(det3), len(det4), len(det5), len(det6), len(det7), len(det8)])
fleet_size = 29388.
found_ratio = (found/fleet_size)

fig, ax1 = plt.subplots(figsize=(12,8))
ax1.plot(range(3,9), found_ratio, linewidth=3)
plt.ylim(0,1)
ax1.set_yticks(np.arange(0, 1.1, .1))
ax1.tick_params('y', colors='b')
ax1.xaxis.grid(True)
ax1.yaxis.grid(True, linestyle='dashed')
ax1.set_xlabel('Detection Number Threashold')
ax1.set_ylabel('Frequency Detected on Footprint', color='b')
plt.title('Bright Limit Detection Results')
ax2 = ax1.twinx()
ax2.set_ylabel(r'Effective Area ($\Omega_A$) [sq. deg.]', color='r')
Omega = np.arange(0,1.1,.1)*41253*(1./8)
ax2.set_yticks(Omega)
ax2.tick_params('y', colors='r')
plt.grid(linestyle='dashed')
plt.tight_layout()
plt.savefig('bright_limit_plot.png', dpi=100)
#plt.show()
plt.close('all')

f = open('bright_limit_statistics.txt','w')
f.write('found on footprint, effective area \n')
for a in found_ratio:
    f.write('%.3f, %.3f \n'%(a,41253*(1./8)*a))
f.close()

# now plot the ra, dec of all the objects with >=4 detections
mcolors = []
for name, hex in matplotlib.colors.cnames.iteritems():
    mcolors.append(name)
for i, ob in enumerate(det4):
    if i % 1000 == 0:
        print i, '/', len(det4)
    co = np.random.choice(mcolors)
    plt.scatter(ob[0]['ra'], ob[0]['dec'], c=co, edgecolor=co)
plt.title('Location of Bright Limit P9s (>=4 det)')
plt.xlabel('ra [deg]')
plt.ylabel('dec [deg]')
plt.axis('equal')
plt.savefig('bright_limit_scatter.png', dpi=100)
#plt.show()
plt.close('all')
    

# now do a histogram plot of the number of detections for everything 
plt.figure(figsize=(12,6))
plt.hist(num_detect, bins=np.arange(1,50,1), align='mid', color='b', label='With 6 Hour Detections', alpha=.5)
plt.hist(count3, bins=np.arange(1,50,1), align='mid', color='r', label='Without 6 Hour Detections', alpha=.5)
plt.title('Histogram of Planet Nine Detections')
plt.xlabel('Number of Detections')
plt.ylabel('Number')
plt.grid()
plt.legend()
plt.savefig('detection_histogram.png', dpi=100)
#plt.show()
plt.close('all')

