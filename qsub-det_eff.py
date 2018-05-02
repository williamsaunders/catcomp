import numpy as np
import os
import glob 

zonelist = np.loadtxt('zones_todo4-30.txt')
zones = np.array(zonelist)
z_chunks = np.array_split(zones, 40)
for i, z in enumerate(z_chunks):
    f = open('z' + str(i).zfill(2) + '.txt', 'w')
    for line in z:
        f.write('%d \n'%line)
    f.close()

for i, z in enumerate(z_chunks):
    zstr = 'z' + str(i).zfill(2) + '.txt'
    cmd = 'qsub -l h_vmem=10G -l des -cwd -b y -N Detection%s python detection_efficiency_folio.py '%str(i) + '/data3/garyb/tno/matcher ' + zstr
    os.system(cmd)
