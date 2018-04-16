import numpy as np
import os
import glob 

zonelist = np.loadtxt('zones_done.txt')
zones = np.array(zonelist)
z1, z2, z3, z4, z5 = np.array_split(zones, 5)
for z, zs in zip([z1,z2,z3,z4,z5], ['z1','z2','z3','z4','z5']):
    f = open(zs + '.txt', 'w')
    for line in z:
        f.write('%d \n'%line)
    f.close()

for z in ['z1','z2','z3','z4','z5']:
    cmd = 'qsub -l h_vmem=10G -l des -cwd -b y -N Detection%s python detection_efficiency_folio.py '%z[-1] + '/data3/garyb/tno/matcher ' + z + '.txt'
    os.system(cmd)
