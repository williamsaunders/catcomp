import argparse
import numpy as np
import os

zones = np.loadtxt('zones_todo4-30.txt')
for zone in zones[120:]:
    z = str(int(zone)).zfill(3)
    path = '/data3/garyb/tno/matcher/'
    zonepath = path + 'zone' + z
    print zonepath
    cmd = 'qsub -l h_vmem=20G -l des -cwd -b y -N Comb%s python combine_tiles_folio.py '%z + zonepath
    os.system(cmd)
    


