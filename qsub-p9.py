import argparse
import numpy as np
import os

sphere_file = 'sphere_points/sphere_points0.csv'
number = 10
for i in range(number):
    cmd = 'qsub -l h_vmem=10G -l des -S /bin/tcsh -cwd -b y -N P9%d python p9population.py sphere_points/sphere_points1.csv --qsub %d --chunk %d' %(i, number, i)
    os.system(cmd)
