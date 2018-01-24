import os 
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('parallel_num', type=int)
args = parser.parse_args()

for i in range(0, args.parallel_num+1):
    os.system('qsub -l h_vmem=5G -cwd -b y -N python%d python detection_efficiency_folio.py --num %d --i %d'%(i, args.parallel_num, i))
