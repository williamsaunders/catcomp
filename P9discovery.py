import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import astropy.io.fits as fits
import sys
import timeit
from discovery import Totals, Discovery

# import P9 simulation results
p9 = fits.getdata('P9simulation_results/P9results.fits')

# identify the total number of objects with at least 1 detection on a DES CCD
total_objects = Totals(p9)

# perform calculation of objects with minimum unique (different date) detections >= (3,4,5,6,7)

