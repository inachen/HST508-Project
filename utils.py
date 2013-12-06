# utils.py
# -------------------------------------------
# util functions for Hi-C 

import sys
import random
import numpy as np
import collections
import copy
from scipy import integrate
import scipy.io as sio
import matplotlib.pyplot as plt
import csv
from constants import *


# loads data for a chromosome 
def load_chr(cell='HIMR',chr=1,corrected=True):

    # file path
    chr_path = DATAPATH + cell + "-all_res40000/" + "chr" + str(chr) + "_chr" + str(chr)
    if corrected:
        chr_path += "_corrected"
    chr_path += ".mat"

    # return corrected array
    if corrected:
        chr_data = sio.loadmat(chr_path)['corrected']

    # sio.savemat('test.mat', {'matrix':chr_data})

    return chr_data

# write data to csv file
def csv_mat(mat, outfile):
    writer = csv.writer(outfile, 'wb')
    for row in chr_data:
       writer.writerow(row)
