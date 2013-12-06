# tads.py
# -------------------------------------------
# testing random stuff for TADs

import sys
import random
import numpy as np
import collections
import copy
from scipy import integrate
import matplotlib.pyplot as plt
import csv
from utils import *
from constants import *

# windows of 50bp on either side of diagonal
def find_windows(chr_mat, size=50):

    window_lst = []

    for i, row in enumerate(chr_mat):
        window = row[max(i-size, 0) : min(i+size+1, len(row))]
        if i-size < 0:
            window = np.concatenate((np.zeros(size - i), window), axis=0)
        if i + i+size+1 > len(row):
            window = np.concatenate(((np.zeros(size+i+1 - len(row))), window), axis=0)
        window_lst.append(window)

    return window_lst

# checks a row vector for 0 reads, default 75% to pass
def check_reads(vec_lst, threshold=0.75):
    lst = []
    
    for vec in vec_lst:
        cutoff = threshold * len(vec)
        if sum(vec != 0) > threshold:
            lst.append(vec)

    # list of passed vectors
    return lst

# plot read quality over region
def plot_read_qual(vec_lst, num):
    qual_lst = []

    # calculate quality
    for vec in vec_lst:
        qual = sum(vec != 0) / float(len(vec))
        qual_lst.append(qual)

    # plot histogram
    fig = plt.figure(0)

    plt.hist(qual_lst, bins=100)
    plt.xlabel('Read Quality (portion of non-zero reads)')
    plt.ylabel('Count')
    plt.title('Chromosome ' + num)
    fig.savefig(OUTPATH + '-chr-' + num + '-qual.png')
    plt.close(0)

# calculate L (sum of left)
def calc_L(vec, size):

    return sum(vec[:size])

# calculate R (sum of R)
def calc_R(vec, size):

    return sum(vec[size:])



# histogram of read quality for all chromosomes
def plot_all_qual():
    for chr_num in range(1, NUM_CHRMS+1):
        chr_mat = load_chr(chr=chr_num)
        window_lst = find_windows(chr_mat)
        plot_read_qual(window_lst, str(chr_num))






