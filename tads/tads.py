# tads.py
# -------------------------------------------
# testing random stuff for TADs

import sys
import random
import numpy as np
import collections
import copy
import scipy.io as sio
from scipy import integrate
import matplotlib.pyplot as plt
import csv
from utils import *
from constants import *
import scipy.signal as sig

# windows of 50bp on either side of diagonal
def find_windows(chr_mat, size=WIN_SIZE):

    window_lst = []

    for i, row in enumerate(chr_mat):
        window = row[max(i-size, 0) : min(i+size+1, len(row))]
        if i-size < 0:
            window = np.concatenate((np.zeros(size - i), window), axis=0)
        if i+size+1 > len(row):
            window = np.concatenate(((np.zeros(size+i+1 - len(row))), window), axis=0)
        window_lst.append(window)

    return window_lst

# checks a row vector for 0 reads, default 75% to pass
def check_reads(vec_lst, threshold=QUAL_THRESH):
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
def calc_L(vec, size=WIN_SIZE):

    return sum(vec[:size])

# calculate R (sum of R)
def calc_R(vec, size=WIN_SIZE):

    return sum(vec[size:])

# calculate L+R
def calc_LpR(L, R):

    return (L + R)

# calculate L-R
def calc_LmR(L, R):

    return (L - R)

# histogram of read quality for all chromosomes
def plot_all_qual():
    for chr_num in range(1, NUM_CHRMS+1):
        chr_mat = load_chr(chr=chr_num)
        window_lst = find_windows(chr_mat)
        plot_read_qual(window_lst, str(chr_num))

# linear plot of measures
def plot_lst(lst, xlabel, ylabel, num):
    fig = plt.figure(0)

    plt.plot(lst)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title('Chromosome ' + str(num))
    fig.savefig(OUTPATH + 'chr' + str(num) + '-' + ylabel + '.png')
    plt.close(0)

# calculates all metric values
def calc_metrics(chr_num):
    # load data
    chr_mat = load_chr(chr=chr_num)
    window_lst = find_windows(chr_mat)

    # calculate L, R, L+R
    length = len(window_lst)
    L_lst = np.zeros(length)
    R_lst = np.zeros(length)
    LpR_lst = np.zeros(length)
    LmR_lst = np.zeros(length)

    for i, vec in enumerate(window_lst):
        L = calc_L(vec)
        R = calc_R(vec)
        L_lst[i] = L
        R_lst[i] = R
        LpR_lst[i] = calc_LpR(L, R)
        LmR_lst[i] = calc_LmR(L, R)

    return (L_lst, R_lst, LpR_lst, LmR_lst)

# finds shifts using LmR
def find_shifts(lst, delta=SHIFT_THRESH):
    data = np.array(lst)

    # get relative min and max
    min_lst = sig.argrelmin(data)[0]
    max_lst = sig.argrelmax(data)[0]

    # store shift (aka boundary) locations
    shift_pts = []

    # pointer for max_lst
    i = 0

    # keep track of confirmed prev min and max
    pmin = min_lst[0]
    pmax = max_lst[0]

    while i < len(max_lst) and min_lst != []:
        xmax = max_lst[i]

        # find first rel min after current rel max
        index_min = np.argmax(min_lst > xmax)
        xmin = min_lst[index_min]

        dmin = lst[xmin]
        dmax = lst[xmax]

        # print 'xmax: ', xmax
        # print 'xmin: ', xmin
        # print 'dmax: ', dmax
        # print 'dmin: ', dmin
        # print 'dff: ', dmax - dmin
        # print 'pmin: ', pmin

        if dmax < (lst[min_lst[index_min - 1]] + DEL_MIN):
            i += 1
        elif (dmax - dmin) > delta:

            # check if prev shift had the same min
            if pmin == dmin:
                if (dmax - dmin) > (pmax - pmin):
                    shift_pts.pop()
                else:
                    i += 1
                    continue

            shift_pts.append((xmin + xmax) / 2)
            i += 1
            pmin = dmin
            pmax = dmax
            # print 'yes'
        else:
            # sets distance limit between max and min
            if (xmin - xmax) > MIN_DIST:
                i += 1

            else:

                if len(min_lst) <= index_min:
                    next_min = min_lst[index_min + 1]

                    # check max's in between

                    next_max_lst = list(enumerate(max_lst))

                    for j, m in next_max_lst[i:]:
                        if m < next_min:
                            if (lst[m] - lst[next_min]) > (lst[xmax] - lst[next_min]):
                                i = j
                                continue

                # discards first part of lst
                min_lst = min_lst[index_min + 1:]

    bounds = np.zeros(len(lst))

    for shift in shift_pts:
        bounds[shift] = 1

    return np.array(shift_pts), bounds

# saves and plots metrics
def process_metrics(L_lst, R_lst, LpR_lst, LmR_lst, chr_num):

    for chr_num in range(1, NUM_CHRMS+1):

        L_lst, R_lst, LpR_lst, LmR_lst = calc_metrics(chr_num)

        # save and plot L, R, L+R, L-R graphs
        sio.savemat(OUTPATH + 'chr' + str(chr_num) + '-LR.mat', {'L': L_lst, 'R': R_lst, 'LpR': LpR_lst, 'LmR': LmR_lst})
        
        plot_lst(L_lst, 'Position', 'L', chr_num)
        plot_lst(R_lst, 'Position', 'R', chr_num)
        plot_lst(LR_lst, 'Position', 'R + L', chr_num)
        plot_lst(LmR_lst, 'Position', 'R - L', chr_num)

# saves and plots shifts for each chromosome
def process_shifts():
    writer = csv.writer(open(OUTPATH + 'LRbounds-output.csv', 'wb'))

    for chr_num in range(1, NUM_CHRMS+1):

        chrm = 'chr' + str(chr_num)

        L_lst, R_lst, LpR_lst, LmR_lst = calc_metrics(chr_num)

        lst = LmR_lst
        shift_lst, bound_lst = find_shifts(lst)

        writer.writerow([chrm])
        writer.writerow(['shifts'])
        writer.writerow(shift_lst)
        writer.writerow(['bounds'])
        writer.writerow(bound_lst)

        f, axarr = plt.subplots(2, sharex=True)
        axarr[0].plot(lst)
        axarr[1].plot(bound_lst)

        f.savefig(OUTPATH + chrm + '-LRbound-py.png')
        sio.savemat(OUTPATH + 'TADbounds_' + chrm + '.mat', {chrm + 'shifts': shift_lst, chrm + 'bounds': bound_lst})
    
chr_num = 1
chrm = 'chr' + str(chr_num)
chr_mat = load_chr(chr=chr_num)
window_lst = find_windows(chr_mat)
windows = list(enumerate(window_lst))
f, axarr = plt.subplots(5, sharex=True)
plt_num = 0
for i, w in windows[1000:1005]:
    axarr[plt_num].plot(w)
    axarr[plt_num].set_xlabel(chrm + '-pos' + str(i))
    plt_num += 1

plt.show()









