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
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
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

# calculate half R
def calc_hR(vec, size=WIN_SIZE):

    return sum(vec[size : int(size * 1.5)])

# calculate half L
def calc_hL(vec, size=WIN_SIZE):

    return sum(vec[size/2 : size])

# calculate quarter R
def calc_qR(vec, size=WIN_SIZE):

    return sum(vec[size : int(size * 1.25)])

# calculate quarter L
def calc_qL(vec, size=WIN_SIZE):
    
    return sum(vec[int(size * 0.75) : size])

# calculate L+R
def calc_LpR(L, R):

    return (L + R)

# calculate L-R
def calc_LmR(L, R):

    return (L - R)

# calculate R-L
def calc_RmL(L, R):

    return (R - L)

# calculate LdR
def calc_LdR(L, R):

    return (L/R)

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
    # LpR_lst = np.zeros(length)
    LmR_lst = np.zeros(length)
    # RmL_lst = np.zeros(length)
    hL_lst = np.zeros(length)
    hR_lst = np.zeros(length)
    qL_lst = np.zeros(length)
    qR_lst = np.zeros(length)
    LdR_lst = np.zeros(length)
    qLdR_lst = np.zeros(length)

    for i, vec in enumerate(window_lst):
        L = calc_L(vec)
        R = calc_R(vec)
        hL = calc_hL(vec)
        hR = calc_hR(vec)
        qL = calc_qL(vec)
        qR = calc_qR(vec)
        L_lst[i] = L
        R_lst[i] = R
        hL_lst[i] = hL
        hR_lst[i] = hR
        qL_lst[i] = qL
        qR_lst[i] = qR
        # LpR_lst[i] = calc_LpR(L, R)
        LmR_lst[i] = calc_LmR(L, R)
        # RmL_lst[i] = calc_RmL(L, R)
        LdR_lst[i] = calc_LdR(L, R)
        qLdR_lst[i] = calc_LdR(qL, qR)

    return (L_lst, R_lst, hL_lst, hR_lst, qL_lst, qR_lst, LmR_lst, LdR_lst, qLdR_lst)

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

# finds shifts using RmL
def find_shifts_RmL(lst, delta=SHIFT_THRESH):
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

    while i < len(min_lst) and max_lst != []:
        xmin = min_lst[i]

        # find first rel min after current rel max
        index_max = np.argmax(max_lst > xmin)
        xmax = max_lst[index_max]

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
            if pmax == dmax:
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

                if len(max_lst) <= index_max:
                    next_max = max_lst[index_max + 1]

                    # check min's in between

                    next_min_lst = list(enumerate(min_lst))

                    for j, m in next_min_lst[i:]:
                        if m < next_max:
                            if (lst[next_max] - lst[m]) > (lst[next_max] - lst[xmin]):
                                i = j
                                continue

                # discards first part of lst
                min_lst = min_lst[index_min + 1:]

    bounds = np.zeros(len(lst))

    for shift in shift_pts:
        bounds[shift] = 1

    return np.array(shift_pts), bounds

# saves and plots metrics
def process_metrics():

    for chr_num in range(1, 2):#NUM_CHRMS+1):

        L_lst, R_lst, hL_lst, hR_lst, qL_lst, qR_lst, LmR_lst, LdR_lst, qLdR_lst = calc_metrics(chr_num)

        # save and plot L, R, L+R, L-R graphs
        sio.savemat(OUTPATH + 'chr' + str(chr_num) + '-LR.mat', 
            {'L': L_lst, 'R': R_lst, 'hL': hL_lst, 'hR': hR_lst, 'qL': qL_lst, 'qR': qR_lst, 'LmR': LmR_lst, 'LdR': LdR_lst, 'qLdR': qLdR_lst})
        
        # plot_lst(L_lst, 'Position', 'L', chr_num)
        # plot_lst(R_lst, 'Position', 'R', chr_num)
        # plot_lst(hL_lst, 'Position', 'hL', chr_num)
        # plot_lst(hR_lst, 'Position', 'hR', chr_num)

# saves and plots shifts for each chromosome
def process_shifts():
    writer = csv.writer(open(OUTPATH + 'LRbounds-output.csv', 'wb'))

    for chr_num in range(1, NUM_CHRMS+1):

        chrm = 'chr' + str(chr_num)

        L_lst, R_lst, LpR_lst, LmR_lst, RmL_lst = calc_metrics(chr_num)

        lst = RmL_lst
        shift_lst, bound_lst = find_shifts(lst)

        writer.writerow([chrm])
        writer.writerow(['shifts'])
        writer.writerow(shift_lst)
        writer.writerow(['bounds'])
        writer.writerow(bound_lst)

        f, axarr = plt.subplots(2, sharex=True)
        axarr[0].plot(lst)
        axarr[1].plot(bound_lst)

        f.savefig(OUTPATH + chrm + '-RmLbound-py.png')
        sio.savemat(OUTPATH + 'TAD-RmL-bounds_' + chrm + '.mat', {chrm + 'shifts': shift_lst, chrm + 'bounds': bound_lst})

# plot window frames
def plot_windows():  
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

chr_num = 1

# metric_lst = zip(L_lst, R_lst, hL_lst, hR_lst)

# plot 

# L_lst, R_lst, hL_lst, hR_lst, qL_lst, qR_lst, LmR_lst, LdR_lst = calc_metrics(chr_num)
# shift_lst, bound_lst = find_shifts(LmR_lst)

# length = len(shift_lst)
# LdR_shift = np.zeros(length)
# LmR_shift = np.zeros(length)
# qL_shift = np.zeros(length)
# qR_shift = np.zeros(length)
# L_shift = np.zeros(length)
# R_shift = np.zeros(length)
# bound_shift = np.zeros(length)

# for i, shift in enumerate(shift_lst):
#     LdR_shift[i] = LdR_lst[shift]
#     L_shift[i] = L_lst[shift]
#     R_shift[i] = R_lst[shift]
#     qL_shift[i] = qL_lst[shift]
#     qR_shift[i] = qR_lst[shift]
#     bound_shift[i] = bound_lst[shift]

# print bound_shift

# fig = plt.figure()
# ax = fig.add_subplot(111)
# colorlist = ['r', 'g', 'b', 'y']
# ax.hist([LdR_lst, LdR_shift], 50, range=[0, 5])
# # ax = fig.add_subplot(122)
# # ax.scatter(qL_shift, qR_shift, c='r')
# plt.title("Test Cluster")
# ax.set_xlabel('LdR')
# ax.set_ylabel('fold')
# # ax.set_zlabel('qR')
# plt.show()

process_metrics()





