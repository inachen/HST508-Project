# seq.py
# -------------------------------------------
# handling sequences in HI-C data

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
# file Josh wrote to process kmers
# import kmer
from sklearn import svm

def read_csv(infile):
    f = open(infile, 'rb')
    reader = csv.reader(f)
    mer_lst = []
    bound_count = []
    int_count = []
    baseline_count = []
    for row in reader:
        mer_lst.append(row[0])
        bound_count.append(int(row[1]))
        int_count.append(int(row[2]))
        baseline_count.append(int(row[3]))

    mers7 = {'mers': mer_lst, 'bounds': bound_count, 'ints': int_count,
        'baselines': baseline_count}
    return mers7

def consolidate(mers7):

    mers6 = {'mers': [], 'bounds': [], 'ints': [], 'baselines': []}

    switch = False

    for j, mer7 in enumerate(mers7['mers']):
        for i, mer6 in enumerate(mers6['mers']):
            m1 = mer6[0][1:]
            m2 = mer6[0][:-1]
            if (m1 in mer7) or (m2 in mer7):
                mers6['mers'][i].append(mer7)
                mers6['bounds'][i] += mers7['bounds'][j]
                mers6['ints'][i] += mers7['ints'][j]
                mers6['baselines'][i] += mers7['baselines'][j]
                switch = True
                break

        if switch:
            switch = False
            continue

        mers6['mers'].append([mer7])
        mers6['bounds'].append(mers7['bounds'][j])
        mers6['ints'].append(mers7['ints'][j])
        mers6['baselines'].append(mers7['baselines'][j])

    outfile = OUTPATH + '6mers.csv'
    writer = csv.writer(open(outfile, 'wb'))
    writer.writerow(mers6.keys())
    for row in zip(mers6['mers'], mers6['bounds'], mers6['ints'], mers6['baselines']):
        writer.writerow(row)

    return outfile

def to_6mers():
    infile = 'data/mer-tracks/bd_vs_intHES7mersFalse.csv'
    mers7 = read_csv(infile)
    outfile = consolidate(mers7)
    kmer.validate_enrichment2(outfile, k=6)

def read_kmer_csv(infile, chr=None):
    # read file
    f = open(infile, 'rb')
    reader = csv.reader(f)

    # kmer vectors
    vec_lst = []

    # identification, 1 for boundary, 0 for interior
    class_lst = []

    chr_lst = []
    gc_lst = []

    rownum = 0

    for row in reader:
        if rownum == 0:
            rownum += 1
            continue
        class_lst.append(int(row[0]))
        chr_lst.append(int(row[1]))
        gc_lst.append(float(row[3]))
        vec_lst.append(map(int, row[4:]))

    if chr != None:
        chr_index = [i for i, y in enumerate(chr_lst) if y == chr]
        chr_class = [class_lst[x] for x in chr_index]
        chr_gc = [gc_lst[x] for x in chr_index]
        chr_vec = [vec_lst[x] for x in chr_index]

        return chr_class, chr_gc, chr_vec

    return class_lst, gc_lst, vec_lst


def svm_kmers(infile):

    gc_max = max(gc_lst)

    class_lst, gc_lst, vec_lst = read_kmer_csv(infile)

    # round to nearest 5
    gc_start = int(5 * round(float(min(gc_lst))/5))

    gc_bins = [x for x in range(gc_start, gc_max) if x % 5 == 0]

    for c, v, g in zip(class_lst, vec_lst, gc_lst):
        if 

    # clf = svm.SVC()

    # num_train = 10000
    # num_test = 5000
    # samples = random.sample(zip(class_lst, vec_lst), num_test + num_train)

    # train_class, train_vec = zip(*samples[:num_train])
    # test_class, test_vec = zip(*samples[num_train:])

    # clf.fit(train_vec, train_class)

    # test_predict = clf.predict(test_vec)
    
    # correct = 0

    # for i, t in enumerate(test_predict):
    #     if t == test_class[i]:
    #         correct += 1

    # return correct/float(len(test_class))

num_trials = 10
# for i in range(num_trials):
svm_kmers('data/mer-tracks/top_7mer_vecsGC.csv')


