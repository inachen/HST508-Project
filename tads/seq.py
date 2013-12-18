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

    writer = csv.writer(open(OUTPATH + '6mers.csv', 'wb'))
    writer.writerow(mers6.keys())
    for row in zip(mers6['mers'], mers6['bounds'], mers6['ints'], mers6['baselines']):
        writer.writerow(row)


infile = 'data/mer-tracks/bd_vs_intHES7mersFalse.csv'
mers7 = read_csv(infile)
consolidate(mers7)