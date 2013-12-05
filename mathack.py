#!/usr/bin/env python

import sys
import numpy as np
import scipy.io as sio

datpath = "/Users/josh/Dropbox/binned_heatmaps/"

def chrdat(cell='HIMR',chr=1,corrected=True):
	thepath = datpath + cell + "-all_res40000/" + "chr" + str(chr) + "_chr" + str(chr)
	if corrected:
		thepath += "_corrected"
	thepath += ".mat"
	return thepath

print chrdat(cell="HES",chr=3)


chrmat = sio.loadmat(chrdat(chr=1))
A = chrmat['corrected']

n = A.shape[0]

realrows = [0]*n

for i in range(n):
	if A[i].sum() > 0:
		realrows[i] = 1

window = 50

#in a fixed window around a bin, how much interaction is to left or right
right_spread = np.zeros(n)
left_spread = np.zeros(n)
for i in range(n):
	if not realrows[i]:
		continue
	if i < window or i > (n-window) or :
		continue

	A[i][(i-window):(i+window)].sum()


print "there are "+ str(sum(realrows)) + "real rows out of " + str(n) + " rows"



