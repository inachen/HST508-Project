#!/usr/bin/env python

import sys
import random
import copy
import numpy as np
import scipy.io as sio
import matplotlib
import matplotlib.pyplot as plt
import csv

datpath = "/Users/josh/Dropbox/binned_heatmaps/"

def chrpath(cell='HIMR',chro=1,corrected=True):
	thepath = datpath + cell + "-all_res40000/" + "chr" + str(chro) + "_chr" + str(chro)
	if corrected:
		thepath += "_corrected"
	thepath += ".mat"
	return thepath

def chrmat(cell='HIMR',chro=1,corrected=True):
	themat = sio.loadmat(chrpath(cell, chro, corrected))
	A = themat['corrected']
	return A

#find rows with few reads
def whatisbad(A):
	n = A.shape[0]
	sparsity = np.zeros(n)
	
	if 0:
		for i in range(n):
			count = 0
			for j in range(n):
				if A[i][j] > 0:
					count += 1
			sparsity[i] = count
	
		plt.hist(sparsity)
		plt.title("regions with reads")
		plt.xlabel("value")
		plt.ylabel("frequency")
		plt.show()
	
	goodrows = [1]*n
	for i in range(n):
		if A[i].sum() == 0:
			goodrows[i] = 0
	
	return goodrows



def directionality(A, goodrows=None, window=50):
	n = A.shape[0]
	if goodrows == None:
		goodrows = [1]*n
		
	#in a fixed window around a bin, how much interaction is to left or right
	right_spread = np.zeros(n)
	left_spread = np.zeros(n)
	for i in range(n):
		if not realrows[i]:
			continue
		if i < window or i > (n-window):
			continue
		A[i][(i-window):(i+window)].sum()

#how does contact frequency scale with distance?
def distance_scaling(A,goodrows):
	n = A.shape[0]	
	counts = np.zeros(n)
	for i in range(n):
		counts[i] = np.trace(A,offset=i)
	
	
	goodrows_mat = np.matrix([goodrows])
	g_by_g = goodrows_mat.T*goodrows_mat
	good_pairs = np.zeros(n)
	for i in range(n):
		good_pairs[i] = np.trace(g_by_g,offset=i)
	
	log_probs = np.zeros(n)
	for i in range(n):
		if counts[i] == 0:
			log_probs[i] = -5
		else:
			log_probs[i] = np.log(counts[i]/good_pairs[i])
	
	return log_probs

def smooth_array(a,half_window):
	n = len(a)
	smootha = np.zeros(n)
	w = half_window
	for i in range(n):
		smootha[i] = np.average(a[max(i-w,0):min(i+w+1,n)])

	return smootha

def finding_distance_function():
	log_probs = sio.loadmat("dist_dependence.mat")
	

def distance_effects():
	print "loading chromosomes..."
	Alist = []
	to_try = 23
	ns = []
	log_probs = []
	
	for i in range(to_try):
		print "on chromosome " + str(i+1)
		A = chrmat(chro=(i+1))
		Alist.append(A)
		ns.append(A.shape[0])
		log_probs.append(distance_scaling(Alist[i],whatisbad(Alist[i])))

	sio.savemat("dist_dependence.mat",{"log_probs": log_probs})
	if 0:
		plt.figure(1)
		for i in range(to_try):
			plt.subplot(5,5,i)
			x = np.arange(ns[i])
			y = log_probs[i]
			plt.plot(x,y)
			plt.title("Scaling with distance")
			plt.xlabel("distance")
			plt.ylabel("log average contacts")
			plt.title("Chromosome " + str(i+1))
		plt.show()
		plt.savefig('foo.png')
		plt.close()
	
	#display all on same graph
	n = max(ns)
	plt.figure(2)
	y = np.zeros((n,to_try))
	for i in range(to_try):
		for j in range(ns[i]):
			y[j][i] = log_probs[i][j]
	plt.plot(np.arange(n), y)
	plt.title("Scaling with distance")
	plt.xlabel("distance")
	plt.ylabel("log average contacts")
	plt.title("Chromosomes all")	
	plt.savefig('all.png',bbox_inches='tight')
	plt.show()

	#things get noisy towards the end. what if only 2/3 of the chromosome contributed?

def main():
	#distance_effects()
	finding_distance_function()
main()



