#!/usr/bin/env python

import sys
import random
import copy
import numpy as np
import scipy.io as sio
import scipy as sp
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.image as img
import csv

# Hi-C data folder
DATAPATH = "/Users/josh/Dropbox/binned_heatmaps/"

# folder to write output
OUTPATH = "/Users/josh/Dropbox/HST508/tracks/"

# number of chromosomes in human cells
NUM_CHRMS = 23

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

#find rows with few reads
def good_rows(A):
	n = A.shape[0]	
	goodrows = [1]*n
	for i in range(n):
		if A[i].sum() == 0:
			goodrows[i] = 0	
	return goodrows

#strips columns and rows which are entirely zero
def clean_mat(A):
	n = A.shape[0]	
	cols_to_keep = []
	for i in range(n):
		if A[i].sum() > 0:
			cols_to_keep += [i]
	return (A[np.ix_(cols_to_keep,cols_to_keep)], cols_to_keep)

#returns log of average bin value at each distance d = i-j
def raw_distance_scaling(A,goodrows):
	n = A.shape[0]	
	counts = np.zeros(n)
	for i in range(n):
		counts[i] = np.trace(A,offset=i)

	goodrows_mat = np.matrix([goodrows])
	g_by_g = goodrows_mat.T*goodrows_mat
	good_pairs = np.zeros(n)
	for i in range(n):
		good_pairs[i] = np.trace(g_by_g,offset=i)
	
	avg_counts = np.zeros(n)
	for i in range(n):
		if counts[i] > 0:
			avg_counts[i] = counts[i]/float(good_pairs[i])
	
	return avg_counts

#the raw scaling has no columns in the first two, and gets ugly at the end.
#so take the average from 2/3 way away to 3/4, and extend to end
def smooth_distance_scaling(A, goodrows):
	n = A.shape[0]
	counts = raw_distance_scaling(A, goodrows)
	
	n1 = (2*n)/3
	n2 = 3*n/4
	print n1
	print n2
	m = np.mean(counts[n1:n2])
	for i in range(n2,n):
		counts[i] = m
	return counts

#normalizes entries so each distance d = |i-j| has the same average number of counts
def distance_normalize(A_in, stash=False, stashspot='nowhere'):
	A = copy.copy(A_in)
	n = A.shape[0]
	goodrows = good_rows(A)
	#rescale
	dist_scale = smooth_distance_scaling(A,goodrows)

	for i in range(n):
		for j in range(n):
			if A[i][j] > 0 and dist_scale[abs(i-j)] > 0:
				A[i][j] = A[i][j]/dist_scale[abs(i-j)]	

	return (A, dist_scale)

def distance_de_normalize(A, dist_scale):

	n = A.shape[0]

	for i in range(n):
		for j in range(n):
			if A[i][j] > 0 and dist_scale[abs(i-j)] > 0:
				A[i][j] = A[i][j]*dist_scale[abs(i-j)]	

	return A


def test_ddn():
	A = load_chr(chr=17)
	#A = np.random.random((10,10))
	B, dist_scale = distance_normalize(A)
	print abs(A - distance_de_normalize(B,dist_scale)).max()

#truncates the top few hits
def trunc(A, high=0.0005):
	lim = np.percentile(A, 100. * (1 - high))
	print "dataset truncated at %lf" % (lim)
	A[A > lim] = lim
	return A

def test_trunc():
	A = np.random.random((7,7))
	print A
	print "truncated the top 30 percent:"
	print trunc(A, high=.3)

def log_for_show(A):
	min_nonzero = min(A[A > 0])
	B = copy.copy(A)
	B[B == 0] = min_nonzero
	return np.log(B)

def test_log_for_show():
	A = np.random.random((5,5))
	A[0][0] = 0
	A[3][3] = 0
	print A
	print np.log(A)
	print log_for_show(A)

def eigen_signal(A):
	#truncate, normalize, truncate
	n = A.shape[0]
	A = trunc(A)
	A, d_scale = distance_normalize(A)
	A = trunc(A)
	#strip zero rows
	B, cols_kept = clean_mat(copy.copy(A))
	#normalize to zero mean
	m = B.shape[0]
	mu = B.mean()
	B = B - np.ones(m)*mu

	#smooth out the matrix by averaging locally
	C = smooth_array(B,2)

	#eigenvector decomposition
	(V, W) = np.linalg.eig(C)
	eig_vecs = W[:,:3]
	eig_vals = V[:3]
	
	#restore missing rows
	full_eig_vecs = np.zeros((n,3))
	full_eig_vecs[np.ix_(cols_kept,range(3))] = eig_vecs
	
	return eig_vals, full_eig_vecs

def test_eigen_signal():
	eig_vals, eig_vecs = eigen_signal(load_chr(chr=17))
	
	x_vals = np.arange(eig_vecs.shape[0])
	plt.plot(x_vals, eig_vecs[:,0], x_vals, eig_vecs[:,1], x_vals, eig_vecs[:,2])
	plt.show()
	
def export_eigen_signals():
	for i in range(NUM_CHRMS):
		export_dict = {}
		eig_vals, eig_vecs = eigen_signal(load_chr(chr=i+1))
		export_dict["evecs"] = eig_vecs
		export_dict["evals"] = eig_vals	
		sio.savemat(OUTPATH + "compartment_eigs_chr" + str(i+1) + ".mat", export_dict)

def eigen_compartments(A):
	plt.subplot(2,2,1)
	plt.imshow(log_for_show(A))
	
	A = trunc(A)
	
	plt.subplot(2,2,2)
	plt.imshow(log_for_show(A))
	
	raw_A = A
	
	A, d_scale = distance_normalize(A)
	A = trunc(copy.copy(A))
	
	plt.subplot(2,2,3)
	plt.imshow(log_for_show(A))
	
	#strip zero rows
	B, cols_kept = clean_mat(copy.copy(A))
	m = B.shape[0]

	#normalize to zero mean
	mu = B.mean()
	B = B - np.ones(m)*mu
	plt.subplot(2,2,4)
	plt.imshow(B)
	plt.show()
	
	#eigenvector decomposition
	#jesus it wants to do all of them. i guess i should just dump the matrix into matlab at that point.

	#smooth out the matrix by averaging locally
	C = smooth_array(B,2)
	m = C.shape[0]

	print "about to do linear algebra"
	(V, W) = np.linalg.eig(C)
	eig_vec = W[:,np.argmax(V)]
	eig_val = V[np.argmax(V)]
	
	top_3_mat = np.zeros((m,m))
	for i in range(3):
		eig_vec_mat = np.matrix([W[:,i]])
		top_3_mat += V[i]*eig_vec_mat.T*eig_vec_mat
	
	
	plt.subplot(2,2,1)
	plt.imshow(B)
	plt.subplot(2,2,2)
	plt.imshow(B-top_3_mat)
	
	#first restore mean
	top_3_mat += np.ones(m)*mu
	#restore missing columns
	T3 = np.zeros((raw_A.shape[0],raw_A.shape[0]))
	T3[np.ix_(cols_kept,cols_kept)] = top_3_mat
	A_eig_and_dist = distance_de_normalize(T3, d_scale)
	
	plt.subplot(2,2,3)
	plt.imshow(log_for_show(raw_A))	
	plt.subplot(2,2,4)
	plt.imshow(log_for_show(abs(raw_A - A_eig_and_dist)))
	
	plt.show()
	#eig_vec_mat = np.matrix([eig_vec])
	#eig_vec_mat = eig_val*eig_vec_mat.T*eig_vec_mat
	#plt.subplot(1,1,1)
	#plt.imshow(eig_vec_mat)
	#plt.show()
	plt.subplot(1,1,1)
	plt.plot(np.arange(m),eig_vec)
	plt.show()
	import pdb; pdb.set_trace()

def test_smooth_array():
	a = np.arange(25).reshape((5,5))
	print a
	print smooth_array(a,1)
	
def smooth_array(a,half_window):
	w = half_window
	if len(a.shape) == 1:
		print "hullo"
		n = len(a)
		smootha = np.zeros(n)
		for i in range(n):
			smootha[i] = np.average(a[max(i-w,0):min(i+w+1,n)])

	#use displacement, and fix the corners	
	if len(a.shape) == 2:
		n = a.shape[0]
		smootha = np.zeros((n + 2*w, n + 2*w))
		for i in range(2*w+1):
			for j in range(2*w+1):
				smootha += np.pad(a, ((i,2*w-i),(j,2*w-j)),mode='edge')
		smootha = smootha[w:n+w,w:n+w]
		smootha = smootha / float(4*w*w + 4*w + 1)
		
		'''smootha = copy.copy(a)[2*w:,2*w:]
		for i in range(2*w):
			smootha += a[i:-2*w+i,i:-2*w+i]
		smootha = smootha / float(2*w+1)
		'''
	return smootha

def main():
	#test_trunc()
	#test_log_for_show()
	eigen_compartments(load_chr(chr=17))
	#export_eigen_signals()
	#test_ddn()
	#test_smooth_array()
	#distance_effects()
	#finding_distance_function()

main()



