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
OUTPATH = "output/"

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
def distance_normalize(A):
	n = A.shape[0]
	counts = np.zeros(A.shape[0])
	for i in range(A.shape[0]):
		counts[i] = np.average(np.diag(A,k=i))

	n = A.shape[0]
	goodrows = good_rows(A)
	#rescale
	dist_scale = smooth_distance_scaling(A,goodrows)

	for i in range(n):
		for j in range(n):
			if A[i][j] > 0 and dist_scale[abs(i-j)] > 0:
				A[i][j] = A[i][j]/dist_scale[abs(i-j)]	
	
	n = A.shape[0]
	counts = np.zeros(A.shape[0])
	for i in range(A.shape[0]):
		counts[i] = np.average(np.diag(A,k=i))

	return A

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

def eigen_compartments(A):
	plt.subplot(2,2,1)
	plt.imshow(log_for_show(A))
	
	A = trunc(A)
	
	plt.subplot(2,2,2)
	plt.imshow(log_for_show(A))
	
	A = distance_normalize(A)
	A = trunc(A)
	
	plt.subplot(2,2,3)
	plt.imshow(log_for_show(A))
	
	#strip zero rows
	(B, cols_kept) = clean_mat(A)
	m = B.shape[0]

	#normalize to zero mean
	B = B - np.ones(m)*B.mean()
	plt.subplot(2,2,4)
	plt.imshow(B)
	plt.show()
	
	#eigenvector decomposition
	#jesus it wants to do all of them. i guess i should just dump the matrix into matlab at that point.

	#before doing eigenvectors/values, should we do an average?
	

	print "about to do linear algebra"
	(V, W) = np.linalg.eig(B)
	eig_vec = W[np.argmax(V)]
	eig_val = V[np.argmax(V)]
	
	eig_val_mat = np.matrix([eig_val])
	eig_vec_mat = eig_val*eig_val_mat.T*eig_val_mat
	plt.subplot(2,2,2)
	plt.imshow(eig_vec_mat)
	plt.show()
	plt.subplot(1,1,1)
	plt.plot(np.arange(m),eig_vec)
	plt.show()
	import pdb; pdb.set_trace()
	
	



def smooth_array(a,half_window):
	n = len(a)
	smootha = np.zeros(n)
	w = half_window
	for i in range(n):
		smootha[i] = np.average(a[max(i-w,0):min(i+w+1,n)])
	return smootha

def finding_distance_function():
	#something funky happened with storage
	qp = sio.loadmat("dist_dependence.mat")["log_probs"][0]
	lp = []
	m = len(qp)
	for i in range(m):
		lp.append(qp[i][0])
	ns = []
	for i in range(m):
		ns.append(len(lp[i]))

	#for global average, just count the first 2/3 of chromosome length
	n_max = max(ns)
	n_max = n_max*2/3
	avg_lp = np.zeros(n_max)
	for i in range(n_max):
		k = 0
		s = 0.0
		for j in range(m):
			if i <= ns[m]*2/3:
				k += 1
				s += lp[j][i]
			avg_lp = s/float(k)

	plt.plot(np.arange(avg_))


#compute diagonal sum

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

	lp_dict = {}
	for i in range(to_try):
		lp_dict["dist_dep_chr" + str(i+1)] = log_probs[i]
	sio.savemat("dist_deps.mat",lp_dict)

	#sio.savemat("dist_dependence.mat",{"log_probs": log_probs})
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
	#test_trunc()
	#test_log_for_show()
	eigen_compartments(load_chr(chr=17))
	#distance_effects()
	#finding_distance_function()

main()



