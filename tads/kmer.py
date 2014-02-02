#!/usr/bin/env python
#
# INSTRUCTIONS FOR USE:
# call program as follows:
#  python kmer.py <FASTA> <BINS>
#  
#
#TODO
#
#INA has new TAD boundary regions
#    they are pretty smeared
#    nonetheless they can be used to pull halfway points between tads
#
#lr_bias
#     use this to call tad boundaries as left boundary or right boundary
#     use that to find midpoints of tads
#
#     build quantitative index. for each kmer, what is its freq correlation with lr_bias?
#
#validate enrichment
#			information theory: what is the distribution of the top mers, histogram in each?
#     if the kmers were poisson distributed (within each family), with mean just depending on average counts in that family, then how many of them would you need to pinpoint a guy against the background. how many would you expect to get right?
#     total counts for the key mers in boundary vs. control
#	    compare to svm on basis of top 100 most frequent
#





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
import heapq

# FASTA datapath
FASTAPATH = "/Users/josh/Documents/HST508-Project/data/hg/"

# BIN locations
BINPATH = "/Users/josh/Dropbox/HST508/tracks/"
BINSIZE = 40000

# folder to write output
OUTPATH = "/Users/josh/Dropbox/HST508/tracks/"

def readBias():
	return readTrack(BINPATH + 'lr_bias.csv')

#returns a dictionary, with a list of bins for each chromosome
def readBins(filename):
	"""reads in a CSV file containing 40kb bins. for now, assuming there is only one type of binlist per file"""
	reader = csv.reader(open(filename, 'rU'))
	chr_num = 0
	bins = {}
	for row in reader:
		if len(row) == 1:
			if row[0][0:3] == 'chr':
				chr_num = int(row[0].lstrip('chr'))
		if len(row) > 1:
				assert chr_num not in bins.keys()
				bins[chr_num] = [int(x) for x in row]
	return bins

def readChrom(chr_num):
	return readSeq(FASTAPATH + 'chr' + str(chr_num) + '.fa')

def readSeq(filename):
	"""reads in a FASTA sequence"""
	seq = []
	with open(filename) as stream:
		for line in stream:
			if line.startswith(">"):
				continue
			seq.append(line.rstrip())

	return "".join(seq)

def writeTrack(track_dict, filename):
	
	writer = csv.writer(open(OUTPATH + filename + '.csv', 'wb'))

	for chr_num in track_dict.keys():
		writer.writerow(['chr'+str(chr_num)])
		writer.writerow(['lr_bias'])
		writer.writerow(track_dict[chr_num])

def readTrack(filename):
	reader = csv.reader(open(filename, 'rU'))
	chr_num = 0
	track = {}
	for row in reader:
		if len(row) == 1:
			if row[0][0:3] == 'chr':
				chr_num = int(row[0].lstrip('chr'))
		if len(row) > 1:
				assert chr_num not in track.keys()
				track[chr_num] = [float(x) for x in row]
	return track
	
#returns indices of the nonzero entries of track
def track_to_bins(track):
	bins = {}
	for chr_num in track.keys():
		bins[chr_num] = [i for i in range(len(track[chr_num])) if track[chr_num][i] > 0]
	return bins
	
#returns bins halfway between ones fed in
def midbins(bins,thresh=50):
	mid_bins = {}
	for chr_num in bins.keys():
		mid_bins[chr_num] = []
		init_bins = bins[chr_num]
		nbins = len(init_bins)
		for i in range(nbins-1):
			if init_bins[i+1] == init_bins[i] + 1:
				continue
			elif init_bins[i+1] > init_bins[i] + thresh:
				continue
			else:
				gap = init_bins[i+1] - init_bins[i]
				#take middle 20 percent
				to_take = -(-gap/5)
				mid = (init_bins[i+1] + init_bins[i])/2
				mid_bins[chr_num].extend(range(mid + (-to_take)/2, mid + to_take/2))
	return mid_bins
		
#how much of the TAD boundaries consist of repeats?
def repeats(bin_file, shift = 0):
	tad_bins = readBins(BINPATH+bin_file)
	bin_size = BINSIZE
	upper_ct = 0
	lower_ct = 0
	N_ct = 0
	for i in range(1, 23):
		print "starting chromosome " + str(i)
		seq = readChrom(i)
		bins = tad_bins[i]
		for b in bins:
			seq_start = bin_size*(b + shift)
			seq_end = bin_size + seq_start
			if seq_end > len(seq):
				continue
			for i in range(seq_start,seq_end):
				if seq[i] == 'N':
					N_ct += 1
				elif seq[i].isupper():
					upper_ct += 1
				elif seq[i].islower():
					lower_ct += 1
				else:
					print "junk sequence"
					sys.exit(1)
	print "In all "+ bin_file.rstrip('.csv') + ", there are %d N's, %d uppers, and %d lowers." % (N_ct, upper_ct, lower_ct)
		
#takes in a list of which chromosomes to check, and a list of bin dictionaries bins_dicts[j][i] corresponds to family j on chrom i
#returns a dictionary with keys kmers and values tuples of counts for each family
def kmers_all(which_chr, bins_dicts, bin_size, k, keep_repeats=True):
	n = len(bins_dicts)
	kmer_dict = {}
	
	for i in which_chr:
		seq = readChrom(i)
		for j in range(n):
			bins = bins_dicts[j][i]
			kd = kmers_in_bins(seq, bins, bin_size, k, shift=0, keep_repeats=keep_repeats)
			for mer, v in kd.iteritems():
				if mer not in kmer_dict:
					kmer_dict[mer] = [0]*n
				kmer_dict[mer][j] += v
	return kmer_dict

def kmers_in_bins(seq, bins, bin_size,k,shift=0,keep_repeats=True):
	kmer_dict = {}
	seq_len = len(seq)
	
	#ignore repeat flags
	if keep_repeats:
		seq = seq.upper()
	
	print "counting %d-mers in %d bins of width %d" % (k, len(bins), bin_size)
	for b in bins:
		seq_start = bin_size*(b+shift)
		seq_end = bin_size + seq_start - k + 1
		
		if seq_end > seq_len or seq_start < 0:
			continue
		for i in range(seq_start,seq_end):
			if i == seq_start:
				key = seq[i:i+k]
			#strip first character add last
			else:
				key = key[1:] + seq[i+k]
			if not key.isupper():
				continue
			if key.find('N') > -1:
				continue
			if key in kmer_dict:
				kmer_dict[key] += 1
			else:
				kmer_dict[key] = 1
	return kmer_dict
	

def kmer_test():
	kd = kmers_in_bins(readChrom(2), readBins(BINPATH + 'badbins.csv')[2],40000,2,shift=0,keep_repeats=False)
	print kd

#enrichment in a quantitative track
#computed as correlation between mer-count and track value
#note that different mers are massively correlated, both because of sequence overlap
#and because of gc-content
#NOTDONE!
def q_track_enrichment(which_chr, track, k, keep_repeats=True):
	#need mean and variance for each track, together with their inner product
	mers_seen = {}
	mer_q_dict = {}
	mer_mer_dict = {}
	mer_dict = {}
	
	q_sum = 0
	q_q_sum = 0
	n_bins = 0
	
	badbins = readBins(BINPATH + "badbins.csv")
	for chr_num in which_chr:
		print "on chromosome " + str(chr_num)
		seq = readChrom(chr_num)
		if keep_repeats:
			seq = seq.toupper()
		seq_len = len(seq)
		nbins = len(seq)/BINSIZE + 1
		bins = range(nbins)
		for b in bins:
			if b in badbins[chr_num]:
				continue
			
			q = track[chr_num][b]
			q_sum += q
			q_q += q*q
			n_bins += 1
			
			seq_start = BINSIZE*b
			seq_end = BINSIZE + seq_start - k + 1
			for i in range(seq_start,seq_end):
				if i == seq_start:
					key = seq[i:i+k]
				else:
					key = key[1:] + seq[i+k]
				if not key.isupper():
					continue
				if key.find('N') > -1:
						continue
				if key not in mers_seen:
					mers_seen[key] = 1
				if key in mers_seen:
					
					
					kmer_dict[key] += 1
				else:
					kmer_dict[key] = 1
			
		track = [0]*nbins
		
	for b in range(nbins):
		print "on bin " + str(b)
		track[b] = count_mers(top_tad_mers, k, b, seq)
	

#get the bins that geoff uses
def geoff_bins():

	return 1

#compare tad boundaries, in Ina's new dataset, versus tad interiors
def tad_bd_vs_int(k=7, keep_repeats=False):
	bou_bins = track_to_bins(readTrack(BINPATH+'LdR-tads-HES/LRbounds-output.csv'))
	int_bins = midbins(bou_bins)
	
	#control
	con_bins = {}
	bas_bins = {}
	
	for chr_num in bou_bins.keys():
		con_bins[chr_num] = [x + 25 for x in bou_bins[chr_num]]
	bas_bins = midbins(con_bins)
	
	bad_bins = readBins(BINPATH + "badbins.csv")
	
	#clean out bad bins from each
	bou_bins = remove_overlap(bad_bins, bou_bins)
	int_bins = remove_overlap(bad_bins, int_bins)
	con_bins = remove_overlap(bad_bins, con_bins)
	bas_bins = remove_overlap(bad_bins, bas_bins)
	
	kd = kmers_all(range(2,23), [bou_bins, int_bins, con_bins, bas_bins], BINSIZE, k,keep_repeats=keep_repeats)

	top_bou = heapq.nlargest(5, kd.iteritems(), key=lambda x: x[1][0])
	top_int = heapq.nlargest(5, kd.iteritems(), key=lambda x: x[1][1])
	print "top five mers for tad boundaries"
	print top_bou
	print "top five mers for tad interiors"
	print top_int
	
	#top hits, corrected
	top_bou = heapq.nlargest(20, kd.iteritems(), key=lambda x: float(x[1][0])/float(x[1][2]+x[1][3]+1))
	top_int = heapq.nlargest(20, kd.iteritems(), key=lambda x: float(x[1][1])/float(x[1][2]+x[1][3]+1))
	print "top five mers for tad boundaries rel baseline"
	print top_bou
	print "top five mers for tad interiors rel baseline"
	print top_int
	
	writer = csv.writer(open(OUTPATH + 'bd_vs_intHES' + str(k) +  'mers' + str(keep_repeats) + '.csv', 'wb'))
	for k, v in kd.iteritems():
		writer.writerow([k,v[0],v[1],v[2]])
	
def remove_overlap(scrub, new):
	for chr_num in new.keys():
		for b in new[chr_num]:
			if b in scrub[chr_num]:
				new[chr_num].remove(b)
	return new
	
#read in boundaries, and find frequency of each kmer in tads, and two control shifts 
def tads_vs_shifts(k=7, shift=13,keep_repeats=True):
	tad_bins = readBins(BINPATH + 'lr_bounds.csv')
	tad_bins_R = {}
	tad_bins_L = {}
	for chr_num in tad_bins.keys():
		tad_bins_R[chr_num] = [x + shift for x in tad_bins[chr_num]]
		tad_bins_L[chr_num] = [x - shift for x in tad_bins[chr_num]]
	
	#remove any which happen to be tad_bins in the first place, and which are bad
	bad_bins = readBins(BINPATH + "badbins.csv")
	for chr_num in tad_bins.keys():
		for b in tad_bins_R[chr_num]:
			if b in tad_bins[chr_num] or b in bad_bins[chr_num]:
				tad_bins_R[chr_num].remove(b)
		for b in tad_bins_L[chr_num]:
			if b in tad_bins[chr_num] or b in bad_bins[chr_num]:
				tad_bins_L[chr_num].remove(b)
	
	kd = kmers_all(range(1,23), [tad_bins, tad_bins_R, tad_bins_L], BINSIZE, k,keep_repeats=keep_repeats)
	
	#what are the top hits, uncorrected?

	top_ten_tad = heapq.nlargest(5, kd.iteritems(), key=lambda x: x[1][0])
	top_ten_tad_R = heapq.nlargest(5, kd.iteritems(), key=lambda x: x[1][1])
	print "top five mers for tads"
	print top_ten_tad
	print "top five mers for control"
	print top_ten_tad_R
	
	#top hits, corrected
	top_ten_tad = heapq.nlargest(20, kd.iteritems(), key=lambda x: float(x[1][0])/float(x[1][2]+1))
	top_ten_tad_R = heapq.nlargest(20, kd.iteritems(), key=lambda x: float(x[1][1])/float(x[1][2]+1))
	print "top five mers for tads rel baseline"
	print top_ten_tad
	print "top five mers for control rel baseline"
	print top_ten_tad_R

	writer = csv.writer(open(OUTPATH + 'tad' + str(k) +  'mers' + str(keep_repeats) + '.csv', 'wb'))
	for k, v in kd.iteritems():
		writer.writerow([k,v[0],v[1],v[2]])


#can the enrichments found by comparing tad boundaries to neighbors be used to distinguish them from such?
def validate_enrichment2(infile, k=7):
	Nmers = 100
	shift=13
	
	#get top 7mers
	#this is a hack, using the old read_kmers
	(kmers, tad_enrichments, control_enrichments) = read_kmers(infile)
	top_bou_mers_and_counts = heapq.nlargest(Nmers, kmers, key=lambda x: float(x[1])/float(x[2]+1))
	top_int_mers_and_counts = heapq.nlargest(Nmers, kmers, key=lambda x: float(x[2])/float(x[3]+1))
	top_tad_mers = frozenset([t[0] for t in top_bou_mers_and_counts])
	top_int_mers = frozenset([t[0] for t in top_int_mers_and_counts])
	
	bou_bins = track_to_bins(readTrack(BINPATH+'LdR-tads-HES/LRbounds-output.csv'))
	int_bins = midbins(bou_bins)
	
	#control
	con_bins = {}
	bas_bins = {}
	
	for chr_num in bou_bins.keys():
		con_bins[chr_num] = [x + 25 for x in bou_bins[chr_num]]
	bas_bins = midbins(con_bins)
	
	bad_bins = readBins(BINPATH + "badbins.csv")
	
	#clean out bad bins from each
	bou_bins = remove_overlap(bad_bins, bou_bins)
	int_bins = remove_overlap(bad_bins, int_bins)
	con_bins = remove_overlap(bad_bins, con_bins)
	bas_bins = remove_overlap(bad_bins, bas_bins)
	
	tad_bins = bou_bins
	control_bins = int_bins
	
	#iterate through the genome, counting the number of top_mers for tads, subtracting those for int
	which_chr = range(2,23)
	tad_counts = []
	control_counts = []
	for i in which_chr:
		print "on chromosome " + str(i)
		seq = readChrom(i)
		for b in tad_bins[i]:
			#tad_counts.append(count_mers(top_tad_mers, k, b, seq) - count_mers(top_int_mers, k, b, seq))
			tad_counts.append(count_mers(top_tad_mers, k, b, seq))
		for b in control_bins[i]:
			#control_counts.append(count_mers(top_tad_mers, k, b, seq) - count_mers(top_int_mers, k, b, seq))
			control_counts.append(count_mers(top_tad_mers, k, b, seq))
		
	print "average count in tad boundaries is: " + str(float(sum(tad_counts))/len(tad_counts))	
	print "average count in control is: " + str(float(sum(control_counts))/len(control_counts))
			
	bmax = max([max(tad_counts), max(control_counts)])
	bmin = min([min(tad_counts), min(control_counts)])
	plt.hist(tad_counts, bins=list(np.linspace(bmin,bmax,40)), normed=1, histtype='step',color='blue',label='tad')
	plt.grid(True)
	
	plt.hist(control_counts, bins=np.linspace(bmin,bmax,40), normed=1, histtype='step',color='red',label='control')
	plt.legend()
	plt.title('Histogram of counts of TADmers')
	plt.show()
	
	writer = csv.writer(open(OUTPATH + 'tad_enrich_valid_Ina7_HES2'+ str(Nmers) + '.csv', 'wb'))
	writer.writerow(['tad_counts'])
	writer.writerow(tad_counts)
	writer.writerow(['control_counts'])
	writer.writerow(control_counts)
	
	draw_ROC(tad_counts,control_counts)
	
def validate_enrichment(infile):
	Nmers = 100
	k = 7
	
	#get top 7mers
	(kmers, tad_enrichments, control_enrichments) = read_kmers(infile)
	top_tad_mers_and_counts = heapq.nlargest(Nmers, kmers, key=lambda x: float(x[1])/float(x[3]+1))
	top_tad_mers = frozenset([t[0] for t in top_tad_mers_and_counts])
	
	#get bins, and comparisons
	tad_bins = readBins(BINPATH + 'lr_bounds.csv')
	tad_bins_R = {}
	tad_bins_L = {}
	for chr_num in tad_bins.keys():
		tad_bins_R[chr_num] = [x + shift for x in tad_bins[chr_num]]
		tad_bins_L[chr_num] = [x - shift for x in tad_bins[chr_num]]
	
	#remove any which happen to be tad_bins in the first place, and which are bad
	bad_bins = readBins(BINPATH + "badbins.csv")
	for chr_num in tad_bins.keys():
		for b in tad_bins_R[chr_num]:
			if b in tad_bins[chr_num] or b in bad_bins[chr_num]:
				tad_bins_R[chr_num].remove(b)
		for b in tad_bins_L[chr_num]:
			if b in tad_bins[chr_num] or b in bad_bins[chr_num]:
				tad_bins_L[chr_num].remove(b)
	
	control_bins = tad_bins_R
	
	#iterate through the genome, counting the number of top_mers
	which_chr = range(1,23)
	tad_counts = []
	control_counts = []
	for i in which_chr:
		print "on chromosome " + str(i)
		seq = readChrom(i)
		for b in tad_bins[i]:
			tad_counts.append(count_mers(top_tad_mers, k, b, seq))
		for b in control_bins[i]:
			control_counts.append(count_mers(top_tad_mers, k, b, seq))
		
	print "average count in tad boundaries is: " + str(float(sum(tad_counts))/len(tad_counts))	
	print "average count in control is: " + str(float(sum(control_counts))/len(control_counts))
			
	bmax = max([max(tad_counts), max(control_counts)])
	bmin = min([min(tad_counts), min(control_counts)])
	plt.hist(tad_counts, bins=list(np.linspace(bmin,bmax,40)), normed=1, histtype='step',color='blue',label='tad')
	plt.grid(True)
	
	plt.hist(control_counts, bins=np.linspace(bmin,bmax,40), normed=1, histtype='step',color='red',label='control')
	plt.legend()
	plt.title('Histogram of counts of TADmers')
	plt.show()
	
	writer = csv.writer(open(OUTPATH + 'tad_enrich_valid_Ina7_HES'+ str(Nmers) + '.csv', 'wb'))
	writer.writerow(['tad_counts'])
	writer.writerow(tad_counts)
	writer.writerow(['control_counts'])
	writer.writerow(control_counts)
	
	draw_ROC(tad_counts,control_counts)

	
def draw_ROC(dat1, dat2):
	dmin = min(min(dat1),min(dat2))
	dmax = max(max(dat1),max(dat2))
	bins = 200
	tot1 = len(dat1)
	tot2 = len(dat2)
	
	y = []
	x = []
	for i in np.linspace(dmin,dmax,bins):
		y.append(len([d for d in dat1 if d > i])/float(tot1)*100)
		x.append(len([d for d in dat2 if d > i])/float(tot2)*100)
	
	plt.plot(x,y,np.linspace(0,100,bins),np.linspace(0,100,bins))
	plt.show()
	
	diff = [y[i] - x[i] for i in range(bins)]
	i = np.argmax(diff)
	
	print "max of true positive - false positive is %f - %f = %f" % (y[i],x[i],y[i]-x[i]) 
	
def count_mers(kmers, k, bin, seq):
	count = 0
	seq_start = BINSIZE*bin
	seq_end = BINSIZE + seq_start - k + 1
	if seq_end > len(seq):
		return 0
		
	#loop over all k-mers, skipping repeats
	for i in range(seq_start,seq_end):
		if i == seq_start:
			key = seq[i:i+k]
		else:
			key = key[1:] + seq[i+k]
		if key in kmers:
			count += 1

	return count

#read a file of kmer enrichments. take the top 100 of them kmer
#for 7mers, there appears to be much headroom for enrichments above 1.25
#return a track with the counts in each region
def top_kmers_track():
	infile = OUTPATH + "tad7mersFalse.csv"
	Nmers = 1000
	k = 7
	
	(kmers, tad_enrichments, control_enrichments) = read_kmers(infile)

	top_tad_mers_and_counts = heapq.nlargest(Nmers, kmers, key=lambda x: float(x[1])/float(x[3]+1))
	top_tad_mers = frozenset([t[0] for t in top_tad_mers_and_counts])
	
	seq = readChrom(1)
	seq_len = len(seq)
	nbins = len(seq)/BINSIZE + 1
	track = [0]*nbins
	for b in range(nbins):
		print "on bin " + str(b)
		track[b] = count_mers(top_tad_mers, k, b, seq)
	
	top_count = max(track)
	tad_bins = readBins(BINPATH + 'lr_bounds.csv')
	tad_track = [0]*nbins
	for i in tad_bins[1]:
		tad_track[i] = top_count/2
	plt.plot(np.arange(len(track)),track,np.arange(len(track)),tad_track)
	plt.show()
	
	writeTrack({1:tad_track},'1000_7mers_enriched_in_tad_boundaries')

def read_kmers(infile):
	reader = csv.reader(open(infile, 'rU'))
	
	kmers = []
	tad_enrichments = []
	control_enrichments = []
	totals = [0,0,0]
	for row in reader:
		mer = row[0]
		tad = int(row[1])
		control = int(row[2])
		baseline = int(row[3])
		tad_enrichment = float(tad)/(1+float(baseline))
		control_enrichment = float(control)/(1+float(baseline))
		#ignore N's for now
		if mer.find('N') > -1:
			continue
		
		totals[0] += tad
		totals[1] += control
		totals[2] += baseline
		kmers.append([mer, tad, control, baseline])
		tad_enrichments.append(tad_enrichment)
		control_enrichments.append(control_enrichment)
	
	totals = [float(x) for x in totals]
	tad_enrichments = [x * totals[2]/totals[0] for x in tad_enrichments]
	control_enrichments = [x * totals[2]/totals[1] for x in control_enrichments]

	return (kmers, tad_enrichments, control_enrichments)

def kmer_analysis(infile):
	reader = csv.reader(open(infile, 'rU'))
	
	kmers = []
	tad_enrichments = []
	control_enrichments = []
	totals = [0,0,0]
	for row in reader:
		mer = row[0]
		tad = int(row[1])
		control = int(row[2])
		baseline = int(row[3])
		tad_enrichment = float(tad)/(1+float(baseline))
		control_enrichment = float(control)/(1+float(baseline))
		#ignore N's for now
		if mer.find('N') > -1:
			continue
		
		totals[0] += tad
		totals[1] += control
		totals[2] += baseline
		kmers.append([mer, tad, control, baseline])
		tad_enrichments.append(tad_enrichment)
		control_enrichments.append(control_enrichment)
	
	totals = [float(x) for x in totals]
	tad_enrichments = [x * totals[2]/totals[0] for x in tad_enrichments]
	control_enrichments = [x * totals[2]/totals[1] for x in control_enrichments]

	#print len(tad_enrichments)
	#print len(control_enrichments)

	plt.hist(tad_enrichments, bins=list(np.linspace(0,2,50)), histtype='step',color='blue',label='tad')
	plt.grid(True)
	
	plt.hist(control_enrichments, bins=np.linspace(0,2,50), histtype='step',color='red',label='control')
	plt.legend()
	plt.title('Histogram of enrichments')
	plt.show()

	print len([x for x in tad_enrichments if x > 1.25])

	plt.hist([x[1]  for x in kmers], bins=np.linspace(0,100000,300))
	plt.show()
	
def main():
	
	# NOTE to users:
  #   If you do not want to use the command line, comment out the command line
  #   parsing code and hard-code the input filenames.
  #
  # For example, use the following:
  # file1 = "human-hoxa-region.fa"
  # file2 = "mouse-hoxa-region.fa"
  # plotfile = "dotplot.jpg"

  # parse command-line arguments
	if len(sys.argv) == 2:
		print "you must call program as:  "
		print "   python kmer.py <FASTA> <BINS>"
		print "   PLOT FILE may be *.ps, *.png, *.jpg"
		sys.exit(1)
	fasta_file = sys.argv[1]
	bin_file = sys.argv[2]


	file1 = DATAPATH
	

	print "reading sequence"
	seq1 = readSeq(file1)

	# length of hash key
	kmerlen = 7

	print "hashing seq1..."
	for i in xrange(len(seq1) - kmerlen + 1):
		#key = seq1[i:i+kmerlen]
		key = ''.join([seq1[i] for i in range(i,i+kmerlen) if i % skip != 0])
		#key = list_mod(key,skip)
		lookup.setdefault(key, []).append(i)
	
#repeats('lr_bounds.csv',0)	
#kmer_test()
#readBins(BINPATH + "badbins.csv")
#tads_vs_shifts(k=9,shift=13,keep_repeats=False)
#tads_vs_shifts(shift=13,keep_repeats=False)
#kmer_analysis(OUTPATH + "tad7mersFalse.csv")
#top_kmers_track()
#validate_enrichment(OUTPATH + "tad7mersFalse.csv")
validate_enrichment2(OUTPATH + "bd_vs_intHES7mersFalse.csv")
#tad_bd_vs_int()
#kmer_analysis(OUTPATH + "bd_vs_intHES7mersFalse.csv")