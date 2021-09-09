import numpy as np
import matplotlib.pyplot as plt
import scipy as scp
from math import factorial
import pysam

def averageQuality(input_file):
	quality_list = []
	for read in input_file.fetch(until_eof = True):
		quality_list += read.query_qualities.tolist()
		a=read.query_qualities.tolist()
		b=read.query_qualities
		#print a
		#print "xxxxxx"
		#print b
		#break
	return sum(quality_list)/len(quality_list)

def dictionaryBED(bed_file):
	dic_BED={}
	file=open(bed_file)
	for line in file:
		col=line.split("\t")
		dic_BED[col[-1][:-2]] = ":".join(col[:-1])
		print dic_BED[col[-1][:2]]
	return dic_BED

def percentTargetAlignedReads(input_file, dic_BED):
	numAlignedReads=0
	for i in dic_BED.keys():
		numAlignedReads += input_file.count(region = dic_BED[i])
	return float(numAlignedReads)/input_file.mapped

def percentTargetAlignedBases(input_file, dic_BED):
	numAlignedBases = 0
	totalAlignedBases = 0
	for i in dic_BED.keys():
		for pileupcolumn in input_file.pileup(region = dic_BED[i]):
			numAlignedBases += pileupcolumn.n
	for pileupcolumn in input_file.pileup():
		totalAlignedBases += pileupcolumn.n
	return float(numAlignedBases)/totalAlignedBases

def baseCoverage(input_file, dic_BED):
	dic_cov = {}
	count_base_position = 0
	all_target_coverage = []
	total_length = 0
	
	for i in dic_BED.keys():
		dic_cov[i] = []
		for pileupcolumn in input_file.pileup(region = dic_BED[i]):
			dic_cov[i].append(pileupcolumn.n)
			total_length += 1
		all_target_coverage += dic_cov[i]
	mean_coverage_all_target = sum(all_target_coverage)/total_length
	for l in all_target_coverage:
		if l > 0.2*mean_coverage_all_target:
			count_base_position += 1
	uniformityCoverage = float(count_base_position)/total_length
	return mean_coverage_all_target, uniformityCoverage, dic_cov

def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')


#plt.show()
