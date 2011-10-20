"""
Various useful functions.
"""
import numpy as np

def parse_ids(arg_str):
	t_ids = arg_str.split(',')
	ids = []
	for s in t_ids:
		if '-' in s:
			id_range = map(int, s.split('-'))
		        for id in range(id_range[0], id_range[1] + 1):
		        	ids.append(id)
		else:
			ids.append(int(s))
	return ids


def getRanks(values, withTies=True):
	"""
	returns ranks (large values w. large ranks)
	"""
	srt_vals = zip(values[:], range(len(values)))
	srt_vals.sort()
	srt_ranks = np.zeros(len(srt_vals))
	if withTies:
		i = 0
		while i < len(srt_vals):
			curVal = srt_vals[i][0]
			j = 1
			while i + j < len(srt_vals) and srt_vals[i + j][0] == curVal:
				j += 1
			rank = i + (j + 1) / 2.0
			max_i = i + j
			while i < max_i:
				srt_ranks[i] = rank
				i += 1

	else:
		srt_ranks = np.arange(len(srt_vals))

	ranks = np.zeros(len(srt_vals))

	for i, (val, val_index) in enumerate(srt_vals):
		ranks[val_index] = srt_ranks[i]
	return ranks


def _kw_get_ranks_(values, withTies=True):
	"""
	returns ranks (large values w. large ranks)
	"""
	srt_vals = zip(values[:], range(len(values)))
	srt_vals.sort()
	srt_ranks = np.zeros(len(srt_vals))
	group_counts = []
	if withTies:
		i = 0
		while i < len(srt_vals):
			curVal = srt_vals[i][0]
			j = 1
			while i + j < len(srt_vals) and srt_vals[i + j][0] == curVal:
				j += 1
			rank = i + (j + 1) / 2.0
			max_i = i + j
			group_counts.append(j)
			while i < max_i:
				srt_ranks[i] = rank
				i += 1

	else:
		srt_ranks = np.arange(len(srt_vals))
		group_counts = [1] * len(srt_vals)

	ranks = np.zeros(len(srt_vals))

	for i, (val, val_index) in enumerate(srt_vals):
		ranks[val_index] = srt_ranks[i]
	return ranks, np.array(group_counts)


def kruskal_wallis(snps, phenVals, useTieCorrection=True, verbose=True):
	from scipy import stats
	if verbose:
		import time
		s1 = time.time()
	assert len(snps[0]) == len(phenVals), "SNPs and phenotypes are not equal length."

	ranks, group_counts = _kw_get_ranks_(phenVals)
	assert len(group_counts) == len(set(ranks)), 'Somethings wrong..'
	tieCorrection = 1.0
	if useTieCorrection:
		n_total = len(ranks)
		ones_array = np.repeat(1.0, len(group_counts))
		s = np.sum(group_counts * (group_counts * group_counts - ones_array))
		s = s / (n_total * (n_total * n_total - 1.0))
		tieCorrection = 1.0 / (1 - s)

	n = len(phenVals)
	ds = np.zeros(len(snps))
	c = 12.0 / (n * (n + 1))
	for i, snp in enumerate(snps):
		ns = np.bincount(snp)
		rs = np.array([0.0, np.sum(snp * ranks)])
		rs[0] = np.sum(ranks) - rs[1]
		nominator = 0.0
		mean_r = (n + 1) / 2.0
		for j in range(2):
			v = (rs[j] / ns[j]) - mean_r
			nominator += ns[j] * v * v
		ds[i] = (c * nominator) * tieCorrection

	ps = stats.chi2.sf(ds, 1)

	#print ps
	if verbose:
		secs = time.time() - s1
		if secs > 60:
			mins = int(secs) / 60
			secs = secs - mins * 60
			print 'Took %d mins and %f seconds.' % (mins, secs)
		else:
			print 'Took %f seconds.' % (secs)
	return {"ps":ps, "ds":ds}



def plotHist(x, y):
	pass


class Queue:
	"A first-in-first-out queue."
	def __init__(self, items=None):
		self.start = 0
		self.A = items or []

	def __len__(self):
		return len(self.A) - self.start

	def append(self, item):
		self.A.append(item)

	def extend(self, items):
		self.A.extend(items)

	def pop(self):
		A = self.A
		item = A[self.start]
		self.start += 1
		if self.start > 100 and self.start > len(A) / 2:
			del A[:self.start]
			self.start = 0
		return item


def valListToStrList(l):
	newl = []
	for val in l:
		newl.append(str(val))
	return newl


#def calcVar(l):
#	mean = sum(l) / float(len(l))
#	var = 0
#	for i in range(0, len(l)):
#		var = var + (l[i] - mean) * (l[i] - mean)
#	var = var / float(len(l) - 1)
#	return var


def transposeDoubleLists(l):
	return map(list, zip(*l))

#def calcSD(l):
#	mean = sum(l) / float(len(l))
#	var = 0
#	for i in range(0, len(l)):
#		var = var + (l[i] - mean) * (l[i] - mean)
#	var = var / float(len(l) - 1)
#	import math
#	return math.sqrt(var)

def calcQuantiles(numbers):
	"""
	Calculates the 1st quantile.
	"""
	import math
	numbers.sort()

	if len(numbers) % 2 == 0:
		#Even
		ns1 = numbers[0:len(numbers) / 2]
		ns2 = numbers[len(numbers) / 2:len(numbers)]
		median = (numbers[len(numbers) / 2 - 1] + numbers[len(numbers) / 2]) / 2.0
	else:
		#Odd
		ns1 = numbers[0:(len(numbers) - 1) / 2]
		ns2 = numbers[(len(numbers) + 1) / 2:len(numbers)]
		median = numbers[(len(numbers) - 1) / 2]

	if len(ns1) % 2 == 0:
		#Even
		q1 = (ns1[len(ns1) / 2 - 1] + ns1[len(ns1) / 2]) / 2.0
		q3 = (ns2[len(ns2) / 2 - 1] + ns2[len(ns2) / 2]) / 2.0
	else:
		q1 = ns1[(len(ns1) - 1) / 2]
		q3 = ns2[(len(ns2) - 1) / 2]

	return (q1, median, q3)

def calcIQR(numbers):
	"""
	Calculates the inter-quantile range.
	"""
	quantiles = calcQuantiles(numbers)
	print quantiles
	return quantiles[2] - quantiles[0]


def bin_counts(values, bins):
	counts = [0] * (len(bins) - 1)
	out_of_bounds_values = []
	for value in values:
		i = 0
		while i + 1 < len(bins) and value > bins[i + 1]:
			i += 1
		if i + 1 < len(bins) and bins[i] < value <= bins[i + 1]:
			counts[i] += 1
		else:
			out_of_bounds_values.append(value)
	if len(out_of_bounds_values):
		print "These values were OOB:", out_of_bounds_values
	return counts


def r_list_2_dict(rList):
	pyDict = {}
	for name, value in zip([i for i in rList.getnames()], [i for i in rList]):
	    	if len(value) == 1:
	    		pyDict[name] = value[0]
	    	else:
	    		pyDict[name] = [i for i in value]
	return pyDict



anti_decoder = {1:0, 0:1}
def anti_binary_snp(snp):
	return [anti_decoder[nt] for nt in snp]

if __name__ == '__main__':
	l = [1, 2, 3, 3, 3, 2, 3, 4, 1, 0, 3, 3, 3]
	ranks = getRanks(l)
	print ranks
	print "Done!"


def snp_kinship_correlation(binary_snps, k, num_perm=0, permutation_snps=None, kinship_type='emma.R'):
	"""
	Performs the Mantel's test between the snps and the kinship matrix
	
	Currently only one kinship matrix is supported
	"""
	#Construct a distance matrix using the snps.	

	import scipy as sp
	import scipy.spatial as spat
	num_snps = len(binary_snps)
	num_lines = len(binary_snps[0])
	if kinship_type == 'emma.R':
		import linear_models as lm
		snps_k = lm.calc_kinship_old(binary_snps)
		snp_dists = []
		for i in range(num_lines):
			#c = i*((num_lines-1)-(i+1)/2)-1
			for j in range(i + 1, num_lines):
				snp_dists.append(snps_k[i, j])
		snp_dists = sp.array(snp_dists)
	elif kinship_type == 'minkowski':
		binary_snps = sp.mat(binary_snps).T
		snp_dists = spat.distance.pdist(binary_snps, 'minkowski') / num_snps
	k_dists = []#sp.ones(num_lines*(num_lines-1)/2)
	for i in range(num_lines):
		#c = i*((num_lines-1)-(i+1)/2)-1
		for j in range(i + 1, num_lines):
			k_dists.append(k[i, j])
	k_dists = sp.array(k_dists)
	corr = sp.corrcoef(snp_dists, k_dists)[1, 0]
	pval = None
	if num_perm:
		import random
		assert permutation_snps != None, 'Currently SNPs are needed for permutation.'
		print "Starting Mantel's permutation test."
		raise NotImplementedError
#		for i in range(num_perm):
#			random.sample()
	return {'corr':corr, 'pval':pval}



def kolmogorov_smirnov_test(self, pvalues):
	"""
	Calculates the kolmogorov smirnov test, with the p-values
	"""
        pass
