"""
This python library aims to do two things.
1. Offer general wrapper classes around SNPs datasets.
2. Offer basic functions which can aid analysis of the SNPs.

Bjarni Vilhjalmsson, bvilhjal@usc.edu
"""

import sys, warnings
import pdb
import env
from itertools import *
from bisect import bisect

try:
	import scipy as sp
except Exception, err_str:
	print 'scipy is missing:', err_str


IUPAC_alphabet = ['A', 'C', 'G', 'T', '-', 'Y', 'R', 'W', 'S', 'K', 'M', 'D', 'H', 'V', 'B', 'X', 'N']

"""
Marker type suggestions:
 'Unknown'
 'SNP'
 'Indel'
 'Other'
  ... 
  etc.
  
"""

def get_haplotypes(snps, num_accessions, count_haplotypes=False):
	curr_haplotypes = map(tuple, sp.transpose(sp.array(snps).tolist()))
	hap_set = list(set(curr_haplotypes))
	hap_hash = {}
	for j, h in enumerate(hap_set):
		hap_hash[h] = j
	new_haplotypes = sp.zeros(num_accessions, dtype='int8')
	for j, h in enumerate(curr_haplotypes):
		new_haplotypes[j] = hap_hash[h]
	if count_haplotypes:
		hap_counts = [curr_haplotypes.count(hap) for hap in hap_set]
		return hap_set, hap_counts, new_haplotypes
	else:
		return new_haplotypes



#def coordinateSnpsAndPhenotypeData(phed, p_i, sd, onlyBinarySNPs=True, data_format='binary'):
#	"""
#	1. Remove accessions which are not represented in either of the two datasets
#	2. Order the data in same way.
#	3. Remove monomorphic SNPs
#	"""
#	import bisect
#	print "Coordinating SNP and Phenotype data."
#	ets = phed.phen_dict[p_i]['ecotypes']
#	#Checking which accessions to keep and which to remove.
#	common_ets = list(set(sd.accessions).union(set(ets)))
#	common_ets.sort()
#
#	sd_indices_to_keep = []
#	for i, acc in enumerate(sd.accessions):
#		b_i = bisect.bisect_left(common_ets, acc)
#		if b_i < len(common_ets) and common_ets[b_i] == acc:
#			sd_indices_to_keep.append(i)
#	pd_indices_to_keep = []
#	for i, acc in enumerate(ets):
#		b_i = bisect.bisect_left(common_ets, acc)
#		if b_i < len(common_ets) and common_ets[b_i] == acc:
#			pd_indices_to_keep.append(i)
#
#
#	#Filter accessions which do not have the phenotype value (from the genotype data).
#	self.filter_accessions_indices()
#	print ""
#	print numAcc - len(accIndicesToKeep), "accessions removed from genotype data, leaving", \
#		len(accIndicesToKeep), "accessions in all."
#
#
#	print "Filtering phenotype data."
#	phed.removeAccessions(phenAccIndicesToKeep) #Removing accessions that don't have genotypes or phenotype values
#
#	#Ordering accessions according to the order of accessions in the genotype file
#	accessionMapping = []
#	i = 0
#	for acc in snpsds[0].accessions:
#		if acc in phed.accessions:
#			accessionMapping.append((phed.accessions.index(acc), i))
#			i += 1
#	phed.orderAccessions(accessionMapping)
#
#
#	if data_format == 'binary':
#		total_num = 0
#		removed_num = 0
#		for snpsd in snpsds:
#			total_num += len(snpsd.snps)
#			removed_num += snpsd.onlyBinarySnps()
#		print 'Removed %d non-binary SNPs out of %d SNPs' % (removed_num, total_num)


class _SnpsData_(object):
	"""
	05/11/2008 yh. add chromosome
	An abstract superclass.
	"""
	def __init__(self, snps, positions, baseScale=None, accessions=None, arrayIds=None, chromosome=None,
			alignment_positions=None, id=None, marker_types=None, alphabet=None, associated_positions=None):
		self.snps = snps #list[position_index][accession_index]
		self.positions = positions #list[position_index]
		if accessions:
			self.accessions = accessions #list[accession_index]
			#self._convert_to_tg_ecotypes_()
		if arrayIds:
			self.arrayIds = arrayIds #list[accession_index]
		self.chromosome = chromosome
		self.alignment_positions = alignment_positions
		self.id = id
		self.marker_types = marker_types #Where do these markers come frome, what type are they?  Useful for later analysis.
		self.alphabet = None
		self.missingVal = None
		self.associated_positions = associated_positions



	def mergeDataUnion(self, snpsd, priority=1, unionType=2, verbose=False):
		"""
		Merges two data, where the resulting data has the union of the SNPs or/and accessions

		unionType=
		1: All SNPs are included in the dataset  
		2: All accessions are included in the merged dataset 
		3: Both all SNPs and accessions are included in the merged dataset 
		
		It doesn't deal with multiple identical ecotypes/accession
		"""
		if verbose:
			print "Merging datas: id_1 =", self.id, "and id_2 =", snpsd.id
			print "Number of snps:", len(self.snps), "and", len(snpsd.snps)

		if self.id and snpsd.id:
			self.id = str(self.id) + "_merged_" + str(snpsd.id)

		#Find new accession indices
		newAccessions = []
		accessionsIndices = []
		commonAccessions = list(set(self.accessions).intersection(set(snpsd.accessions)))
		commonPositions = list(set(self.positions).intersection(set(snpsd.positions)))
		allAccessions = list(set(self.accessions).union(set(snpsd.accessions)))

		if unionType == 2 or unionType == 3:
			newAccessions = allAccessions
		elif unionType == 1:
			newAccessions = self.accessions

		for acc in newAccessions:
			index1 = -1
			index2 = -1
			if self.accessions.count(acc):
				index1 = self.accessions.index(acc)
			if snpsd.accessions.count(acc):
				index2 = snpsd.accessions.index(acc)
			accessionsIndices.append([index1, index2])


		if verbose:
			print "Number of common accessions:", len(commonAccessions)
			print "Total number of accessions:", len(allAccessions)
			print "Number of common Snps:", len(commonPositions)
			print "Only in 1st data set", list(set(self.accessions).difference(set(commonAccessions)))
			print "Only in 2st data set", list(set(snpsd.accessions).difference(set(commonAccessions)))

		snpErrorRate = []
		newSnps = []
		newPositions = []
		new_marker_types = []
		i = 0
		j = 0
		while i <= len(self.positions) and j <= len(snpsd.positions):
			if i < len(self.positions):
				pos1 = self.positions[i]
			if j < len(snpsd.positions):
				pos2 = snpsd.positions[j]
			if i < len(self.positions) and pos1 < pos2:
				newPositions.append(pos1)
				if self.marker_types and snpsd.marker_types:
					new_marker_types.append(self.marker_types[i])
				newSnp = []
				oldSnp = self.snps[i]
				for index in accessionsIndices:
					if index[0] != -1:
						newSnp.append(oldSnp[index[0]])
					else:
						newSnp.append(self.missingVal)
				newSnps.append(newSnp)
				i = i + 1
			elif j < len(snpsd.positions) and pos2 < pos1:
				if unionType == 1 or unionType == 3:
					newPositions.append(pos2)
					if self.marker_types and snpsd.marker_types:
						new_marker_types.append(snpsd.marker_types[j])
					newSnp = []
					oldSnp = snpsd.snps[j]
					for index in accessionsIndices:
						if index[1] != -1:
							newSnp.append(oldSnp[index[1]])
						else:
							newSnp.append(self.missingVal)
					newSnps.append(newSnp)
				j = j + 1
			elif i < len(self.positions) and j < len(snpsd.positions) and pos1 == pos2:
				counts = 0
				fails = 0
				for index in accessionsIndices:
					if index[0] != -1 and index[1] != -1:
						snp1 = self.snps[i][index[0]]
						snp2 = snpsd.snps[j][index[1]]
						if snp1 != self.missingVal and snp2 != self.missingVal:
							counts += 1
							if snp1 != snp2:
								fails = fails + 1
				error = 0
				if counts > 0:
					error = float(fails) / float(counts)
				snpErrorRate.append(error)

				newPositions.append(pos1)
				if self.marker_types and snpsd.marker_types:
					print "len(self.marker_types):", len(self.marker_types), ", i:", i
					print "len(snpsd.marker_types):", len(snpsd.marker_types), ", j:", j
					if self.marker_types[i] != snpsd.marker_types[j]:
						warn_str = "Different types of markers were merged.  Marker type 1: " + str(self.marker_types[i]) + ", Marker type 2: " + str(snpsd.marker_types[j]) + " ."
						warnings.warn(warn_str)
						new_marker_types.append(str(self.marker_types[i]) + "," + str(snpsd.marker_types[j]))
					else:
						new_marker_types.append(self.marker_types[i])

				newSnp = []
				oldSnp1 = self.snps[i]
				oldSnp2 = snpsd.snps[j]
				for index in accessionsIndices:
					if index[0] != -1 and oldSnp1[index[0]] != self.missingVal and priority == 1:
						newSnp.append(oldSnp1[index[0]])
					elif index[0] != -1 and oldSnp1[index[0]] != self.missingVal and priority == 2:
						if index[1] != -1:
							if oldSnp2[index[1]] == self.missingVal:
								newSnp.append(oldSnp1[index[0]])
							else:
								newSnp.append(oldSnp2[index[1]])
						else:
							newSnp.append(oldSnp1[index[0]])
					elif index[1] != -1:
						newSnp.append(oldSnp2[index[1]])
					else:
						newSnp.append(self.missingVal)
				newSnps.append(newSnp)
				i = i + 1
				j = j + 1

			else:
				# One pointer has reached the end and the end and the other surpassed it, i.e. we only need to copy the remaining one..
				while i < len(self.positions):
					newPositions.append(self.positions[i])
					newSnp = []
					oldSnp = self.snps[i]
					for index in accessionsIndices:
						if index[0] != -1:
							newSnp.append(oldSnp[index[0]])
						else:
							newSnp.append(self.missingVal)
					newSnps.append(newSnp)
					i = i + 1

				while j < len(snpsd.positions):
					if unionType == 1 or unionType == 3:
						newPositions.append(snpsd.positions[j])
						newSnp = []
						oldSnp = snpsd.snps[j]
						for index in accessionsIndices:
							if index[1] != -1:
								newSnp.append(oldSnp[index[1]])
							else:
								newSnp.append(self.missingVal)
						newSnps.append(newSnp)
					j = j + 1

				break

		if verbose:
			if snpErrorRate :
				print "Mean Snp Error:", sum(snpErrorRate) / float(len(snpErrorRate))
			print "Number of SNPs in merged data:", len(newSnps)
			print "Number of accessions in merged data:", len(newAccessions)


		if not (self.marker_types and snpsd.marker_types):
			new_marker_types = None
		self.marker_types = new_marker_types
		self.snps = newSnps
		self.positions = newPositions
		self.accessions = newAccessions
		self.arrayIds = None




	def mergeDataIntersection(self, snpsd, priority=1, intersectionType=2, verbose=False):
		"""
		Merges two data, where the resulting data has the intersection of the SNPs or/and accessions

		intersectionType=
		1: Only common SNPs are included in the dataset  (not implemented)
		2: Only common accessions are included in the merged dataset 
		3: Only common SNPs and accessions are included in the merged dataset (not implemented)
		
		It doesn't deal with multiple identical ecotypes/accession
		"""

		if verbose:
			print "Merging datas: id_1 =", self.id, "and id_2 =", snpsd.id
			print "Number of snps:", len(self.snps), "and", len(snpsd.snps)

		if self.id and snpsd.id:
			self.id = str(self.id) + "_merged_" + str(snpsd.id)

		#Find new accession indices
		newAccessions = []
		accessionsIndices = []
		commonAccessions = list(set(self.accessions).intersection(set(snpsd.accessions)))
		commonPositions = list(set(self.positions).intersection(set(snpsd.positions)))
		allAccessions = list(set(self.accessions).union(set(snpsd.accessions)))

		if intersectionType == 2 or intersectionType == 3:
			newAccessions = commonAccessions
		elif intersectionType == 1:
			newAccessions = self.accessions

		for acc in newAccessions:
			index1 = -1
			index2 = -1
			if self.accessions.count(acc):
				index1 = self.accessions.index(acc)
			if snpsd.accessions.count(acc):
				index2 = snpsd.accessions.index(acc)
			accessionsIndices.append([index1, index2])

		if verbose:
			print "Number of common accessions:", len(commonAccessions)
			print "Total number of accessions:", len(allAccessions)
			print "Number of common Snps:", len(commonPositions)
			print "Only in 1st data set", list(set(self.accessions).difference(set(commonAccessions)))
			print "Only in 2st data set", list(set(snpsd.accessions).difference(set(commonAccessions)))

		snpErrorRate = []
		newSnps = []
		newPositions = []
		i = 0
		j = 0
		while i <= len(self.positions) and j <= len(snpsd.positions):
			if i < len(self.positions):
				pos1 = self.positions[i]
			if j < len(snpsd.positions):
				pos2 = snpsd.positions[j]
			if i < len(self.positions) and pos1 < pos2:
				if intersectionType == 2:
					if self.marker_types and snpsd.marker_types:
						new_marker_types.append(self.marker_types[i])
					newPositions.append(pos1)
					newSnp = []
					oldSnp = self.snps[i]
					for index in accessionsIndices:
						if index[0] != -1:
							newSnp.append(oldSnp[index[0]])
						else:
							newSnp.append(self.missingVal)
					newSnps.append(newSnp)
				i = i + 1
			elif j < len(snpsd.positions) and pos2 < pos1:
				j = j + 1
			elif i < len(self.positions) and j < len(snpsd.positions) and pos1 == pos2:
				counts = 0
				fails = 0
				for index in accessionsIndices:
					if index[0] != -1 and index[1] != -1:
						snp1 = self.snps[i][index[0]]
						snp2 = snpsd.snps[j][index[1]]
						if snp1 != self.missingVal and snp2 != self.missingVal:
							counts += 1
							if snp1 != snp2:
								fails = fails + 1
				error = 0
				if counts > 0:
					error = float(fails) / float(counts)
				snpErrorRate.append(error)

				newPositions.append(pos1)
				if self.marker_types and snpsd.marker_types:
					if self.marker_types[i] != snpsd.marker_types[j]:
						warnings.warn("Different types of markers were merged.  Marker type 1: "
							+ str(self.marker_types[i]) + ", Marker type 2: " + str(snpsd.marker_types[j]) + " .")
						new_marker_types.append(self.marker_types[i])
					else:
						new_marker_types.append(self.marker_types[i])

				newSnp = []
				oldSnp1 = self.snps[i]
				oldSnp2 = snpsd.snps[j]
				for index in accessionsIndices:
					if index[0] != -1 and oldSnp1[index[0]] != self.missingVal and priority == 1:
						newSnp.append(oldSnp1[index[0]])
					elif index[0] != -1 and oldSnp1[index[0]] != self.missingVal and priority == 2:
						if index[1] != -1:
							if oldSnp2[index[1]] == self.missingVal:
								newSnp.append(oldSnp1[index[0]])
							else:
								newSnp.append(oldSnp2[index[1]])
						else:
							newSnp.append(oldSnp1[index[0]])
					elif index[1] != -1:
						newSnp.append(oldSnp2[index[1]])
					else:
						newSnp.append(self.missingVal)
				newSnps.append(newSnp)
				i = i + 1
				j = j + 1

			else:
				# One pointer has reached the end and the end and the other surpassed it, i.e. we only need to copy the remaining one..
				while i < len(self.positions):
					if intersectionType == 2:
						newPositions.append(self.positions[i])
						newSnp = []
						oldSnp = self.snps[i]
						for index in accessionsIndices:
							if index[0] != -1:
								newSnp.append(oldSnp[index[0]])
							else:
								newSnp.append(self.missingVal)
						newSnps.append(newSnp)
					i = i + 1

				while j < len(snpsd.positions):
					j = j + 1

				break



		if verbose:
			if len(snpErrorRate):
				print "Mean Snp Error:", sum(snpErrorRate) / float(len(snpErrorRate))
			print "Number of SNPs in merged data:", len(newSnps)
			print "Number of accessions in merged data:", len(newAccessions)

		if not (self.marker_types and snpsd.marker_types):
			new_marker_types = None
		self.snps = newSnps
		self.positions = newPositions
		self.accessions = newAccessions
		self.arrayIds = None

	def merge_data(self, sd, union_accessions=True, allow_multiple_markers=True, error_threshold=0.2):
		"""
		Merges data, allowing multiple markers at a position. (E.g. deletions and SNPs.)
		However it merges markers which overlap to a significant degree.
		"""
		if union_accessions:
			new_accessions = list(set(self.accessions).union(set(sd.accessions)))
		else:
			new_accessions = self.accessions
		acc_map = []
		for acc in new_accessions:
			try:
				ai1 = self.accessions.index(acc)
			except:
				ai1 = -1
			try:
				ai2 = sd.accessions.index(acc)
			except:
				ai2 = -1
			acc_map.append((ai1, ai2))
		#print acc_map
		#pdb.set_trace()
		index_dict = {}
		j = 0
		last_pos = sd.positions[j]
		for i, pos in enumerate(self.positions):
			curr_pos = last_pos
			while j < len(sd.positions) and curr_pos < pos:
				j += 1
				if j < len(sd.positions):
					curr_pos = sd.positions[j]
			last_pos = curr_pos
			index_list = []
			while j < len(sd.positions) and curr_pos == pos:
				index_list.append(j)
				j += 1
				if j < len(sd.positions):
					curr_pos = sd.positions[j]
			if index_list:
				index_dict[i] = index_list


		indices_to_skip = set()
		new_snps = []
		merge_count = 0
		for i, snp1 in enumerate(self.snps):
			new_snp = []
			if i in index_dict:
				index_list = index_dict[i]
				for j in index_list:
					error_count = 0
					t_count = 0
					snp2 = sd.snps[j]
					for (ai1, ai2) in acc_map:
						if ai1 != -1 and ai2 != -1:
							if snp1[ai1] != snp2[ai2] and snp1[ai1] != self.missingVal\
										and snp2[ai2] != self.missingVal:
								error_count += 1
							t_count += 1
					if t_count > 0 and error_count / float(t_count) < error_threshold:
						#print "Merge error is %f"%(error_count/float(t_count))
						new_snp = []
						for (ai1, ai2) in acc_map:
							if ai1 != -1 and ai2 != -1:
								if snp1[ai1] != self.missingVal:
									new_snp.append(snp1[ai1])
								else:
									new_snp.append(snp2[ai2])
							elif ai1 == -1:
								new_snp.append(snp2[ai2])
							else:
								new_snp.append(snp1[ai1])
						merge_count += 1
						indices_to_skip.add(j)
					elif t_count > 0:
						print "Not merging since error is %f." % (error_count / float(t_count))
						new_snp = []
						for (ai1, ai2) in acc_map:
							if ai1 != -1:
								new_snp.append(snp1[ai1])
							else:
								new_snp.append(self.missingVal)

			else:
				for (ai1, ai2) in acc_map:
					if ai1 != -1:
						new_snp.append(snp1[ai1])
					else:
						new_snp.append(self.missingVal)
				#print snp1,new_snp
			if new_snp == []:
				raise Exception
			new_snps.append(new_snp)
		new_positions = self.positions

		for j in range(len(sd.snps)):
			if not j in indices_to_skip:
				snp2 = sd.snps[j]
				new_snp = []
				for (ai1, ai2) in acc_map:
					if ai2 != -1:
						new_snp.append(snp2[ai2])
					else:
						new_snp.append(self.missingVal)
				if new_snp == []:
					raise Exception
				new_snps.append(new_snp)
				new_positions.append(sd.positions[j])



		pos_snp_list = zip(new_positions, new_snps)
		pos_snp_list.sort()
		r = map(list, zip(*pos_snp_list))
		self.positions = r[0]
		self.snps = r[1]
		self.accessions = new_accessions
		if len(self.snps) != len(self.positions):
			raise Exception
		print "Merged %d SNPs!" % (merge_count)
		print "Resulting in %d SNPs in total" % len(self.snps)



	def merge_data_old(self, snpsd, union_accessions=True):
		"""
		Merges data, allowing multiple markers at a position. (E.g. deletions and SNPs.)
		"""
		if union_accessions:
			new_accessions = list(set(self.accessions).union(set(snpsd.accessions)))
		else:
			raise NotImplementedError

		acc_map = []
		for acc in new_accessions:
			try:
				ai1 = self.accessions.index(acc)
			except:
				ai1 = -1
			try:
				ai2 = snpsd.accessions.index(acc)
			except:
				ai2 = -1
			acc_map.append((ai1, ai2))

		self.snps.extend(snpsd.snps)
		self.positions.extend(snpsd.positions)
		origin_list = [0] * len(self.snps) + [1] * len(snpsd.snps)
		pos_snp_list = zip(self.positions, self.snps, origin_list)
		pos_snp_list.sort()
		new_snps = []
		for i, (pos, snp, origin) in enumerate(pos_snp_list):
			new_snp = []
			for ai in acc_map:
				nt = snp[ai[origin]]
				if nt != -1:
					new_snp.append(nt)
				else:
					new_snp.append(self.missingVal)
			new_snps.append(new_snp)
		self.snps = new_snps
		self.accessions = new_accessions



#	def mergeData(self,snpsd, priority=1, debug=0):
#		"""
#
#		Merges two RawSnpsData objects.
#
#		If snps disagree, then the snps from the object called from is used.		
#		"""
#		sys.stderr.write("Merging datas Number of snps: %s vs %s ..."%(len(self.snps),len(snpsd.snps)))
#		# Find common accession indices
#		accessionsIndices = []
#		commonAccessions = list(set(self.accessions).intersection(set(snpsd.accessions)))
#		allAccessions = list(set(self.accessions).union(set(snpsd.accessions)))
#
#		for i in range(0,len(self.accessions)):
#			acc1 = self.accessions[i]
#			for k in range(0,len(snpsd.accessions)):
#				acc2 = snpsd.accessions[k]
#				if acc1==acc2:
#					accessionsIndices.append([i,k])
#   
#		if debug:
#			sys.stderr.write("Common accessions: %s\n"% len(commonAccessions))
#			sys.stderr.write("All accessions: %s\n"%len(allAccessions))
#			#print "Only in 1st data set", list(set(self.accessions).difference(set(commonAccessions)))
#			#print "Only in 2st data set", list(set(snpsd.accessions).difference(set(commonAccessions)))
#			print snpsd.accessions
#			print len(snpsd.accessions), len(list(set(snpsd.accessions)))
#			
#		commonSnpsPos = []
#		snpErrorRate = []
#		goodSnpsCounts = []
#		i = 0
#		j = 0 
#		while i <= len(self.positions) and j <= len(snpsd.positions):
#			if i < len(self.positions):
#				pos1 = self.positions[i]
#			if j < len(snpsd.positions):
#				pos2 = snpsd.positions[j] 
#			if i < len(self.positions) and pos1 < pos2: #Do nothing
#				i = i+1
#			elif j < len(snpsd.positions) and pos2 < pos1:  #Do nothing
#				j = j+1
#			elif i < len(self.positions) and j < len(snpsd.positions) and pos1==pos2:
#				commonSnpsPos.append(pos1)
#				counts = 0
#				fails = 0
#				for index in accessionsIndices:
#					snp1 = self.snps[i][index[0]]
#					snp2 = snpsd.snps[j][index[1]]
#					if snp1!=self.missingVal and snp2!=self.missingVal:
#						counts += 1
#						if snp1!=snp2:
#							fails = fails+1
#
#					if self.snps[i][index[0]] == self.missingVal or priority==2 and snpsd.snps[j][index[1]] != self.missingVal:
#						self.snps[i][index[0]]=snpsd.snps[j][index[1]]		   
#
#				goodSnpsCounts.append(counts)
#				error = 0
#				if counts>0:
#					error = float(fails)/float(counts)
#				snpErrorRate.append(error)										
#				i = i+1
#				j = j+1
#				action ="pos2=pos1"
#			else: 
#				# One pointer has reached the end so we're done...
#				break
#		sys.stderr.write("In all %s common snps found.\n"%len(commonSnpsPos))
#		sys.stderr.write("In all %s common accessions found.\n"%len(commonAccessions))
#		sys.stderr.write("Mean Snp Error: %s.\n"%(sum(snpErrorRate)/float(len(snpErrorRate))) )



	def extend(self, snpsd, union_accessions=True):
		"""
		Extends the marker data with another snpsd, but assumes they are non overlapping 
		and spatially ordered.  
		"""
		if union_accessions:
			new_accessions = list(set(self.accessions).union(set(snpsd.accessions)))
		else:
			raise NotImplementedError

		acc_map = []
		for acc in new_accessions:
			try:
				ai1 = self.accessions.index(acc)
			except:
				ai1 = -1
			try:
				ai2 = snpsd.accessions.index(acc)
			except:
				ai2 = -1
			acc_map.append((ai1, ai2))

		new_snps = []
		new_positions = []
		for i, snp in enumerate(self.snps):
			new_snp = []
			for (ai1, ai2) in acc_map:
				if ai1 == -1:
					new_snp.append(self.missingVal)
				else:
					new_snp.append(snp[ai1])
			new_positions.append(self.positions[i])

		for i, snp in enumerate(snpsd.snps):
			new_snp = []
			for (ai1, ai2) in acc_map:
				if ai2 == -1:
					new_snp.append(self.missingVal)
				else:
					new_snp.append(snp[ai2])
			new_positions.append(snpsd.positions[i])
		self.snps = new_snps
		self.positions = new_positions
		self.accessions = new_accessions


	def sample_snps(self, random_fraction):
		"""
		Discards SNPs, leaving a random fraction of them untouched.
		"""
		import random
		snp_pos_list = random.sample(zip(self.positions, self.snps), int(random_fraction * len(self.snps)))
		snp_pos_list.sort()
		(self.positions, self.snps) = map(list, zip(*snp_pos_list))


	def scalePositions(self, baseScale):
		for i in range(0, len(self.positions)):
			self.positions[i] = int(self.positions[i] * baseScale)
		self.baseScale = baseScale


	def addSnp(self, snp):
		self.snps.append(snp)


	def addPos(self, position):
		self.positions.append(position)


	def _convert_to_tg_ecotypes_(self):
		import phenotypeData as pd
		e_dict = pd._getEcotype2TgEcotypeDict_()
		new_ecotypes = []
		conversion_count = 0
		for et in self.accessions:
			if int(et) != e_dict[int(et)]:
				conversion_count += 1
			new_ecotypes.append(str(e_dict[int(et)]))
		self.accessions = new_ecotypes
		if conversion_count:
			print conversion_count, "ecotypes were converted to tg_ecotypes."


	def removeAccessions(self, accessions, arrayIds=None):
		"""
		Removes accessions from the data.
		"""
		accPair = set()
		for i in range(0, len(accessions)):
			if arrayIds:
				key = (accessions[i], arrayIds[i])
			else:
				key = accessions[i]
			accPair.add(key)
		newAccessions = []
		newArrayIds = []
		for i in range(0, len(self.snps)):
			snp = self.snps[i]
			newSnp = []
			for j in range(0, len(self.accessions)):
				if arrayIds:
					key = (self.accessions[j], self.arrayIds[j])
				else:
					key = self.accessions[j]
				if not key in accPair:
					newSnp.append(snp[j])
					if i == 0:
						newAccessions.append(self.accessions[j])
						if arrayIds:
							newArrayIds.append(self.arrayIds[j])
			self.snps[i] = newSnp
		self.accessions = newAccessions
		if arrayIds:
			self.arrayIds = newArrayIds


	def remove_accessions(self, accessions_to_keep, use_accession_names=False):
		"""
		Uses accession names!
		"""

		if use_accession_names:
			import phenotypeData as pd
			ad = pd._getAccessionToEcotypeIdDict_(accessions_to_keep)
			ecotypes = []
			for acc in accessions_to_keep:
				ecotypes.append(ad[acc])
		else:
			ecotypes = accessions_to_keep

		indices_to_keep = []
		str_accessions = map(str, self.accessions)
		for e in ecotypes:
			try:
				indices_to_keep.append(str_accessions.index(str(e)))
			except:
				pass

		self.removeAccessionIndices(indices_to_keep)


	def removeAccessionIndices(self, indicesToKeep):
		"""
		Removes accessions from the data.
		"""
		newAccessions = []
		newArrayIds = []
		for i in indicesToKeep:
			newAccessions.append(self.accessions[i])
			if self.arrayIds:
				newArrayIds.append(self.arrayIds[i])
		for i in range(len(self.snps)):
			snp = self.snps[i]
			newSnp = []
			for j in indicesToKeep:
				newSnp.append(snp[j])
			self.snps[i] = newSnp
		self.accessions = newAccessions
		if self.arrayIds:
			#print "removeAccessionIndices: has array IDs: self.arrayIds =",self.arrayIds
			self.arrayIds = newArrayIds
			#print "len(self.arrayIds):",len(self.arrayIds)
		#pdb.set_trace()
		#print "len(self.accessions):",len(self.accessions)


	def filterMonoMorphicSnps(self):
		"""
		05/12/08 yh. add no_of_monomorphic_snps_removed
		Removes SNPs from the data which are monomorphic.
		"""
		newPositions = []
		newSnps = []
		for i in range(0, len(self.positions)):
			count = 0
			for nt in self.alphabet:
				if nt in self.snps[i]:
					count += 1
			if count > 1:
				newSnps.append(self.snps[i])
				newPositions.append(self.positions[i])
		numRemoved = len(self.positions) - len(newPositions)
		self.no_of_monomorphic_snps_removed = numRemoved
		self.snps = newSnps
		self.positions = newPositions
		return numRemoved

	def onlyBinarySnps(self):
		"""
		Removes all but binary SNPs.  (I.e. monomorphic, tertiary and quaternary alleles SNPs are removed.)
		"""
		newPositions = []
		newSnps = []
		for i in range(0, len(self.positions)):
			count = 0
			for nt in self.alphabet:
				if nt in self.snps[i]:
					count += 1
			if count == 2:
				newSnps.append(self.snps[i])
				newPositions.append(self.positions[i])
		numRemoved = len(self.positions) - len(newPositions)
		self.no_of_nonbinary_snps_removed = numRemoved
		self.snps = newSnps
		self.positions = newPositions
		print "Removed %d non-binary SNPs, leaving %d SNPs in total." % (numRemoved, len(self.snps))

		return numRemoved



	def orderAccessions(self, accessionMapping):
		newAccessions = [None for i in range(0, len(accessionMapping))]
		for (i, j) in accessionMapping:
			newAccessions[j] = self.accessions[i]

		newSnps = []
		for pos in self.positions:
			newSnps.append([self.missingVal for i in range(0, len(accessionMapping))])

		for (i, j) in accessionMapping:
			for posIndex in range(0, len(self.positions)):
				newSnps[posIndex][j] = self.snps[posIndex][i]

		self.accessions = newAccessions
		self.snps = newSnps

		if self.arrayIds:
			newArrayIds = [None for i in range(0, len(self.arrayIds))]
			for (i, j) in accessionMapping:
				newArrayIds[j] = self.arrayIds[i]
			self.arrayIds = newArrayIds




	def filterRegion(self, startPos, endPos):
		"""
		Filter all SNPs but those in region.
		"""
		i = 0
		while i < len(self.snps) and self.positions[i] < startPos:
			i += 1

		newSnps = []
		newPositions = []
		while i < len(self.snps) and self.positions[i] < endPos:
			newSnps.append(self.snps[i])
			newPositions.append(self.positions[i])
			i += 1
		self.snps = newSnps
		self.positions = newPositions

	def get_region_snpsd(self, start_pos, end_pos):
		"""
		Constructs and returns a new object based on the prev. using SNPs only in the given region.
		"""
		import copy
		snpsd = copy.deepcopy(self) #Clone
		i = 0
		while i < len(self.snps) and self.positions[i] < start_pos:
			i += 1

		new_snps = []
		new_positions = []
		while i < len(self.snps) and self.positions[i] < end_pos:
			new_snps.append(self.snps[i])
			new_positions.append(self.positions[i])
			i += 1
		snpsd.positions = new_positions
		snpsd.snps = new_snps
		return snpsd



	def get_snp_phen_pair(self, snp_index, phen_data, phen_index=None, missingVal=None):
		"""
		returns a pair of snp and phenotype values, removing NA's
		"""
		if not missingVal:
			missingVal = self.missingVal
		snp = self.snps[snp_index]
		phen_vals = phen_data.getPhenVals(phen_index, noNAs=False)
#		print phen_data.accessions
#		print self.accessions
#		print len(phen_vals),len(phen_data.accessions)
		acc_set = list(set(self.accessions).intersection(set(phen_data.accessions)))
		acc_set.sort() #Not really necessary...
		new_snp = []
		new_phen_vals = []
		for acc in acc_set:
			si = self.accessions.index(acc)
			pi = phen_data.accessions.index(acc)
			if snp[si] != missingVal and phen_vals[pi] != missingVal:
				new_snp.append(snp[si])
				new_phen_vals.append(phen_vals[pi])
		return (new_snp, new_phen_vals)


	def get_mafs(self):
		"""
		Returns MAFs and MARFs
		
		Assumes that this is a binary allele.
		"""
		mafs = []
		marfs = []
		for snp in self.snps:
			missing_count = list(snp).count(self.missingVal)
			num_nts = len(snp) - missing_count
			nts = set(snp)
			if missing_count:
				nts.remove(self.missingVal)
			c = list(snp).count(nts.pop()) / float(num_nts)
			if c > 0.5:
				c = 1.0 - c
			marfs.append(c)
			mafs.append(int(c * num_nts))
		return {"mafs":mafs, "marfs":marfs}





	def remove_snps(self, remove_indices):
		"""
		Filters all SNPs in the indices list.
		"""
		remove_indices.sort(reverse=True)
		for i in remove_indices:
			del self.snps[i]
			del self.positions[i]
			if self.alignment_positions:
				del self.alignment_positions[i]


	def filter_snp_indices(self, indices):
		"""
		Filters all SNPs but those in the indices list.
		"""
		new_snps = []
		new_positions = []
		new_alignm_positions = []
		for i in indices:
			new_snps.append(self.snps[i])
			new_positions.append(self.positions[i])
		self.snps = new_snps
		self.positions = new_positions
		if self.alignment_positions:
			new_alignm_positions = []
			for i in indices:
				new_alignm_positions.append(self.alignment_positions[i])
			self.alignment_positions = new_alignm_positions


	def filter_accessions_by_NAs(self, maxMissing):
		"""
		Removes the accessions in the snpsd if it has NA-rates greater than maxMissing. 
		"""
		sys.stderr.write("Calculating NA rate ...")
		missingCounts = self.accessionsMissingCounts()

		indices_to_keep = []
		num_snps = float(len(self.snps))
		for i, mc in enumerate(missingCounts):
			if (mc / num_snps) <= maxMissing:
				indices_to_keep.append(i)
		self.removeAccessionIndices(indices_to_keep)


	def filter_na_snps(self, max_na_rate=0.0):
		new_snps = []
		new_positions = []
		num_snps = len(self.snps)
		e_num = float(len(self.accessions))
		for i, snp in enumerate(self.snps):
			if snp.count(self.missingVal) / e_num <= max_na_rate:
				new_snps.append(snp)
				new_positions.append(self.positions[i])
		self.snps = new_snps
		self.positions = new_positions
		print "Removed %d SNPs, leaving %d in total." % (num_snps - len(self.snps), len(self.snps))



	def accessionsMissingCounts(self):
		"""
		Returns a list of accessions and their missing value rates.
		"""
		missingCounts = [0] * len(self.accessions)

		for i in range(0, len(self.positions)):
			for j in range(0, len(self.accessions)):
				if self.snps[i][j] == self.missingVal:
					missingCounts[j] += 1

		return missingCounts

	def __str__(self):
		st = "Number of accessions: " + str(len(self.accessions)) + "\n"
		st += "Number of SNPs: " + str(len(self.snps)) + "\n"
		uniqueAccessions = list(set(self.accessions))
		if len(uniqueAccessions) < len(self.accessions):
			st += "Not all accessions are unique. \n" + "Number of unique accessions: " + str(len(uniqueAccessions)) + "\n"
			for acc in uniqueAccessions:
				count = self.accessions.count(acc)
				if count > 1:
					st += acc + " has " + str(count) + " occurrences.\n\n"
		return st


	def countMissingSnps(self):
		"""
		Returns a list of accessions and their missing value rates.
		"""
		missingCounts = [0.0 for acc in self.accessions]
		totalCounts = [0.0 for acc in self.accessions]
		for i in range(0, len(self.positions)):
			for j in range(0, len(self.accessions)):
				totalCounts[j] += 1
				if self.snps[i][j] == self.missingVal:
					missingCounts[j] += 1

		return [missingCounts[i] / totalCounts[i] for i in range(len(self.accessions))], missingCounts


	def print_ecotype_info_to_file(self, filename):
		"""
		Prints info on the accessions/ecotype IDs to a file.
		"""
		import phenotypeData as pd
		eid = pd._getEcotypeIdInfoDict_()
		f = open(filename, 'w')
		f.write((", ".join(["Ecotype_id", "Native_name", "Stock_parent", "latitude", "longitude", "country"])) + "\n")
		for e in self.accessions:
			l = [e] + list(eid[int(e)])
			l = map(str, l)
			f.write((", ".join(l)) + "\n")
		f.close()





class Marker():
	"""
	A class for marker info. 
	"""
	def __init__(self, position, chromosome, accessions=None, alleles=None, type=None, marker_data=None):
		self.position = position
		self.chromosome = chromosome
		self.accessions = accessions
		self.alleles = alleles
		if not alleles and marker_data:
			self._getAllele_(marker_data)


	def _getAllele_(self, marker_data):
		chr = self.chromosome - 1
		snpsd = snpsds[chr]
		self.snpsdIndex = -1
		self.alleles = None
		try:
			self.marker_data_index = marker_data.positions.index(self.positions)
			self.alleles = marker_data.snps[self.marker_data_index]

		except Exception:
			print "The corresponding allele was not found in the marker data."




class MarkerData(_SnpsData_):
	"""
	A class to encompass general marker data, as well as handle 
	representative markers, which represent other markers, through LD.
	"""

	def __init__(self, *args, **kwargs):
		linked_markers = kwargs.pop("linked_markers", None)
		_SnpsData_.__init__(self, *args, **kwargs) #Calling the parent's init fun.





class RawDecoder(dict):
	def __missing__(self, key):
		return 'NA'
	def __init__(self, initdict={}):
		for letter in ['A', 'C', 'G', 'T']:
			self[letter] = letter
		for letter in initdict:
			self[letter] = initdict[letter]



class RawSnpsData(_SnpsData_):
	"""
	05/11/2008 yh. give default values to all initial arguments so that it could be called without any arguments.
	Similar to the SnpsData class, except it deals with bases (A,C,G,T), instead of allele numbers (0s and 1s).
	
	Alphabet: A,C,G,T, and NA

	"""
	#alphabet = ['A','C','G','T']
	#alphabet = ['A','C','G','T','-']




	def __init__(self, snps=None, positions=None, baseScale=None, accessions=None, arrayIds=None,
			chromosome=None, callProbabilities=None, alignment_positions=None, id=None,
			marker_types=None, missing_val='NA'):
		self.snps = snps
		self.positions = positions
		self.accessions = accessions
		#if accessions:
		#	self._convert_to_tg_ecotypes_()

		self.arrayIds = arrayIds
		self.callProbabilities = []  #list[position_index][accession_index]
		if callProbabilities:
			self.callProbabilities = callProbabilities
		self.alignment_positions = alignment_positions
		self.missingVal = missing_val
		self.alphabet = IUPAC_alphabet
		self.marker_types = marker_types #Where do these markers come frome, what type are they?  Useful for later analysis.
		self.id = id
		self.chromosome = chromosome


	def writeToFile(self, filename, chromosome, withArrayId=False):
		"""
		Writes data to a file.  1 file for each chromosome.

		WARNING OLD, outdated!
		"""
		outStr = ""
		fieldStrings = ["Chromosome", "Positions"]
		for acc in self.accessions:
			fieldStrings.append(str(acc))
		outStr += ", ".join(fieldStrings) + "\n"
		for i in range(0, len(self.positions)):
			outStr += str(chromosome) + ", " + str(self.positions[i]) + ", " + ", ".join(self.snps[i]) + "\n"
		f = open(filename, "w")
		f.write(outStr)
		f.close()



	def impute_data_region(self, reference_snpsd, window_count=100, window_size=None, verbose=False, **kwargs):
		"""
		Impute any NA SNPs.
		"""
		import ImputeSNPs as imp
		import tempfile, os, copy
		if window_size:
			start_pos = max(0, self.positions[0] - window_size)
			end_pos = self.positions[-1] + window_size
			snpsd = reference_snpsd.get_region_snpsd(start_pos, end_pos)
		else:
			snpsd = copy.deepcopy(reference_snpsd)
		if (not window_count) and window_size:
			i = 0
			while i < len(snpsd.positions) and snpsd.positions[i] < self.positions[0]:
				i += 1
			suffix_count = i
			i = 1
			while i >= len(snpsd.positions) and snpsd.positions[len(snpsd.positions) - i] > self.positions[-1]:
				i += 1
			window_count = max(suffix_count + 1, i)

		snpsd.mergeDataUnion(self, priority=2, unionType=3, verbose=verbose)
		snpsd.na_ambigious_calls()
		snpsd.onlyBinarySnps()
		print snpsd.snps
		print snpsd.accessions
		print len(snpsd.accessions), len(set(snpsd.accessions))
		tmpFile1 = tempfile.mkstemp()
		os.close(tmpFile1[0])
		tmpFile2 = tempfile.mkstemp()
		os.close(tmpFile2[0])

		pdb.set_trace()
		if verbose:
			print "Preparing data in", tmpFile1[1]
		imp.writeAsNputeFile(snpsd, tmpFile1[1])
		imp.checkNputeFile(tmpFile1[1])
		nputeCmd = "python " + imp.path_NPUTE + "NPUTE.py -m 0 -w " + str(window_count) + " -i " + str(tmpFile1[1]) + " -o " + str(tmpFile2[1])
		if verbose:
			print "Imputing data..."
			print nputeCmd
		os.system(nputeCmd)
		if verbose:
			print "Imputation done!"

		snpsd = imp.readNputeFile(tmpFile2[1], snpsd.accessions, snpsd.positions)
		self.mergeDataUnion(snpsd, priority=2)
		os.remove(tmpFile1[1])
		os.remove(tmpFile2[1])



	def impute_data(self, verbose=True):
		"""
		Impute any NAs in the data.
		"""
		import ImputeSNPs as imp
		import tempfile, os
		tmpFile1 = tempfile.mkstemp()
		os.close(tmpFile1[0])
		tmpFile2 = tempfile.mkstemp()
		os.close(tmpFile2[0])

		if verbose:
			print "Preparing data in", tmpFile1[1]
		imp.writeAsNputeFile(self, tmpFile1[1])
		imp.checkNputeFile(tmpFile1[1])
		nputeCmd = "python " + imp.path_NPUTE + "NPUTE.py -m 0 -w 30 -i " + str(tmpFile1[1]) + " -o " + str(tmpFile2[1])
		if verbose:
			print "Imputing data..."
			print nputeCmd
		os.system(nputeCmd)
		if verbose:
			print "Imputation done!"
		pdb.set_trace()
		snpsd = imp.readNputeFile(tmpFile2[1], self.accessions, self.positions)
		self.snps = snpsd.snps
		os.remove(tmpFile1[1])
		os.remove(tmpFile2[1])


	def compareWith(self, snpsd, withArrayIds=0, verbose=True, heterozygous2NA=False):
		"""
		05/10/2008 len(commonSnpsPos) could be zero
		This function performs QC on two datasets.

		Requires accessions to be defined.

		withArrayIds = 0 (no array IDs), =1 the object called from has array IDs, =2 both objects have array IDs 
		"""
		basicAlphabet = ['-', 'A', 'C', 'G', 'T']

		if verbose:
			print "Comparing datas"
			print "Number of snps:", len(self.snps), "and", len(snpsd.snps)
		# Find common accession indices
		accessionsIndices = []
		accessionErrorRate = []
		accessionCallRates = [[], []]
		accessionCounts = []
		commonAccessions = []
		arrayIds = []
		for i in range(0, len(self.accessions)):
			acc = self.accessions[i]
			if snpsd.accessions.count(acc):
				j = snpsd.accessions.index(acc)
				accessionsIndices.append([i, j])
				accessionErrorRate.append(0)
				accessionCounts.append(0)
				accessionCallRates[0].append(0)
				accessionCallRates[1].append(0)
				commonAccessions.append(acc)
				if withArrayIds > 0:
					if withArrayIds == 1:
						aId = self.arrayIds[i]
					elif withArrayIds == 2:
						aId = (self.arrayIds[i], snpsd.arrayIds[j])
					arrayIds.append(aId)
		commonSnpsPos = []
		snpErrorRate = []
		snpCallRate = [[], []]
		goodSnpsCounts = []
		totalCounts = 0
		totalFails = 0
		i = 0
		j = 0
		while i <= len(self.positions) and j <= len(snpsd.positions):
			if i < len(self.positions):
				pos1 = self.positions[i]
			if j < len(snpsd.positions):
				pos2 = snpsd.positions[j]
			if i < len(self.positions) and pos1 < pos2:
				i = i + 1
			elif j < len(snpsd.positions) and pos2 < pos1:
				j = j + 1
			elif i < len(self.positions) and j < len(snpsd.positions) and pos1 == pos2:
				commonSnpsPos.append(pos1)
				fails = 0
				counts = 0
				missing1 = 0
				missing2 = 0
				for k in range(0, len(accessionsIndices)):
					accIndices = accessionsIndices[k]
					snp1 = self.snps[i][accIndices[0]]
					snp2 = snpsd.snps[j][accIndices[1]]

					if heterozygous2NA:
						if snp1 in basicAlphabet and snp2 in basicAlphabet:
							accessionCounts[k] += 1
							counts += 1
							if snp1 != snp2:
								fails = fails + 1
								accessionErrorRate[k] += 1
						else:
							if snp1 == self.missingVal:
								accessionCallRates[0][k] += 1
								missing1 += 1
							if snp2 == self.missingVal:
								accessionCallRates[1][k] += 1
								missing2 += 1
					else:
						if snp1 != self.missingVal and snp2 != self.missingVal:
							accessionCounts[k] += 1
							counts += 1
							if snp1 != snp2:
								fails = fails + 1
								accessionErrorRate[k] += 1
						else:
							if snp1 == self.missingVal:
								accessionCallRates[0][k] += 1
								missing1 += 1
							if snp2 == self.missingVal:
								accessionCallRates[1][k] += 1
								missing2 += 1
				goodSnpsCounts.append(counts)
				error = 0
				totalCounts += counts
				totalFails += fails
				if counts > 0:
					error = float(fails) / float(counts)
					snpErrorRate.append(error)
				snpCallRate[0].append(missing1 / float(len(accessionsIndices)))
				snpCallRate[1].append(missing2 / float(len(accessionsIndices)))
				i = i + 1
				j = j + 1
			else:
				# One pointer has reached and the end and the other surpassed it.
				break

		for i in range(0, len(accessionErrorRate)):
			if accessionCounts[i] > 0:
				accessionErrorRate[i] = accessionErrorRate[i] / float(accessionCounts[i])
			no_of_common_snps_pos = len(commonSnpsPos)
			if no_of_common_snps_pos > 0:	#05/08/10 yh
				accessionCallRates[0][i] = accessionCallRates[0][i] / float(no_of_common_snps_pos)
				accessionCallRates[1][i] = accessionCallRates[1][i] / float(no_of_common_snps_pos)
			else:
				accessionCallRates[0][i] = 0
				accessionCallRates[1][i] = 1

		print "In all", len(snpErrorRate), "common snps found"
		print "In all", len(commonAccessions), "common accessions found"
		print "Common accessions IDs :", commonAccessions
		print "Common SNPs positions :", commonSnpsPos
		print "Accessions error rates", accessionErrorRate
		print "Average Accession SNP Error:", sum(accessionErrorRate) / float(len(accessionErrorRate))
		print "SNP error rates", snpErrorRate
		print "Average Snp Error:", sum(snpErrorRate) / float(len(snpErrorRate))


		naCounts1 = [0] * len(accessionsIndices)
		for i in range(0, len(self.positions)):
			for k in range(0, len(accessionsIndices)):
				accIndices = accessionsIndices[k]
				snp = self.snps[i][accIndices[0]]
				if snp == self.missingVal:
					naCounts1[k] += 1

		naCounts2 = [0] * len(accessionsIndices)
		for i in range(0, len(snpsd.positions)):
			for k in range(0, len(accessionsIndices)):
				accIndices = accessionsIndices[k]
				snp = snpsd.snps[i][accIndices[1]]
				if snp == self.missingVal:
					naCounts2[k] += 1


		return [commonSnpsPos, snpErrorRate, commonAccessions, accessionErrorRate, accessionCallRates, arrayIds, accessionCounts, snpCallRate, [naCounts1, naCounts2], [totalCounts, totalFails]]

	def getSnpsData(self, missingVal= -1, reference_ecotype='6909', only_non_binary=True):
		"""
		Returns a SnpsData object correspoding to this RawSnpsData object.

		Note that some of the SnpsData attributes are a shallow copy of the RawSnpsData obj.

		Reference ecotype is set to be the 0 allele. 
		"""
		decoder = {self.missingVal:missingVal} #Might cause errors somewhere???!!!

		if reference_ecotype in self.accessions:
			ref_i = self.accessions.index(reference_ecotype)
		else:
			ref_i = 0
			import warnings
			warnings.warn("Given reference ecotype %s wasn't found, using %s as 0-reference." % (reference_ecotype, self.accessions[ref_i]))
		snps = []
		positions = []
		num_lines = len(self.accessions)
		for i, snp in enumerate(self.snps):
			alphabet = self.alphabet[:]

			if snp[ref_i] != 'NA':
				decoder[snp[ref_i]] = 0
				alphabet.remove(snp[ref_i])
				k = 1
			else:
				k = 0
			for nt in alphabet:
				if nt in snp:
					decoder[nt] = k
					k = k + 1
			if only_non_binary:
				if k == 2:
					positions.append(self.positions[i])
					new_snp = sp.zeros(num_lines, dtype='int8')
					for j, nt in enumerate(snp):
						new_snp[j] = decoder[nt]
					snps.append(new_snp)
				else:
					continue
			else:
				if k > 2:
					max1 = 0
					maxnt1 = ''
					max2 = 0
					maxnt2 = ''
					for nt in alphabet:
						c = snp.count(nt)
						if c > max1:
							max1 = c
							maxnt1 = nt
						elif c > max2:
							max2 = c
							maxnt2 = nt
						decoder[nt] = missingVal
					decoder[maxnt1] = 0
					decoder[maxnt2] = 1
				new_snp = []
				for nt in snp:
					new_snp.append(decoder[nt])
				snps.append(new_snp)

		if only_non_binary:
			print 'Removed %d non-binary SNPs out of %d, when converting to binary SNPs.'\
				% (len(self.positions) - len(positions), len(self.positions))
		else:
			positions = self.positions

		assert len(snps) == len(positions), 'Somthing odd with the lengths.'
		return SNPsData(snps, positions, accessions=self.accessions, marker_types=self.marker_types, missing_val=missingVal)


	def filterBadSnps(self, snpsd, maxNumError=0):
		"""
		05/12/08 add no_of_snps_filtered_by_mismatch
		Filters snps with high rate mismatches, when compared against another snpsd.
		"""

		newSnps = []
		newPositions = []
		sys.stderr.write("Comparing datas. Number of snps: %s vs %s. \n" % (len(self.snps), len(snpsd.snps)))
		# Find common accession indices
		accessionsIndices = []
		commonAccessions = []
		for i in range(0, len(self.accessions)):
			acc = self.accessions[i]
			if snpsd.accessions.count(acc):
				j = snpsd.accessions.index(acc)
				accessionsIndices.append([i, j])
				commonAccessions.append(acc)


		commonSnpsPos = []
		i = 0
		j = 0
		while i <= len(self.positions) and j <= len(snpsd.positions):
			if i < len(self.positions):
				pos1 = self.positions[i]
			if j < len(snpsd.positions):
				pos2 = snpsd.positions[j]
			if i < len(self.positions) and pos1 < pos2:
				newSnps.append(self.snps[i])
				newPositions.append(self.positions[i])
				i = i + 1
			elif j < len(snpsd.positions) and pos2 < pos1:
				j = j + 1
			elif i < len(self.positions) and j < len(snpsd.positions) and pos1 == pos2:
				commonSnpsPos.append(pos1)
				fails = 0
				counts = 0
				for k in range(0, len(accessionsIndices)):
					accIndices = accessionsIndices[k]
					snp1 = self.snps[i][accIndices[0]]
					snp2 = snpsd.snps[j][accIndices[1]]
					if snp1 != self.missingVal and snp2 != self.missingVal:
						counts += 1
						if snp1 != snp2:
							fails = fails + 1
				error = 0
				if counts > 0:
					error = float(fails) / float(counts)
				if error <= maxNumError:
					newSnps.append(self.snps[i])
					newPositions.append(self.positions[i])
				i = i + 1
				j = j + 1
			else:
				# One pointer has reached and the end and the other surpassed it.
				break
		self.no_of_snps_filtered_by_mismatch = len(self.snps) - len(newSnps)
		sys.stderr.write('%s SNPs were filtered\n' % (self.no_of_snps_filtered_by_mismatch))
		self.snps = newSnps
		self.positions = newPositions


	def convertBadCallsToNA(self, minCallProb=0):
		"""
		Converts all base calls with call prob. lower than the given one to NAs.
		"""
		totalCount = len(self.positions) * len(self.snps[0])
		count = 0
		for i in range(0, len(self.positions)):
			for j in range(0, len(self.snps[i])):
				if self.callProbabilities[i][j] < minCallProb:
					self.snps[i][j] = self.missingVal
					count += 1
		fractionConverted = count / float(totalCount)
		return fractionConverted





	def getSnpsNArates(self):
		"""
		Returns a list of accessions and their missing value rates.
		"""
		naRates = []
		for i in range(0, len(self.positions)):
			missingCounts = 0
			totalCounts = 0
			for j in range(0, len(self.accessions)):
				totalCounts += 1
				if self.snps[i][j] == self.missingVal:
					missingCounts += 1
			naRates.append(missingCounts / float(totalCounts))

		return naRates


	def filterMissingSnps(self, maxNumMissing=0):
		"""
		05/12/2008 add no_of_snps_filtered_by_na
		Removes SNPs from the data which have more than maxNumMissing missing values.
		"""
		warnings.warn("This function is deprecated!!")
		newPositions = []
		newSnps = []
		for i in range(0, len(self.positions)):
			missingCount = 0
			for nt in self.snps[i]:
				if nt == self.missingVal:
					missingCount += 1
			if missingCount <= maxNumMissing:
				newSnps.append(self.snps[i])
				newPositions.append(self.positions[i])
		numRemoved = len(self.positions) - len(newPositions)
		self.no_of_snps_filtered_by_na = numRemoved
		self.snps = newSnps
		self.positions = newPositions
		return numRemoved


	def na_ambigious_calls(self, maxNumMissing=0):
		"""
		NAs any ambiguous calls, i.e. anything but A,C,G,T,-
		"""
		for snp in self.snps:
			for i, nt in enumerate(snp):
				if not nt in ['A', 'C', 'G', 'T', '-', self.missingVal]:
					snp[i] = self.missingVal


	def filterMinMAF(self, minMAF=0):
		"""
		Removes SNPs from the data which have a minor allele count less than the given one.
		"""
		warnings.warn("This function is old, and should be replaced.. please don't use it!!")
		newPositions = []
		newSnps = []
		ntCounts = [0.0] * len(self.alphabet)
		ntl = self.alphabet
		for i in range(0, len(self.positions)):
			snp = self.snps[i]
			totalCount = 0
			for j in range(0, len(self.alphabet)):
				nt = ntl[j]
				c = snp.count(nt)
				ntCounts[j] = c
				totalCount += c
			if totalCount == 0 :
				print "snp:", snp
				print "self.alphabet:", self.alphabet
			maxAF = max(ntCounts)
			maf = 0
			if ntCounts.count(maxAF) < 2:
				for j in range(0, len(self.alphabet)):
					if ntCounts[j] < maxAF and ntCounts[j] > 0:
						if ntCounts[j] > maf:
							maf = ntCounts[j]
			else:
				maf = maxAF

			if minMAF <= maf:
				newSnps.append(self.snps[i])
				newPositions.append(self.positions[i])

		numRemoved = len(self.positions) - len(newPositions)
		self.snps = newSnps
		self.positions = newPositions
		return numRemoved



	def filterMinMARF(self, minMARF=0):
		"""
		Removes SNPs from the data which have a Relative MAF less than the given one.
		"""
		warnings.warn("This function is old, and should be replaced.. please don't use it!!")
		newPositions = []
		newSnps = []
		ntCounts = [0.0] * len(self.alphabet)
		ntl = self.alphabet
		for i in range(0, len(self.positions)):
			snp = self.snps[i]
			totalCount = 0
			for j in range(0, len(self.alphabet)):
				nt = ntl[j]
				c = snp.count(nt)
				ntCounts[j] = c
				totalCount += c
			if totalCount == 0 :
				print "snp:", snp
				print "self.alphabet:", self.alphabet
			for j in range(0, len(self.alphabet)):
				ntCounts[j] = ntCounts[j] / float(totalCount)
			maxAF = max(ntCounts)
			maf = 0
			if ntCounts.count(maxAF) < 2:
				for j in range(0, len(self.alphabet)):
					if ntCounts[j] < maxAF and ntCounts[j] > 0.0:
						if ntCounts[j] > maf:
							maf = ntCounts[j]
			else:
				maf = maxAF

			if minMARF <= maf:
				newSnps.append(self.snps[i])
				newPositions.append(self.positions[i])

		numRemoved = len(self.positions) - len(newPositions)
		self.snps = newSnps
		self.positions = newPositions
		return numRemoved



	def snpsFilterMAF(self, maf_thresholds):
		self.filterMinMARF(maf_thresholds[0])
		if maf_thresholds[1] < 0.5:
			raise Exception





	def mergeIdenticalAccessions(self, accessionIndexList, priority):
		"""
		The priority argument gives particular accessions in the list priority
		if priority is set to 0, then majority (if any) rules.
		"""
		pass


class SNPsData(_SnpsData_):
	"""
	Efficient genotype data, using the numpy class.
	
	An alternative to the old SnpsData class.
	"""
	alphabet = [-1, 0, 1, 2, 3]  #Here -1 is thought to be missing data value.
	def __init__(self, snps, positions, accessions=None, arrayIds=None, chromosome=None,
			alignment_positions=None, id=None, marker_types=None, missing_val= -1):
		self.snps = snps
		self.positions = positions
		self.accessions = accessions
		self.arrayIds = arrayIds
		self.chromosome = chromosome
		self.alignment_positions = alignment_positions
		self.marker_types = marker_types #Where do these markers come frome, what type are they?  Useful for later analysis.
		self.id = id
		self.missingVal = missing_val



	def removeAccessionIndices(self, indicesToKeep):
		"""
		Removes accessions from the data.
		"""
		newAccessions = []
		newArrayIds = []
		for i in indicesToKeep:
			newAccessions.append(self.accessions[i])
			if self.arrayIds:
				newArrayIds.append(self.arrayIds[i])
		num_accessions = len(indicesToKeep)
		for i in range(len(self.snps)):
			self.snps[i] = self.snps[i][indicesToKeep]
#			snp = self.snps[i]
#			newSnp = sp.empty(num_accessions,dtype='int8')
#			for j,k in enumerate(indicesToKeep):
#				newSnp[j] = snp[k]
#			self.snps[i] = newSnp
		self.accessions = newAccessions
		if self.arrayIds:
			#print "removeAccessionIndices: has array IDs: self.arrayIds =",self.arrayIds
			self.arrayIds = newArrayIds
			#print "len(self.arrayIds):",len(self.arrayIds)
		#pdb.set_trace()
		#print "len(self.accessions):",len(self.accessions)


	def onlyBinarySnps(self):
		"""
		Removes all but binary SNPs.  (I.e. monomorphic, tertiary and quaternary alleles SNPs are removed.)
		"""
		new_positions = []
		new_snps = []
		for i, (snp, pos) in enumerate(izip(self.snps, self.positions)):
			if len(sp.unique(snp)) == 2:
				new_snps.append(snp)
				new_positions.append(pos)
		num_removed = len(self.positions) - len(new_positions)
		self.no_of_nonbinary_snps_removed = num_removed
		self.snps = new_snps
		self.positions = new_positions
		#print "Removed %d non-binary SNPs, leaving %d SNPs in total." % (num_removed, len(self.snps))
		return num_removed


	def haplotize(self, snp_window=None, base_window=None):

		if snp_window:

			haplotypes_list = []
			i = 0
			num_accessions = len(self.accessions)
			while i < snp_window:
				cur_snps = self.snps[:i + snp_window + 1]
				haplotypes_list.append(get_haplotypes(cur_snps, num_accessions))
				i += 1
			while i < len(self.snps) - snp_window:
				cur_snps = self.snps[i - snp_window:i + snp_window + 1]
				haplotypes_list.append(get_haplotypes(cur_snps, num_accessions))
				i += 1
				if i % 100 == 0: print i
			while i < len(self.snps):
				cur_snps = self.snps[i - snp_window:]
				haplotypes_list.append(get_haplotypes(cur_snps, num_accessions))
				i += 1
			self.snps = haplotypes_list

		elif base_window:
			raise NotImplementedError
		else:
			raise Exception('Window size missing')


	def get_mafs(self, w_missing=False, binary=True):
		"""
		Returns MAFs and MARFs
		
		(Uses numpy.bincount)
		"""

#		def _get_maf_(snp):
#			l = sp.bincount(sp.unique(snp, False, True)[1])
#			maf = min(l)
#			return maf


		def _get_maf_(snp):
			l = sp.bincount(snp)
			maf = max(l)
			for m in l[1:]:
				if m != 0 and m < maf:
					maf = m
			return maf

		mafs = []
		marfs = []
		num_nts = len(self.snps[0])
		if w_missing:
			for snp in self.snps:
				if self.missingVal in snp:
					missing_count = list(snp).count(self.missingVal)
					num_nts = len(snp) - missing_count
					nts = set(snp)
					nts.remove(self.missingVal)
					c = list(snp).count(nts.pop()) / float(num_nts)
					if c > 0.5:
						c = 1.0 - c
					marfs.append(c)
					mafs.append(int(c * num_nts))
				else:
					l = sp.bincount(snp)
					maf = min(l)
					mafs.append(maf)
					marfs.append(maf / float(num_nts))
		else:
			if binary:
				for snp in self.snps:
					l = sp.bincount(snp)
					maf = min(l)
					mafs.append(maf)
					marfs.append(maf / float(num_nts))
			else:
				for snp in self.snps:
					maf = _get_maf_(snp)
					mafs.append(maf)
					marfs.append(maf / float(num_nts))

		return {"mafs":mafs, "marfs":marfs}




	def filter_mac(self, min_mac=15, w_missing=False):
       		"""
       		Filter minor allele count SNPs.
       		"""
       		new_snps = []
       		new_positions = []
		if w_missing and self.missingVal in snp:
			for snp, pos in izip(self.snps, self.positions):
				missing_count = list(snp).count(self.missingVal)
				num_nts = len(snp) - missing_count
				nts = set(snp)
				nts.remove(self.missingVal)
				mac = list(snp).count(nts.pop())
				if mac > num_nts * 0.5: mac = num_nts - mac

				if mac >= min_mac:
					new_snps.append(snp)
					new_positions.append(pos)
		else:
			for snp, pos in izip(self.snps, self.positions):
				if min(sp.bincount(snp)) >= min_mac:
					new_snps.append(snp)
					new_positions.append(pos)

		print 'Removed %d SNPs out of %d, leaving %d SNPs.' % (len(self.positions) - len(new_positions),
								len(self.positions), len(new_positions))
		self.positions = new_positions
		self.snps = new_snps







class SnpsData(_SnpsData_):
	"""
	A class for SNPs data.  It uses 0, 1 (and 2, 3 if there are more than 2 alleles) to represent SNPs.  
	-1 is used if the allele data is missing.

	This class should really be named marker data, as it is used for more than just SNP markers.

	Contains various functions that aid the analysis of the data.
	"""
	#freqs = [] #list[position_index1][position_index2-position_index1+1] Linkage frequencies. 
	#baseScale = 1 #Scaling for positions
	alphabet = [0, 1, 2, 3]
	def __init__(self, snps, positions, baseScale=None, accessions=None, arrayIds=None, chromosome=None,
			alignment_positions=None, id=None, marker_types=None, missing_val= -1):
		self.snps = snps
		self.positions = positions
		if baseScale:
			self.scalePositions(baseScale)
		self.accessions = accessions
		self.arrayIds = arrayIds
		self.chromosome = chromosome
		self.alignment_positions = alignment_positions
		self.marker_types = marker_types #Where do these markers come frome, what type are they?  Useful for later analysis.
		self.id = id
		self.missingVal = missing_val


	def clone(self):
		"""
		Filter all SNPs but those in region.
		"""
		newSNPs = []
		newPositions = []
		newAccessions = []
		for acc in self.accessions:
			newAccessions.append(acc)


		for i in range(0, len(self.positions)):
			new_snp = []
			for j in range(0, len(self.accessions)):
				new_snp.append(self.snps[i][j])
			newSNPs.append(new_snp)
			newPositions.append(self.positions[i])

		return SnpsData(newSNPs, newPositions, accessions=newAccessions, arrayIds=self.arrayIds, chromosome=self.chromosome)


	def calc_r2_levels(self, window_size=100000, pic_file=None):
		"""
		Returns data on r2 levels of this snps data.
		
		requires calcFreqs function to be defined.  
		"""
		print "calculating r2's"
		freqs_list = self.calcFreqsUnbiased(window_size)
		i = 0
		mean_r2_list = []
		mean_pos_list = []
		while i < len(self.snps) - 1:
			r2_list = []
			j = i + 1
			curr_pos = self.positions[i]
			print "Now at position", curr_pos
			while j < len(self.snps) - 1 and self.positions[j] < curr_pos + window_size:
				j += 1
			mean_pos_list.append((curr_pos + self.positions[j]) / 2.0)
			while i < j:
				for freqs in freqs_list[i]:
					r2_list.append(r2(freqs))
				i += 1
			if len(r2_list):
				mean_r2 = sum(r2_list) / float(len(r2_list))
			else:
				mean_r2 = 0
			mean_r2_list.append(r2_list)

		print zip(mean_r2_list, mean_pos_list)
		return {"mean_r2_list":mean_r2_list, "mean_pos_list":mean_pos_list}


	def _calc_r2_(self, si1, si2):
		"""
		Calculate the r2 between the two SNPs.
		"""
		raise NotImplementedError


	def calc_r2_matrix(self, w_missing=False):
		"""
		Calculate the r2 matrix for all of these SNPs.
		"""
		r2s = sp.ones((len(self.snps), len(self.snps)))
		freqs_list = self.calc_all_freqs(w_missing=w_missing)
		for i in range(len(self.snps)):
			for j in range(i):
				r2s[i, j] = r2(freqs_list[i][j])
				r2s[j, i] = r2s[i, j]
		return r2s


	def calc_r2_with_snp(self, snp):
		"""
		returns a list of r2 values between the snpsd snps and the given snp.
		"""
		freqs = self._calc_freqs_with_snp_(snp)
		return [r2(f) for f in freqs]



	def _calc_freqs_with_snp_(self, snp):
		"""
		returns a list of r2 values between the snpsd snps and the given snp.
		"""
		freqs = 	[]
		for snp1 in self.snps:
			jfreqs = [0.0] * 4   #The frequencies, of a pair of alleles.
			count = 	0.0
			for nt1, nt2 in zip(snp1, snp):
				if nt1 != self.missingVal and nt2 != self.missingVal:
					count += 1.0
					val = nt1 * 2 + nt2
					jfreqs[val] += 1.0
			jfreqs = [f / count for f in jfreqs]
			if not 0.99 < sum(jfreqs) < 1.01:
				raise Exception("freqs, don't sum up!!")
			freqs.append(jfreqs)
		return freqs


	def remove_redundant_snps(self, r2_threshold=0.8, w_missing=False):
		"""
		Removes any redundant SNPs as measured by correlation.  Returns information on what SNPs are correlated.
		"""
		r2s = self.calc_r2_matrix(w_missing=w_missing)
		snps_indices = range(len(self.snps))
		new_snps_indices = []
		redundant_snps_indices_list = []
		redundant_pos_list = []
		while snps_indices:
			si1 = snps_indices.pop()
			new_snps_indices.append(si1)
			redundant_snps_indices = []
			redundant_positions = []
			for si2 in snps_indices[:]:
				if r2s[si1, si2] >= r2_threshold:
					redundant_snps_indices.append(si2)
					redundant_positions.append(self.positions[si2])
					snps_indices.remove(si2)

			redundant_pos_list.append(redundant_positions)
			redundant_snps_indices_list.append(redundant_snps_indices)
		print "In total", len(self.snps), ", thereof", len(self.snps) - len(new_snps_indices), "are redundant."
		self.filter_snp_indices(new_snps_indices)
		self.redundant_positions = redundant_pos_list
		return redundant_pos_list



	def calcFreqs(self, windowSize, innerWindowSize=0):
		"""
		Returns a list of two loci comparison frequencies with in a window.
		"""
		freqs = 	[]
		delta = 0
		if len(self.snps) > 0:
			delta = 1.0 / float(len(self.snps[0]))
		for i in xrange(0, len(self.snps) - 1):
			l = []
			j = i + 1
			while j < len(self.snps) and (self.positions[j] - self.positions[i]) < innerWindowSize :
				j = j + 1
			while j < len(self.snps) and (self.positions[j] - self.positions[i]) <= windowSize:
				jfreqs = [0.0] * 4   #The frequencies, of a pair of alleles.
				snp1 = self.snps[i]
				snp2 = self.snps[j]
				count = 	0
				for k in xrange(0, len(snp1)):
					val = snp1[k] * 2 + snp2[k]
					jfreqs[val] = jfreqs[val] + delta
				l.append(jfreqs)
				j = j + 1
			freqs.append(l)
		self.freqs = freqs
		return self.freqs

	def calc_all_freqs(self, w_missing=False):
		"""
		Returns a matrix of two loci comparison frequencies.
		"""
		if w_missing:
			freqs = 	[]
			for i, snp1 in enumerate(self.snps):
				l = []
				for j in range(i):
					jfreqs = [0.0] * 4   #The frequencies, of a pair of alleles.
					snp2 = self.snps[j]
					count = 	0.0
					for nt1, nt2 in zip(snp1, snp2):
						if nt1 != self.missingVal and nt2 != self.missingVal:
							count += 1.0
							val = nt1 * 2 + nt2
							jfreqs[val] += 1.0
					jfreqs = [f / count for f in jfreqs]
					l.append(jfreqs)
				freqs.append(l)
			return freqs
		else:
			freqs = 	[]
			delta = 0
			if len(self.snps) > 0:
				delta = 1.0 / float(len(self.snps[0]))
			for i, snp1 in enumerate(self.snps):
				l = []
				for j in range(i):
					jfreqs = [0.0] * 4   #The frequencies, of a pair of alleles.
					snp2 = self.snps[j]
					for k in range(len(snp1)):
						val = snp1[k] * 2 + snp2[k]
						jfreqs[val] = jfreqs[val] + delta
					l.append(jfreqs)
				freqs.append(l)

			return freqs

	def calcFreqsUnbiased(self, windowSize, innerWindowSize=0): # Returns a list of two loci comparison frequencies with in a window.
		""" 
		Uses a distribution of frequencies that is not dependent on the lenght of the sequence.  (Alot of data is disregarded.)
		"""
		numPairs = 0  #Counts the number of pairs
		freqs = 	[]
		delta = 0
		if len(self.snps) > 0:
			delta = 1.0 / float(len(self.snps[0]))
		for i in xrange(0, len(self.snps) - 1):
			if (self.positions[len(self.snps) - 1] - self.positions[i]) >= windowSize:
				j = i + 1
				l = []
				while j < len(self.snps) and (self.positions[j] - self.positions[i]) < innerWindowSize :
					j = j + 1
				while  j < len(self.snps) and (self.positions[j] - self.positions[i]) <= windowSize:
					jfreqs = [0.0] * 4
					snp1 = self.snps[i]
					snp2 = self.snps[j]
					count = 	0
					for k in xrange(0, len(snp1)):
						val = snp1[k] * 2 + snp2[k]
						jfreqs[val] = jfreqs[val] + delta
					l.append(jfreqs)
					j = j + 1
					numPairs = numPairs + 1
				freqs.append(l)
			else:
				break
		self.freqs = freqs
		return self.freqs
		#return numPairs


	def calcFreqsSimple(self):
		freqs = 	[]
		delta = 0
		if len(self.snps) > 0:
			delta = 1.0 / float(len(self.snps[0]))
		for i in xrange(0, len(self.snps) - 1):
			l = []
			for j in xrange(i + 1, len(self.snps)):
				jfreqs = [0.0] * 4
				snp1 = self.snps[i]
				snp2 = self.snps[j]
				count = 	0
				for k in xrange(0, len(snp1)):
					val = snp1[k] * 2 + snp2[k]
					jfreqs[val] = jfreqs[val] + delta
				l.append(jfreqs)
			freqs.append(l)
		self.freqs = freqs
		return self.freqs


	def snpsFilter(self):
		"""
		Splits the SNPs up after their allele frequency ..
		"""
		snpsDatas = [SnpsData([], []), SnpsData([], []), SnpsData([], []), SnpsData([], []), SnpsData([], []), SnpsData([], []), SnpsData([], []), SnpsData([], []), SnpsData([], []), SnpsData([], [])]

		l = len(self.snps[0])
		for j in xrange(0, len(self.snps)):
			snp = self.snps[j]
			c = snp.count(0) / float(l)
			"""
			if c>0.5:
			c = 1-c
			if c ==	0.5:
			c =0.499
			"""
			if c == 	1:
				c = 0.99999
			i = int(c * 10)
			snpsDatas[i].addSnp(snp)
			snpsDatas[i].addPos(self.positions[j])
		return snpsDatas


	def snpsFilterMAF(self, mafs):
		"""
		Filters all snps with MAF not in the interval out of dataset.
		"""
		newsnps = []
		newpos = []
		#print self.snps
		l = len(self.snps[0])
		if l == 0:
			print self.snps
		for j in xrange(0, len(self.snps)):
			snp = self.snps[j]
			c = snp.count(0) / float(l)
			if c > 0.5:
				c = 1 - c
			if c > mafs[0] and c <= mafs[1]:
				newsnps.append(snp)
				newpos.append(self.positions[j])
		if len(newsnps) == 0:
			print "Filtered out all snps from", len(self.snps), " what to do what to do?"
		del self.snps
		self.snps = newsnps
		del self.positions
		self.positions = newpos


	def snpsFilterRare(self, threshold=0.1):
		"""Filters all snps with MAF of less than threshold out of dataset."""
		newsnps = []
		newpos = []
		#print self.snps
		l = len(self.snps[0])
		if l == 0:
			print self.snps
		for j in xrange(0, len(self.snps)):
			snp = self.snps[j]
			c = snp.count(0) / float(l)
			if c > 0.5:
				c = 1 - c
			if c > threshold:
				newsnps.append(snp)
				newpos.append(self.positions[j])
		#if len(newsnps)==0:  
			#print "Filtered out all snps from",len(self.snps)," what to do what to do?"
		del self.snps
		self.snps = newsnps
		del self.positions
		self.positions = newpos


	def _genRecombFile(self, filename, windowSize, maxNumPairs):
		n = len(self.snps[0]) #number of individuals/accessions
		numPairs = self.calcFreqsUnbiased(windowSize)
		if numPairs != 0:
			filterProb = float(maxNumPairs) / numPairs
		else:
			return numPairs
		f = open(filename, 'w')
		numPairs = 0
		for i in range(0, len(self.freqs)):
			for j in range(0, len(self.freqs[i])):
				if random.random() <= filterProb:
					numPairs = numPairs + 1
					st = str(i + 1) + " " + str(j + i + 2) + " " + str(self.positions[j + i + 1] - self.positions[i]) + " u "  # u denotes unknown as opposed to ad ancestral derived
					st = st + str(int(self.freqs[i][j][0] * n + 0.5)) + " " + str(int(self.freqs[i][j][1] * n + 0.5)) + " "
					st = st + str(int(self.freqs[i][j][2] * n + 0.5)) + " " + str(int(self.freqs[i][j][3] * n + 0.5)) + " 0 0 0 0 " + str(100 - n) + "\n"
					f.write(st)
		f.close()
		return numPairs

	def estimateRecomb(self, windowSize, maxNumPairs=10000, tempfile1="tmp1", tempfile2="tmp2", meanTract=200, numPoints=50):
		num = self._genRecombFile(tempfile1, windowSize, maxNumPairs)
		if num < 1:
			return [0, 0, 0]
		os.system(homedir + "Projects/programs/Maxhap/maxhap 1 " + homedir + "Projects/programs/Maxhap/h100rho  .0008 10 .1 0.01 500 " + str(numPoints) + " " + str(meanTract) + " < " + tempfile1 + " > " + tempfile2)
		f = open(tempfile2, 'r')
		lines = f.readlines()
		npairs = float(lines[1].split()[2])
		i = 2
		while(lines[i].split()[0] == "Warning:"):
			i = i + 1
		rho = float(lines[i].split()[1])
		ratio = float(lines[i].split()[2])
		f.close()
		return [rho, ratio, npairs]


	def meanAF(self):
		""" Mean allele frequency. """
		if len(self.snps):
			l = float(len(self.snps[0]))
			c = 0
			for j in xrange(0, len(self.snps)):
				snp = self.snps[j]
				snpsc = snp.count(0)
				if snpsc < (l / 2.0):
					c = c + snpsc / l
				else:
					c = c + abs((l / 2.0) - snpsc) / l
			return c / len(self.snps)
		else:
			return 0

	def EHH(self, snp1, snp2):
		""" Calculates the EHH between two SNPs"""
		data = self.snps[snp1:snp2]
		haplotypes = []
		haplotypecount = []
		for i in range(0, len(self.snps[0])):
			haplotype = []
			for j in range(snp1, snp2 + 1):
				haplotype.append(self.snps[j][i])
			if not haplotype in haplotypes:
				haplotypes.append(haplotype)
				haplotypecount.append(1.0)
			else:
				k = haplotypes.index(haplotype)
				haplotypecount[k] = haplotypecount[k] + 1.0
		s = 0.0
		for i in range(0, len(haplotypes)):
			if haplotypecount[i] > 1:
				s = s + haplotypecount[i] * (haplotypecount[i] - 1)
		s = s / (len(self.snps[0]) * (len(self.snps[0]) - 1))
		return s

	def totalEHH(self, windowSize, innerWindowSize):
		""" 
		Lenght indep mean EHH statistics.. (Note: no data filtering!)
		"""
		ehhcount = 0
		ehh = 0.0
		for i in xrange(0, len(self.snps) - 1):
			if (self.positions[len(self.snps) - 1] - self.positions[i]) >= windowSize:
				j = i + 1
				l = []
				while j < len(self.snps) and (self.positions[j] - self.positions[i]) < innerWindowSize :
					j = j + 1
				while  j < len(self.snps) and (self.positions[j] - self.positions[i]) <= windowSize:
					ehh = ehh + self.EHH(i, j)
					ehhcount = ehhcount + 1
					j = j + 1
			else:
				break
		return [ehh, ehhcount]



class SNPsDataSet:
	"""
	A class that encompasses multiple _SnpsData_ chromosomes objects (chromosomes), and can deal with them as a whole.

	This object should eventually replace the snpsdata lists.. 
	"""
	snpsDataList = None
	chromosomes = None
	accessions = None

	def __init__(self, snpsds, chromosomes, id=None, is_binary=None, call_method=None, data_format=None):
		self.snpsDataList = snpsds
		self.chromosomes = chromosomes
		self.accessions = self.snpsDataList[0].accessions
		self.array_ids = self.snpsDataList[0].arrayIds
		self.id = id
		self.is_binary = is_binary
		self.missing_val = snpsds[0].missingVal
		self.call_method = call_method
		self.data_format = data_format
		if not id and snpsds[0].id:
				self.id = id
		for i in range(1, len(self.chromosomes)):
			if self.accessions != self.snpsDataList[i].accessions:
				raise Exception("Accessions (or order) are different between SNPs datas")
		if not is_binary:
			self.is_binary = list(snpsds[0].snps[0]).count(0) or list(snpsds[0].snps[0]).count(1)
		else:
			self.is_binary = is_binary





	def add_to_db(self, short_name, method_description='', data_description='', comment='', **kwargs):
		"""
		Other possible keyword args are parent_id, accession_set_id ,imputed ,unique_ecotype 
		"""
		import dbutils
		conn = dbutils.connect_to_papaya()
		cursor = conn.cursor()


		#Checking whether data is already inserted...
		sql_statement = "SELECT id FROM stock_250k.call_method WHERE short_name='%s';" % (short_name)
		print sql_statement
		cursor.execute(sql_statement)
		row = cursor.fetchone()
		print row
		if row:
			print "Data is already inserted in DB.  File will however be updated."
			call_method_id = int(row[0])
			filename = '/Network/Data/250k/db/dataset/call_method_%d.tsv' % call_method_id
		else:

			#Inserting data
			sql_statement = "INSERT INTO stock_250k.call_method (short_name,method_description,data_description,comment"
			for k in kwargs:
				sql_statement += ',%s' % k
			sql_statement += ") VALUES ('%s','%s','%s','%s'" % (short_name, method_description, data_description, comment)
			for k in kwargs:
				sql_statement += ',%s' % str(kwargs[k])
			sql_statement += ');'
			print sql_statement
			cursor.execute(sql_statement)

			#Getting method_id
			sql_statement = "SELECT id FROM stock_250k.call_method WHERE short_name='%s';" % (short_name)
			print sql_statement
			cursor.execute(sql_statement)
			row = cursor.fetchone()
			call_method_id = int(row[0])
			filename = '/Network/Data/250k/db/dataset/call_method_%d.tsv' % call_method_id

			#Updating the filename
			sql_statement = "UPDATE stock_250k.call_method SET filename='%s' WHERE id=%d" % (filename, call_method_id)
			print sql_statement
			cursor.execute(sql_statement)

			print "Committing transaction (making changes permanent)."
			conn.commit()

		print "Call method id is %d" % call_method_id
		#Generating Yu's file...
		self.write_to_file_yu_format(filename)



		#Generate DB call files...
		self._generate_db_call_files_(call_method=call_method_id, cursor=cursor, conn=conn)

		#Add SNPs to DB?


		cursor.close()
		conn.close()

		return call_method_id




	def _generate_db_call_files_(self, call_method=None, array_ids=None, file_dir='/Network/Data/250k/db/calls/', cursor=None, conn=None):
		import os
		import warnings

		if not array_ids:
			if not self.array_ids:
				raise Exception('Array IDs are missing.')
			else:
				array_ids = self.array_ids
		if not call_method:
			raise Exception("Call method is missing!!")
		chr_pos_snp_list = self.getChrPosSNPList()

		#Create the call method directory
		file_dir = file_dir + 'method_' + str(call_method) + '/'
		if not os.path.lexists(file_dir):
			os.mkdir(file_dir)
		else:
			warnings.warn('Directory already exists: %s' % file_dir)

		#Connect to DB, if needed
		if not cursor:
			import dbutils
			conn = dbutils.connect_to_papaya()
			cursor = conn.cursor()

		for i, aid in enumerate(array_ids):
			print "Inserting genotype files into DB."
			#Checking if it is already in the DB.
			sql_statement = "SELECT id FROM stock_250k.call_info WHERE array_id=%s AND method_id=%d;"\
					% (aid, call_method)
			print sql_statement
			cursor.execute(sql_statement)
			row = cursor.fetchone()
			if row:
				'Already in DB. File will be updated.'
				call_info_id = int(row[0])
			else:
				#Insert info
				sql_statement = "INSERT INTO stock_250k.call_info (array_id, method_id) VALUES (%s,%d);"\
						% (aid, call_method)
				print sql_statement
				cursor.execute(sql_statement)

				sql_statement = "SELECT id FROM stock_250k.call_info WHERE array_id=%s AND method_id=%d;"\
						% (aid, call_method)
				print sql_statement
				cursor.execute(sql_statement)
				row = cursor.fetchone()
				if row:
					call_info_id = int(row[0])
				print "Committing transaction (making changes permanent)."
				conn.commit()

			try:
				#Write to a designated place.
				file_name = file_dir + str(call_info_id) + "_call.tsv"
				sql_statement = "UPDATE stock_250k.call_info \
						 SET filename='%s'\
						 WHERE id=%d" % (file_name, call_info_id)
				print sql_statement
				cursor.execute(sql_statement)

				#Generate file in the right place
				f = open(file_name, 'w')
				f.write('SNP_ID\t%s\n' % aid)
				for (c, p, s) in chr_pos_snp_list:
					f.write('%d_%d\t%s\n' % (c, p, s[i]))
				f.close()
			except Exception, err_str:
				print "Couldn't generate call info file, error message:%s" % err_str
		print "Closing connection."
		#Close connection
		if not cursor:
			cursor.close()
			conn.close()
		print "Remember to copy files to papaya, i.e. everything in directory: %s" % file_dir




	def writeToFile(self, filename, delimiter=", ", missingVal="NA", accDecoder=None,
			withArrayIds=False, decoder=None, callProbFile=None, binary_format=False):
		"""
		Writes data to a file. 
		
		Note that there is no decoder dictionary option here..
		"""

		print "Writing data to file:", filename
		numSnps = 0
		for i in range(0, len(self.chromosomes)):
			numSnps += len(self.snpsDataList[i].positions)


		#outStr = "NumSnps: "+str(numSnps)+", NumAcc: "+str(len(accessions))+"\n"
		if withArrayIds:
			outStr = ", ".join(["-", "-"] + self.snpsDataList[0].arrayIds) + "\n"
		else:
			outStr = ""
		fieldStrings = ["Chromosome", "Positions"]
		if accDecoder:
			for acc in self.snpsDataList[i].accessions:
				fieldStrings.append(str(accDecoder[acc]))
		else:
			for acc in self.snpsDataList[i].accessions:
				fieldStrings.append(str(acc))
		outStr += delimiter.join(fieldStrings) + "\n"
		if binary_format:
			f = open(filename, "wb")
		else:
			f = open(filename, "w")
		f.write(outStr)
		if decoder:
			for chromosome, snpsd in izip(self.chromosomes, self.snpsDataList):
				sys.stdout.write(".")
				sys.stdout.flush()
				for pos, snp in izip(snpsd.positions, snpsd.snps):
					outStr = str(chromosome) + delimiter + str(pos)
					for nt in snp:
						outStr += delimiter + str(decoder[snp])
					outStr += "\n"
					f.write(outStr)
#				for j in range(0,len(self.snpsDataList[i].positions)):
#					outStr =""
#					outStr += str(self.chromosomes[i])+delimiter+str(self.snpsDataList[i].positions[j])
#					for k in range(0, len(self.snpsDataList[0].accessions)):
#						outStr += delimiter+str(decoder[self.snpsDataList[i].snps[j][k]])
#					outStr +="\n"
#					f.write(outStr)
		else:
			for chromosome, snpsd in izip(self.chromosomes, self.snpsDataList):
				sys.stdout.write(".")
				sys.stdout.flush()
				for pos, snp in izip(snpsd.positions, snpsd.snps):
					outStr = '%d%s%d%s%s\n' % (chromosome, delimiter, pos, delimiter,
								delimiter.join(map(str, snp.tolist())))
					f.write(outStr)
#			for i in range(0,len(self.chromosomes)):
#				sys.stdout.write(".")
#				sys.stdout.flush()
#				for j in range(0,len(self.snpsDataList[i].positions)):
#					outStr =""
#					outStr += str(self.chromosomes[i])+delimiter+str(self.snpsDataList[i].positions[j])
#					snp = self.snpsDataList[i].snps[j]
#					if len(snp) != len(self.snpsDataList[0].accessions):
#						print "len(snp):",len(snp),", vs. len(self.snpsDataList[0].accessions):",len(self.snpsDataList[0].accessions)
#						raise Exception("The length didn't match")
#					for nt in snp:
#						outStr += delimiter+str(nt)
#					outStr +="\n"
#					f.write(outStr)
		f.close()
		print ""

		if callProbFile:
			if withArrayIds:
				outStr = "-, -, " + ", ".join(self.snpsDataList[0].arrayIds) + "\n"
			else:
				outStr = ""
			f = open(callProbFile, "w")
			outStr += delimiter.join(fieldStrings) + "\n"
			f.write(outStr)
			f.flush()
			for i in range(0, len(self.chromosomes)):
				outStr = ""
				snpsd = self.snpsDataList[i]
				self.snpsDataList[i] = []
				for j in range(0, len(snpsd.positions)):
					outStr += str(self.chromosomes[i]) + delimiter + str(snpsd.positions[j])
					for k in range(0, len(snpsd.accessions)):
						outStr += delimiter + str(snpsd.callProbabilities[j][k])
					outStr += "\n"
				del snpsd
				f.write(outStr)
				f.flush()
			f.close()



	def write_to_file_yu_format(self, filename):
		"""
		Writes data to a file in Yu's format. (Requires array IDs)
		
		Only works with raw sequence data.. (IUPAC nucleotides)
		"""
		import yu_snp_key as yk
		print "transposing chr_pos_snp_list"
		chr_pos_snp_list = map(list, zip(*self.getChrPosSNPList()))
		chrs = chr_pos_snp_list[0]
		positions = chr_pos_snp_list[1]
		snps = chr_pos_snp_list[2]
		assert len(snps) == len(chrs) == len(positions), "SNPs, chromosomes, and positions not with same lenght"

		print "transposing SNPs"
		snps = map(list, zip(*snps))

		array_ids = self.array_ids
		ecotypes = self.accessions
		assert len(snps) == len(array_ids) == len(ecotypes), "SNP, array IDs, and ecotype IDs not with same lenght"

		print "len(array_ids) len(ecotypes):", len(array_ids), len(ecotypes)

		print "Writing data to file:", filename
		f = open(filename, "w")
		str_list = ['ecotype_id', 'array_id']
		str_list.extend([str(c) + "_" + str(p) for (c, p) in zip(chrs, positions)])
		f.write('\t'.join(str_list) + '\n')
		for ei, ai, nts in zip(ecotypes, array_ids, snps):
			str_list = [str(ei), str(ai)]
			str_list.extend([str(yk.nt_2_number[nt]) for nt in nts])
			f.write('\t'.join(str_list) + '\n')
		f.close()



	def coordinate_w_phenotype_data(self, phend, pid, coord_phen=True):

		"""
		Deletes accessions which are not common, and sorts the accessions, removes monomorphic SNPs, etc.
		"""
#		import bisect
		print "Coordinating SNP and Phenotype data."
		ets = phend.phen_dict[pid]['ecotypes']
		#Checking which accessions to keep and which to remove.
#		common_ets = list(set(self.accessions + ets))
#		common_ets.sort()

		sd_indices_to_keep = set()#[]
		pd_indices_to_keep = []

		for i, acc in enumerate(self.accessions):
			for j, et in enumerate(ets):
				if et == acc:
					sd_indices_to_keep.add(i)
					pd_indices_to_keep.append(j)

#
#		for i, acc in enumerate(self.accessions):
#			if common_ets[bisect.bisect(common_ets, acc) - 1] == acc:
#					sd_indices_to_keep.add(i)
#		for j, et in enumerate(ets):
#			if common_ets[bisect.bisect(common_ets, et) - 1] == et:
#					pd_indices_to_keep.append(j)
		sd_indices_to_keep = list(sd_indices_to_keep)
		sd_indices_to_keep.sort()



		#Filter accessions which do not have phenotype values (from the genotype data).
		print "Filtering genotype data"
		#if len(sd_indices_to_keep) != len(self.accessions):
		self.filter_accessions_indices(sd_indices_to_keep)
		if coord_phen:
			num_values = len(phend.phen_dict[pid]['ecotypes'])
			print "Filtering phenotype data."
			phend.filter_ecotypes(pd_indices_to_keep, pids=[pid]) #Removing accessions that don't have genotypes or phenotype values
			ets = phend.phen_dict[pid]['ecotypes']
			print "Out of %d, leaving %d values." % (num_values, len(ets))
		#Ordering accessions according to the order of accessions in the genotype file

#		if ets != self.accessions:
#			l = zip(ets, range(len(ets)))
#			l.sort()
#			l = map(list, zip(*l))
#			ets_map = l[1]
#			phend.order_ecotypes(ets_map, pids=[pid])


		if self.data_format == 'binary':
			print 'Filtering non-binary SNPs'
			total_num = 0
			removed_num = 0
			for snpsd in self.snpsDataList:
				total_num += len(snpsd.snps)
				removed_num += snpsd.onlyBinarySnps()
			print 'Removed %d non-binary SNPs out of %d SNPs' % (removed_num, total_num)
		return pd_indices_to_keep



	def get_ibs_kinship_matrix(self, debug_filter=1, num_dots=1000, snp_dtype='int8', dtype='single', type='binary'):
		"""
		
		"""
		print 'Starting kinship calculation, it prints %d dots.' % num_dots
		snps = self.getSnps(debug_filter)
		print 'Constructing a SNP array'
		snps_array = sp.array(snps, dtype=snp_dtype)
		print 'Transposing the SNP array'
		snps_array = snps_array.T
		num_lines = len(self.accessions)
		num_snps = float(len(snps))
		print 'Allocating K matrix'
		k_mat = sp.ones((num_lines, num_lines), dtype=dtype)
		num_comp = num_lines * (num_lines - 1) / 2
		comp_i = 0
		print 'Starting calculation'
		for i in range(num_lines):
			for j in range(i):
				comp_i += 1
				k_mat[i, j] = sp.sum(sp.absolute(snps_array[i] - snps_array[j])) / num_snps
				k_mat[j, i] = k_mat[i, j]
				if num_comp >= num_dots and (comp_i + 1) % (num_comp / num_dots) == 0: #Print dots
					sys.stdout.write('.')
					sys.stdout.flush()
		return k_mat


	def get_snp_cov_matrix(self, debug_filter=1, num_dots=100, dtype='single'):
		print 'Starting covariance calculation'
		#Normalizing
		norm_snps_array = self.get_normalized_snps(debug_filter=debug_filter, dtype=dtype)
		accession_means = sp.mean(norm_snps_array, 1)
		x = sp.mat(norm_snps_array.T - accession_means)
		cov_mat = x.T * x / (len(snps) - 1)
		print 'Finished calculating covariance matrix'
		return cov_mat



	def get_normalized_snps(self, debug_filter=1, dtype='single'):
		print 'Normalizing SNPs'
		snps = self.getSnps(debug_filter)
		snps_array = sp.array(snps)
		snps_array = snps_array.T
		#Normalizing
		norm_snps_array = (snps_array - sp.mean(snps_array, 0)) / sp.std(snps_array, 0)
		print 'Finished normalizing them'
		return norm_snps_array



	def convert_2_binary(self):
		"""
		Converts the underlying raw data format to a binary one, i.e. A,C,G,T,NA,etc. are converted to 0,1,-1
		"""
		if self.data_format == 'binary':
			import warnings
			warnings.warn("Data appears to be already in binary format!")
		else:
			snpsd_list = []
			for snpsd in self.snpsDataList:
				snpsd_list.append(snpsd.getSnpsData())
			self.snpsDataList = snpsd_list
			self.data_format = 'binary'
		self.missing_val = self.snpsDataList[0].missingVal


	def haplotize(self, snp_window=None, base_window=None):
		"""
		Converts the data-format to a haplotype numbering format
		"""
		print 'Haplotizing!'
		assert (snp_window or base_window), 'snp_window or base_window arguments are missing.'
		for i, sd in enumerate(self.snpsDataList):
			sd.haplotize(snp_window=snp_window, base_window=base_window)
			print i



	def impute_missing(self, reference_data):
		if not len(self.snpsDataList) == len(reference_data.snpsDataList):
			raise Exception("Reference data for imputation doesn't have equal chromosome number.")
		for i, snpsd in enumerate(self.snpsDataList):
			snpsd.inpute_data_region(reference_data.snpsDataList[i])




	def getSnps(self, random_fraction=None):
		snplist = []
		if random_fraction:
			import random
			for snpsd in self.snpsDataList:
				for snp in snpsd.snps:
					if random.random() < random_fraction:
						snplist.append(snp)
		else:
			for snpsd in self.snpsDataList:
				snplist += snpsd.snps
		return snplist


	def num_snps(self):
		num_snps = 0
		for snpsd in self.snpsDataList:
			num_snps += len(snpsd.snps)
		return num_snps


#	def get_snps(self, random_fraction=None, region=None):
#		snplist = []
#		if random_fraction:
#			import random
#			for snpsd in self.snpsDataList:
#				for snp in snpsd.snps:
#					if random.random() < random_fraction:
#						snplist.append(snp)
#		elif region:
#			chrom, start_pos, end_pos = region
#			snpsd = self.snpsDataList[chrom - 1]
#			ps_gen = izip(snpsd.positions, snpsd.snps)
#			p, s = ps_ge.next()
#			while p < start_pos:
#				p, s = ps_ge.next()
#			#Now found
#			while p < end_pos:
#				snplist.append(s)
#				pos_list = []
#				p, s = ps_ge.next()
#
#		else:
#			for snpsd in self.snpsDataList:
#				snplist += snpsd.snps
#		return snplist


	def get_snp_at(self, chromosome, position):
		"""
		Returns the SNP at the given position, if it exits.
		"""
		c_i = self.chromosomes.index(chromosome)
		sd = self.snpsDataList[c_i]
		i = 0
		while sd.positions[i] < position:
			i += 1
		if sd.positions[i] == position:
			print 'Found the SNP.'
			return sd.snps[i]
		else:
			print "Didn't find the SNP on chromosome %d, at position %d" % (chromosome, position)
			return None


	def getPositions(self):
		poslist = []
		for snpsd in self.snpsDataList:
			for pos in snpsd.positions:
				poslist.append(pos)
		return poslist


	def get_top_correlated_snp(self, snp, r2_threshold=0.5):
		sample_snp_chr_pos_marf = []
		sample_r2s = []
		print_r2s = []
		print_snp_chr_pos_marf = []
		for snpsd, chromosome in zip(self.snpsDataList, self.chromosomes):
			marfs = snpsd.get_mafs()['marfs']
			r2s = snpsd.calc_r2_with_snp(snp)
			for r2, sample_snp, sample_snp_pos, marf in zip(r2s, snpsd.snps, snpsd.positions, marfs):
				if r2 > r2_threshold:
					sample_r2s.append(r2)
					sample_snp_chr_pos_marf.append((sample_snp, chromosome, sample_snp_pos, marf))
				if r2 > 0.1:
					print_r2s.append(r2)
					print_snp_chr_pos_marf.append((sample_snp, chromosome, sample_snp_pos, marf))

		l = zip(print_r2s, print_snp_chr_pos_marf)
		l.sort()
		print l[-10:]
		return sample_snp_chr_pos_marf


	def get_all_snp_w_info(self):
		"""
		Returns the SNPs along with some info..
		"""
		sample_snp_chr_pos_marf = []
		for snpsd, chromosome in zip(self.snpsDataList, self.chromosomes):
			marfs = snpsd.get_mafs()['marfs']
			sample_snp_chr_pos_marf += zip(snpsd.snps, [chromosome] * len(snpsd.snps), snpsd.positions, marfs)
		return sample_snp_chr_pos_marf


	def getChrPosList(self):
		chr_pos_list = []
		for i in range(0, len(self.snpsDataList)):
			snpsd = self.snpsDataList[i]
			chr = i + 1
			for pos in snpsd.positions:
				chr_pos_list.append((chr, pos))
		return chr_pos_list

	def get_chr_list(self):
		chr_list = []
		for c, snpsd in izip(self.chromosomes, self.snpsDataList):
			chr_list.extend([c] * len(snpsd.positions))
		return chr_list


	def get_mafs(self):
		"""
		Returns the mafs and marfs as a dictionary.
		"""
		if not self.data_format:
			binary = True
		else:
			if  self.data_format == 'binary':
				binary = True
			else:
				binary = False

		maf_list = []
		marf_list = []
		for snpsd in self.snpsDataList:
			r = snpsd.get_mafs(binary=binary)
			maf_list.extend(r["mafs"])
			marf_list.extend(r["marfs"])
		print "Finished calculating MAFs."
		return {"mafs":maf_list, "marfs":marf_list}

	def getChrPosSNPList(self):
		chr_pos_snp_list = []
		for i in range(0, len(self.snpsDataList)):
			snpsd = self.snpsDataList[i]
			chr = i + 1
			for j in range(0, len(snpsd.positions)):
				pos = snpsd.positions[j]
				snp = snpsd.snps[j]
				chr_pos_snp_list.append((chr, pos, snp))
		return chr_pos_snp_list


	def get_region_pos_snp_dict(self, chromosome, start_pos=None, end_pos=None):
		"""
		Returns a dict containing a list of positions and snps.
		"""
		positions = []
		snps = []
		for snpsd in self.snpsDataList:
			if snpsd.chromosome == chromosome:
				break
		assert snpsd.chromosome == chromosome, "SNPs data with appropriate chromosomes wasn't found"
		i = 0
		while i < len(snpsd.positions) and snpsd.positions[i] < start_pos:
			i += 1
		if end_pos:
			while i < len(snpsd.positions) and snpsd.positions[i] < end_pos:
				positions.append(snpsd.positions[i])
				snps.append(snpsd.snps[i])
				i += 1
		else:
			positions = snpsd.positions[i:]
			snps = snpsd.snps[i:]
		return {'positions':positions, 'snps':snps}



	def get_pc(self, pc_num=1, random_fraction=0.1):
		"""
		Returns the pc_num'th principal components of the genotype 
		"""
		import random
		import rpy, util

		if not self.is_binary:
			print "Converting the snps data to binary format."
			self.convert_2_binary()

		snps = []
		for sd in self.snpsDataList:
			snps.extend(random.sample(sd.snps, int(random_fraction * len(sd.snps))))
		genotypes = map(list, zip(*snps))


		for genotype in genotypes:
			sd = util.calcSD(genotype)
			for j in range(len(genotype)):
				genotype[j] = genotype[j] / sd

		genotypes = sp.transpose(sp.array(genotypes))
		#print genotypes
		pc = rpy.r.princomp(genotypes)
		pc_sorted = zip(list(pc["scores"][pc_num - 1]), self.accessions)
		pc_sorted.sort()
		print pc_sorted
		pc = list(pc["scores"][pc_num - 1])
		self.pc = pc
		return pc




	def updateRegions(self, regionList):
		"""
		Deprecated 11/11/08 - Bjarni
		"""
		c_i = 0
		i = 0
		rl_i = 0 #region list index
		while c_i < len(self.chromosomes) and rl_i < len(regionList):
			region = regionList[rl_i]
			snpsd = self.snpsDataList[c_i]
			cp1 = (c_i + 1, snpsd.positions[i])
			cp_start = (region.chromosome, region.startPos)
			while cp1 < cp_start:
				if i < len(snpsd.positions):
					i += 1
				else:
					c_i += 1
					snpsd = self.snpsDataList[c_i]
					i = 0
				cp1 = (c_i + 1, snpsd.positions[i])
			cp_end = (region.chromosome, region.endPos)
			while cp1 <= cp_end:
				"""Update current region!"""
				region.snps.append(snpsd.snps[i])
				region.snps_indices.append((c_i, i))

				i += 1
				if i < len(snpsd.positions):
					cp1 = (c_i + 1, snpsd.positions[i])
				else:
					c_i += 1
					i = 0
					break

			rl_i += 1

	def filter_na_snps(self, max_na_rate=0.0):
		for snpsd in self.snpsDataList:
			snpsd.filter_na_snps(max_na_rate=max_na_rate)

	def remove_snps_indices(self, indices_to_remove):
		remove_indices = [[] for i in range(len(self.snpsDataList))]
		indices_to_remove.sort()
		offset = 0
		i = 0
		for snpsd_i, snpsd in enumerate(self.snpsDataList):
			max_i = offset + len(snpsd.snps)
			i2r = indices_to_remove[i]
			while i < len(indices_to_remove) and i2r < max_i:
					remove_indices[snpsd_i].append(i2r - offset)
					i += 1
					i2r = indices_to_remove[i]
			offset = max_i
		for i, snpsd in enumerate(self.snpsDataList):
			snpsd.remove_snps(remove_indices[i])


	def filter_na_accessions(self, max_na_rate=0.2, verbose=False):
		accessions_na_counts = [0.0 for acc in self.accessions]
		total_snps = 0
		for snpsd in self.snpsDataList:
			for i, c in enumerate(snpsd.countMissingSnps()[1]):
				accessions_na_counts[i] += c
			total_snps += float(len(snpsd.snps))
		accessions_na_counts = [accessions_na_counts[i] / total_snps for i in range(len(accessions_na_counts))]
		#print accessions_na_counts
		acc_to_keep = []
		for i, na_rate in enumerate(accessions_na_counts):
			if na_rate <= max_na_rate:
				acc_to_keep.append(self.accessions[i])
		#print len(acc_to_keep)
		self.filter_accessions(acc_to_keep)



	def filter_maf_snps(self, maf, maf_ub=1):
		for snpsd in self.snpsDataList:
			snpsd.snpsFilterMAF([maf, maf_ub])


	def filter_mac_snps(self, mac_threshold=15):
		for snpsd in self.snpsDataList:
			snpsd.filter_mac(mac_threshold)


	def filter_monomorphic_snps(self):
		for snpsd in self.snpsDataList:
			snpsd.filterMonoMorphicSnps()


	def filter_accessions(self, accessions_to_keep, use_accession_names=False):
		assert len(accessions_to_keep) != 0, "Can't remove all ecotypes."
		if use_accession_names:
			import phenotypeData as pd
			ad = pd._getAccessionToEcotypeIdDict_(accessions_to_keep)
			ecotypes = []
			for acc in accessions_to_keep:
				ecotypes.append(str(ad[acc]))
		else:
			ecotypes = accessions_to_keep
		num_accessions = len(self.accessions)
		acc_indices_to_keep = []
		for et in ecotypes:
			try:
				i = self.accessions.index(et)
				acc_indices_to_keep.append(i)
			except:
				continue
		#pdb.set_trace()
		#acc_indices_to_keep.sort()
		self.filter_accessions_indices(acc_indices_to_keep)


	def filter_accessions_indices(self, acc_indices_to_keep):
		num_accessions = len(self.accessions)
		for i, snpsd in enumerate(self.snpsDataList):
			#print i, len(snpsd.accessions)
			snpsd.removeAccessionIndices(acc_indices_to_keep)
		self.accessions = self.snpsDataList[0].accessions
		self.array_ids = self.snpsDataList[0].arrayIds
		print "Removed %d accessions, leaving %d in total." % (num_accessions - len(acc_indices_to_keep), len(acc_indices_to_keep))


	def filter_for_countries(self, country_codes, complement=False):
		import phenotypeData as pd
		ei_dict = pd._getEcotypeIdInfoDict_()
		acc_indices_to_keep = []
		if complement:
			for i, ei in enumerate(self.accessions):
				if ei_dict[int(ei)][4] not in country_codes:
					acc_indices_to_keep.append(i)
		else:
			 for i, ei in enumerate(self.accessions):
				if ei_dict[int(ei)][4] in country_codes:
					acc_indices_to_keep.append(i)
		self.filter_accessions_indices(acc_indices_to_keep)

		for ei in self.accessions:
			print ei, ei_dict[int(ei)]


	def sample_snps(self, random_fraction):
		"""
		Samples a random fraction of the SNPs.
		"""

		for snpsd in self.snpsDataList:
			snpsd.sample_snps(random_fraction)


	def get_region_snpsd(self, chr, start_pos=None, end_pos=None):
		"""
		Modifies original object, beware!
		"""
		c_i = self.chromosomes.index(chr)
		print chr, start_pos, end_pos
		snpsd = self.snpsDataList[c_i]
		if start_pos != None:
			new_snps = []
			new_positions = []
			i = 0
			while i < len(snpsd.positions) - 1 and snpsd.positions[i] < start_pos:
				i += 1
			while i < len(snpsd.positions) - 1 and snpsd.positions[i] < end_pos:
				new_snps.append(snpsd.snps[i])
				new_positions.append(snpsd.positions[i])
				i += 1
			snpsd.positions = new_positions
			snpsd.snps = new_snps
		return snpsd

	def merge_snps_data(self, sd):
		"""
		Merges data using _SnpsData_.merge_data
		"""
		if self.chromosomes != sd.chromosomes:
			raise Exception
		else:
			self.new_snps_data_list = []
			for sd1, sd2, chromosome in zip(self.snpsDataList, sd.snpsDataList, self.chromosomes):
				print "Merging data on chromosome %s." % (str(chromosome))
				sd1.merge_data(sd2)
				self.new_snps_data_list.append(sd1)
			self.accessions = self.new_snps_data_list[0].accessions
			self.snpsDataList = self.new_snps_data_list




def readSNPsDataSetFile(delim=","):
	"""
	Read data file and return a SNPsDataSet object.
	"""
	pass

def readSNPsDataSetAccessions(datafile, delim=","):
	f = open(datafile, 'r')
	line = f.readline().split(delim)
	arrayIDs = []
	if line[0] != "Chromosome":
		for aID in line[2:]:
			arrayIDs.append(aID.strip())
		line = f.readline().split(delim)
	f.close()
	accessions = []
	for acc in line[2:]:
		accessions.append(acc.strip())
	return (accessions, arrayIDs)





def writeRawSnpsDatasToFile(filename, snpsds, chromosomes=[1, 2, 3, 4, 5], deliminator=", ", missingVal="NA", accDecoder=None, withArrayIds=False, callProbFile=None):
	"""
	2008-05-17
		for callProbFile, modify it to output a row immediately. no memory hoarding, untested.
	2008-05-17
		modify it to output a row immediately. no memory hoarding
		fix a slight bug. arrayIds is outputted always using ', ' as delimiter
	Writes data to a file. 
	"""
	import sys, csv
	sys.stderr.write("Writing data to file: %s ..." % filename)
	numSnps = 0
	for i in range(0, len(chromosomes)):
		numSnps += len(snpsds[i].positions)

	accessions = snpsds[0].accessions
	for i in range(1, len(chromosomes)):
		if accessions != snpsds[i].accessions:
			raise Exception("Accessions are different between SNPs datas")


	decoder = RawDecoder()
	decoder[snpsds[0].missingVal] = missingVal

	print "NumSnps: " + str(numSnps) + ", NumAcc: " + str(len(accessions)) + "\n"
	#writer = csv.writer(open(filename, 'w'), delimiter=deliminator)
	if withArrayIds:
		#writer.writerow(['-', '-']+snpsds[0].arrayIds)
		outStr = deliminator.join(['-', '-']) + deliminator + deliminator.join(snpsds[0].arrayIds) + "\n"
	else:
		outStr = ""
	fieldStrings = ["Chromosome", "Positions"]
	if accDecoder:
		for acc in snpsds[i].accessions:
			fieldStrings.append(str(accDecoder[acc]))
	else:
		for acc in snpsds[i].accessions:
			fieldStrings.append(str(acc))
	#writer.writerow(fieldStrings)
	f = open(filename, 'w')
	outStr += deliminator.join(fieldStrings) + "\n"
	f.write(outStr)
	import util  #Used to convert val list to stringlist.
	for i in range(0, len(chromosomes)):
		for j in range(0, len(snpsds[i].positions)):
			outStr = str(chromosomes[i]) + deliminator + str(snpsds[i].positions[j]) + deliminator
			snp = util.valListToStrList(snpsds[i].snps[j])
			outStr += deliminator.join(snp) + "\n"
			f.write(outStr)
			f.flush()
	f.close()

	if callProbFile:
		outStr = ""
		if withArrayIds:
			outStr = deliminator.join(["-", "-"] + snpsds[0].arrayIds) + "\n"
		f = open(callProbFile, "w")
		outStr += deliminator.join(fieldStrings) + "\n"
		f.write(outStr)
		f.flush()
		for i in range(0, len(chromosomes)):
			snpsd = snpsds[i]
			for j in range(0, len(snpsd.positions)):
				outStr = str(chromosomes[i]) + deliminator + str(snpsd.positions[j])
				for k in range(0, len(snpsd.accessions)):
					outStr += deliminator + str(snpsd.callProbabilities[j][k])
				outStr += "\n"
				f.write(outStr)
				f.flush()
			del snpsd

		f.close()

	sys.stderr.write("Done.\n")

def getMAF(snp, alphabet=[0, 1]):
	counts = []
	for letter in alphabet:
		counts.append(snp.count(letter))
	maf = min(counts)
	marf = maf / float(sum(counts))
	return (marf, maf)

def estimateRecomb(snpsdList, baseNum, filterProb, id):
	rho = 0
	npairs = 0
	for i in range(0, len(snpsdList)):
		snpsd = snpsdList[i]
		tmp1 = "tmp" + id + "1"
		tmp2 = "tmp" + id + "2"
		(rho2, npairs2) = snpsd.estimateRecomb(baseNum, filterProb, tmp1, tmp2)
		rho = rho + rho2 * npairs2
		npairs = npairs + npairs2
	rho = rho / float(npairs)
	print "rho: " + str(rho) + ", npairs: " + str(npairs)
	return rho

def D(freqs):
	""" Returns the D' LD measure """
	p1 = freqs[1] + freqs[3]  #Always positive (if no trivial SNPs are allowed).
	p2 = freqs[2] + freqs[3]
	Dmax = min(p1 * (1 - p2), p2 * (1 - p1))
	Dmin = -min(p1 * p2, (1 - p2) * (1 - p1))
	D = freqs[3] - p1 * p2
	if D >= 0.0:
		return D / Dmax
	else:
		return D / Dmin

def r2(freqs):
	f1 = freqs[1] + freqs[3]
	f2 = freqs[2] + freqs[3]
	D = freqs[3] - f1 * f2
	divisor = f1 * f2 * (1 - f1) * (1 - f2)
	if divisor != 0:
		return D * D / divisor
	else:
		return - 1



def r2listAll(data, windowSize, nbins=10):
	sums = {1:[0.0] * nbins, 2:[0.0] * nbins, 3:[0.0] * nbins, 4:[0.0] * nbins, 5:[0.0] * nbins, 6:[0.0] * nbins, 7:[0.0] * nbins, 8:[0.0] * nbins, "mean":[0.0] * nbins}
	counts = {1:[0] * nbins, 2:[0] * nbins, 3:[0] * nbins, 4:[0] * nbins, 5:[0] * nbins, 6:[0] * nbins, 7:[0] * nbins, 8:[0] * nbins, "tot":[0] * nbins}
	for snpsd in data:
		fsnpd = snpsd.snpsFilter()
		for k in [1, 2, 3, 4, 5, 6, 7, 8]:
			freqs = 	fsnpd[k].calcFreqs(windowSize)
			for i in xrange(0, len(freqs)):
				for j in xrange(0, len(freqs[i])):
					r = r2(freqs[i][j])
					if r != 	 -1:
						bin = int((fsnpd[k].positions[j + i + 1] - fsnpd[k].positions[i]) * nbins / (windowSize + 0.01))
						#print fsnpd[k].positions[j+i+1], fsnpd[k].positions[i]
						#print bin						
						counts[k][bin] = counts[k][bin] + 1
						sums[k][bin] = sums[k][bin] + 	r

	for k in [1, 2, 3, 4, 5, 6, 7, 8]:
		for i in xrange(0, nbins):
			if counts[k][i] != 0:
				sums["mean"][i]	 = sums["mean"][i] + sums[k][i]
				sums[k][i] = float(sums[k][i]) / float(counts[k][i])
				counts["tot"][i] = counts["tot"][i] + counts[k][i]
	for i in xrange(0, nbins):
		if counts["tot"][i] != 0:
			sums["mean"][i] = float(sums["mean"][i]) / float(counts["tot"][i])

	return (sums)



def DlistAll(snpsdlist, windowSize, nbins=10):
	sums = {1:[0.0] * nbins, 2:[0.0] * nbins, 3:[0.0] * nbins, 4:[0.0] * nbins, 5:[0.0] * nbins, 6:[0.0] * nbins, 7:[0.0] * nbins, 8:[0.0] * nbins, "mean":[0.0] * nbins}
	counts = {1:[0] * nbins, 2:[0] * nbins, 3:[0] * nbins, 4:[0] * nbins, 5:[0] * nbins, 6:[0] * nbins, 7:[0] * nbins, 8:[0] * nbins, "tot":[0] * nbins}
	for snpsd in snpsdlist:
		fsnpd = snpsd.snpsFilter()
		for k in [1, 2, 3, 4, 5, 6, 7, 8]:
			freqs = 	fsnpd[k].calcFreqs(windowSize)
			for i in xrange(0, len(freqs)):
				for j in xrange(0, len(freqs[i])):
					r = D(freqs[i][j])
					if r != 	 -1:
						bin = int((fsnpd[k].positions[j + i + 1] - fsnpd[k].positions[i]) * nbins / (windowSize + 0.01))
						counts[k][bin] = counts[k][bin] + 1
						sums[k][bin] = sums[k][bin] + r

	for k in [1, 2, 3, 4, 5, 6, 7, 8]:
		for i in xrange(0, nbins):
			if counts[k][i]	 != 0:
				sums["mean"][i]	 = sums["mean"][i] + sums[k][i]
				sums[k][i] = float(sums[k][i]) / float(counts[k][i])
				counts["tot"][i] = counts["tot"][i] + counts[k][i]
	for i in xrange(0, nbins):
		if counts["tot"][i] != 0:
			sums["mean"][i]	 = float(sums["mean"][i]) / float(counts["tot"][i])

	return (sums)


def DlistAll2(snpsdlist, windowSize, nbins=10):
	sums = {0:[0.0] * nbins, 1:[0.0] * nbins, 2:[0.0] * nbins, 3:[0.0] * nbins, 4:[0.0] * nbins, "mean":[0.0] * nbins}
	counts = {0:[0] * nbins, 1:[0] * nbins, 2:[0] * nbins, 3:[0] * nbins, 4:[0] * nbins, "tot":[0] * nbins}
	for snpsd in snpsdlist:
		fsnpdata = snpsd.snpsFilter()
		for k in [0, 1, 2, 3, 4]:
			fsnpsd1 = fsnpdata[k]
			fsnpsd2 = fsnpdata[9 - k]
			for fsnpd in [fsnpsd1, fsnpsd2]:
				freqs = 	fsnpd.calcFreqs(windowSize)
				for i in xrange(0, len(freqs)):
					for j in xrange(0, len(freqs[i])):
						r = D(freqs[i][j])
						if r != 	 -1:
							bin = int((fsnpd.positions[j + i + 1] - fsnpd.positions[i]) * nbins / (windowSize + 0.01))
							counts[k][bin] = counts[k][bin] + 1
							sums[k][bin] = sums[k][bin] + r


	for k in [0, 1, 2, 3, 4]:
		for i in xrange(0, nbins):
			if counts[k][i] != 0:
				sums["mean"][i]	 = sums["mean"][i] + sums[k][i]
				sums[k][i] = float(sums[k][i]) / float(counts[k][i])
				counts["tot"][i] = counts["tot"][i] + counts[k][i]
	for i in xrange(0, nbins):
		if counts["tot"][i] != 0:
			sums["mean"][i]	 = float(sums["mean"][i]) / float(counts["tot"][i])

	return (sums)

def r2listAll2(snpsdlist, windowSize, nbins=10):
	sums = {0:[0.0] * nbins, 1:[0.0] * nbins, 2:[0.0] * nbins, 3:[0.0] * nbins, 4:[0.0] * nbins, "mean":[0.0] * nbins}
	counts = {0:[0] * nbins, 1:[0] * nbins, 2:[0] * nbins, 3:[0] * nbins, 4:[0] * nbins, "tot":[0] * nbins}
	for snpsd in snpsdlist:
		fsnpdata = snpsd.snpsFilter()
		for k in [0, 1, 2, 3, 4]:
			fsnpsd1 = fsnpdata[k]
			fsnpsd2 = fsnpdata[9 - k]
			for fsnpd in [fsnpsd1, fsnpsd2]:
				freqs = 	fsnpd.calcFreqs(windowSize)
				for i in xrange(0, len(freqs)):
					for j in xrange(0, len(freqs[i])):
						r = r2(freqs[i][j])
						if r != 	 -1:
							bin = int((fsnpd.positions[j + i + 1] - fsnpd.positions[i]) * nbins / (windowSize + 0.01))
							counts[k][bin] = counts[k][bin] + 1
							sums[k][bin] = sums[k][bin] + r


	for k in [0, 1, 2, 3, 4]:
		for i in xrange(0, nbins):
			if counts[k][i] != 0:
				sums["mean"][i]	 = sums["mean"][i] + sums[k][i]
				sums[k][i] = float(sums[k][i]) / float(counts[k][i])
				counts["tot"][i] = counts["tot"][i] + counts[k][i]
	for i in xrange(0, nbins):
		if counts["tot"][i] != 0:
			sums["mean"][i]	 = float(sums["mean"][i]) / float(counts["tot"][i])

	return (sums)





#def writeRFormat(plotData, list = [0,1,2,3,4,"mean"],windowSize=20000,ylab="",lab=""):
#	delta = float(windowSize)/float(len(plotData[list[0]]))
#	st = "xv <- seq("+str(delta/2.0)+","+str(windowSize)+","+str(delta)+")\n"
#	maxlim = 0.0
#	minlim = 1.0
#	for k in list:
#		data = plotData[k]
#		st = st+"d"+str(k)+" <- c("
#		for i in range(0,len(data)-1):
#			st = st + str(data[i])+", "
#			maxlim = max(maxlim,data[i])
#			minlim = min(minlim,data[i])
#		st = st+str(data[len(data)-1])+")\n"
#		maxlim = max(maxlim,data[len(data)-1])
#		minlim = min(minlim,data[len(data)-1])
#	i = 1
#	for k in list[0:len(list)-1]:
#		st = st+"plot(y=d"+str(k)+", x=xv, ylim = c("+str(minlim)+","+str(maxlim)+"), col = "+str(i)+', type = "b",ylab="",xlab="")\npar(new=T)\n'
#		i = i + 1	
#	st = st+"plot(y=d"+str(list[len(list)-1])+", x=xv, ylim = c("+str(minlim)+","+str(maxlim)+'), type = "b",ylab="'+ylab+'",xlab="Bases", col = '+str(i)+', main="'+lab+'")\n'
#	return st

def write_accessions_info_file(filename, file_250k_data="/Users/bjarnivilhjalmsson/Projects/Data/250k/250K_t43_192.csv"):
	import dataParsers
	snps_accessions = dataParsers.parseCSVDataAccessions("/Users/bjarnivilhjalmsson/Projects/Data/250k/250K_data_t43_081009.csv")
	print len(snps_accessions)
	import phenotypeData
	e_dict = phenotypeData._getEcotypeIdToStockParentDict_()
	f = open(filename, "w")
	f.write("ecotype_id, accession_name, stock_parent\n")
	for e_id in snps_accessions:
		try:
			(acc, sp) = e_dict[int(e_id)]
			acc = unicode(acc, "latin-1")
			#print acc
			f.write(e_id + ", " + acc + ", " + sp + "\n")
		except Exception, err_str:
			print err_str, e_id
			#acc = unicode(acc,"utf-16")
			f.write(e_id + ", " + repr(acc) + ", " + sp + "\n")

	f.close()


def get_FRI_data_slice():
	import phenotypeData as pd
	import dataParsers as dp
	phend = pd.readPhenotypeFile('/Users/bjarnivilhjalmsson/Projects/Data/phenotypes/phen_all_051710.tsv')
	phend.removePhenotypeIDs([43])
	sd = dp.parse_snp_data_region('/Users/bjarnivilhjalmsson/Projects/Data/250k/250K_t52.csv', 4, 0, 1600000, format=0)
	sd = dp.parse_snp_data_region('/Users/bjarnivilhjalmsson/Projects/Data/250k/250K_t52.csv', 4, 0, 1600000, format=0)
	coordinateSnpsAndPhenotypeData(phend, 43, [sd])
	phen_vals = map(str, phend.getPhenVals(43))
	f = open('/tmp/fri_region_chr4_100kb-700kb.csv', 'w')
	f.write('chromosome,position,' + ','.join(sd.accessions) + '\n')
	f.write(',,' + ','.join(phen_vals) + '\n')
	for pos, snp in zip(sd.positions, sd.snps):
		snp = map(str, snp)
		f.write('4,' + str(pos) + ',' + ','.join(snp) + '\n')
	f.close()


def get_AW_common_dataset():
	import phenotypeData as pd
	filename = "/Users/bjarnivilhjalmsson/Projects/Data/phenotypes/phen_wilzcek_wo_OF_NS06_060210.tsv"
	phed = pd.readPhenotypeFile(filename)
	sd_t54 = dataParsers.parse_snp_data('/Users/bjarnivilhjalmsson/Projects/Data/250k/250K_t54.csv')#,filter=0.001)	
	sd_t54.filter_accessions(phed.accessions)
	cmid = sd_t54.add_to_db('common_wilczek_54', parent_id=54 , imputed=1 , unique_ecotype=1)
	sd_t54.writeToFile('/Users/bjarnivilhjalmsson/Projects/Data/250k/250K_t' + str(cmid) + '.csv', withArrayIds=True)


def write_out_01_dataset():
	sd_t75 = dataParsers.parse_snp_data(env.env['data_dir'] + '250K_t75.csv', format='binary')
	file_name = '/Users/bjarnivilhjalmsson/Projects/Data/250k/250K_t75.csv.binary'
	sd_t54.writeToFile(file_name)

def get_JB_datasets():
	pass


def get_call_method_dataset_file(call_method_id, binary_format=False):
	if binary_format:
		return env.env['data_dir'] + '250K_t' + str(call_method_id) + '.csv.binary'
	else:
		return env.env['data_dir'] + '250K_t' + str(call_method_id) + '.csv'

def get_call_method_kinship_file(call_method_id):
	return env.env['data_dir'] + 'kinship_matrix_cm' + str(call_method_id) + '.pickled'


def get_bergelssons_region_datasets():
	import phenotypeData as pd
	first_192 = pd._getFirst192Ecotypes_()
	print first_192
	print len(first_192)
#	sd_t54 = dataParsers.parse_snp_data('/Users/bjarnivilhjalmsson/Projects/Data/250k/250K_t54.csv')	
#	sd_t54.filter_for_countries(['UK'],complement=True)
#	cm_id = sd_t54.add_to_db('cm54_only_non_UK',method_description='',data_description='',comment='')
#	sd_t54.writeToFile('/Users/bjarnivilhjalmsson/Projects/Data/250k/250K_t'+str(cm_id)+'.csv')

if __name__ == "__main__":
	import dataParsers
#	d2010_file = "/Users/bjarnivilhjalmsson/Projects/Data/2010/2010_073009.csv"
#	d2010_sd = dataParsers.parse_snp_data(d2010_file,id="2010_data")
#	d250k_file = "/Users/bjarnivilhjalmsson/Projects/Data/250k/250K_t43_192.csv"
#	d250k_sd = dataParsers.parse_snp_data(d250k_file)
#	d2010_sd.impute_missing(d250k_sd)
#	d2010_sd.writeToFile("/tmp/test.csv")
#       get_AW_common_dataset()
#	write_out_01_dataset()
	get_bergelssons_region_datasets()





