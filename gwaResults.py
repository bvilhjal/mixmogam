"""
Contains classes to handle results of GWAS

020210
- TODO:
	- Result object should be coupled to a marker data, in a simple manner.  
	  E.g. with a pointer, and/or a list marking what markers are in the result data. 

"""


import pdb, gc
import dbutils
import csv
import math
import itertools as it
try:
	import scipy as sp
	import scipy.stats as st
except Exception, err_str:
	print 'scipy is missing:', err_str
import random
import sys
import os
import cPickle
import env
import tair_converter as tc
#A dictionary for loaded results.. to avoid reloading. 
#Use carefully to avoid memory leaks!



class ResultType(object):
	def __init__(self, resultType=None, fileType=None, datasetName=None, resultDir=None, mafCutoff=0, logTransform=None, name=None):
		self.resultType = resultType
		self.fileType = fileType
		self.datasetName = datasetName
		self.resultDir = resultDir
		self.mafCutoff = mafCutoff
		self.logTransform = logTransform
		if self.logTransform == None:
			self.logTransform = self.fileType == ".pvals"
		self.name = name
		if not name and datasetName:
			self.name = resultType + "_" + datasetName
		if resultDir and datasetName:
			self.filenamePrefix = resultDir + resultType + "_" + datasetName + "_"


	def getFileName(self, phed, phenotypeID, secondRun=False):
		print phenotypeID
		if secondRun:
			filename = self.filenamePrefix + str(phed.getPhenotypeName(phenotypeID)) + ".sr" + self.fileType
		else:
			filename = self.filenamePrefix + str(phed.getPhenotypeName(phenotypeID)) + self.fileType
		return filename

	def __str__(self):
		return self.resultType + "_" + self.datasetName


class Region(object):

	def __init__(self, chromosome, startPos, endPos, snps=None, snpsd_indices=None, snpsInfo=None, name=""):
		self.name = name
		self.chromosome = chromosome
		self.startPos = startPos
		self.endPos = endPos
		self.size = endPos - startPos
		self.snps = []
		if snps:
			self.snps = snps
		self.snpsd_indices = []  #Indicies in other data structures.
		if snpsd_indices:
			self.snpsd_indices = snpsd_indices
		self.snpsInfo = {}       #E.g. a dictionary of information about the snps.  
		if snpsInfo:
			self.snpsInfo = snpsInfo

	def __cmp__(self, other):
		return cmp((self.chromosome, self.startPos), (other.chromosome, other.startPos))

	def __str__(self):
		return "Chr.:" + str(self.chromosome) + ", start pos.:" + str(self.startPos) + ", end pos.:" + str(self.endPos)

	def get_chr_pos_str(self):
		return str(self.chromosome) + "_" + str(self.startPos) + "_" + str(self.endPos)

	def overlapping(self, region):
		if self.chromosome != region.chromosome:
			return False
		else:
			if self.startPos < region.startPos:
				return self.endPos >= region.startPos
			else:
				return region.endPos >= self.startPos

	def merge(self, region):
		if not self.overlapping(region):
			raise Exception
		new_start = min(self.startPos, region.startPos)
		new_stop = max(self.endPos, region.endPos)
		self.startPos = new_start
		self.endPos = new_stop
		self.size = self.endPos - self.startPos


#fri_region_small = Region(4, 220000, 320000, name="FRI_small")
#flc_region_small = Region(5, 3130000, 3230000, name="FLC_small")
#
#fri_region = Region(4, 120000, 550000, name="FRI_large")
#flc_region = Region(5, 2100000, 3300000, name="FLC_large")
#
#fri_flc_regions = [fri_region_small, fri_region, flc_region_small, flc_region]

def getRegions(regionSet, window=[25000, 25000]):
	"""
	Converts a set of crh_pos into a list of regions objects.
	"""

	res_ls = list(regionSet)
	res_ls.sort()

	oldPos = 0
	countRegions = 0
	chr = -1
	regions = []
	curPosList = []
	for i in range(0, len(res_ls)):
		pos = res_ls[i][1]
		if chr != res_ls[i][0]:
			if len(curPosList):
				regions.append(curPosList)
			curPosList = []
			countRegions += 1
			chr = res_ls[i][0]
		elif pos - oldPos > sum(window):
			if len(curPosList):
				regions.append(curPosList)
			curPosList = []
			countRegions += 1

		curPosList.append((chr, pos))
		oldPos = pos

	print countRegions, len(regions)

	regionList = []
	for region in regions:
		chr = region[0][0]
		positions = []
		for (chr, pos) in region:
			positions.append(pos)
		regionList.append(Region(chr, min(positions) - window[0], max(positions) + window[1]))
	regionList.sort()
	return regionList


class Result(object):
	"""
	Contains information on the result.  (The old object is renamed Result_old, for now..)  Poised to cause problems?
	"""
	pickle_attributes = ['snp_results', 'keys', 'orders', 'ranks', 'phen_id', 'result_type', 'name',
				'accessions', 'chromosome_ends']

	def __init__(self, result_file=None, snp_results=None, scores=None, snps_data=None, accessions=None, name=None,
		     result_type=None, phen_id=None, positions=None, chromosomes=None, mafs=None, macs=None,
		     snps=None, auto_pickling_on=True, keys=None, **snp_results_info):
		"""
		A simple init function.
		snps_data is assumed to match the results, (in term of positions, etc)
		
		accessions (when used) should fit the results!
		"""

		self.snp_results = {}  #Main data container
		if snp_results:
			self.snp_results = snp_results
		self.phen_id = phen_id
		self.result_type = result_type
		self.name = name
		if keys:
			self.keys = keys
		else:
			self.keys = []
		self.accessions = accessions
		self.chromosome_ends = []
		self.orders = None
		self.ranks = None
		self.auto_pickling_on = auto_pickling_on


		if result_file:
			pickle_file = result_file + '.pickled'
			if self.auto_pickling_on and os.path.isfile(pickle_file):
				try:
					with open(pickle_file) as f:
						d = cPickle.load(f)
					for attr in self.pickle_attributes:
						setattr(self, attr, d[attr])
					pickle_failed = False
				except Exception, err_str:
					print 'Loading pickle file failed:', pickle_file
					self._load_result_(result_file)
					pickle_failed = True

			else:
				self._load_result_(result_file)
		else:

			if chromosomes != None:
				self.keys.append('chromosomes')
				self.snp_results['chromosomes']	 = chromosomes
			if positions != None:
				self.keys.append('positions')
				self.snp_results['positions'] = positions
			if scores != None:
				self.keys.append('scores')
				self.snp_results['scores'] = scores
			if mafs != None:
				self.keys.append('mafs')
				self.snp_results['mafs'] = mafs
			if macs != None:
				self.keys.append('macs')
				self.snp_results['macs'] = macs

		if snps_data:
			self._load_snps_(snps_data)
			if 'mafs' not in self.snp_results:
				self._load_mafs_(snps_data)
		else:
			if snps:
				self.keys.append('snps')
				self.snp_results['snps'] = snps



		#Adding various info...
		if snp_results_info:
			for info in snp_results_info:
				self.keys.append(info)
				self.snp_results[info] = snp_results_info[info]

		if result_file and self.auto_pickling_on:
			if not os.path.isfile(pickle_file) or pickle_failed:
				d = {}
				for attr in self.pickle_attributes:
					d[attr] = getattr(self, attr)
				with open(pickle_file, 'wb') as f:
					cPickle.dump(d, f)





	def _load_result_(self, result_file, try_delims=[",", "\t", " ", ]):

		with open(result_file, "r") as f:
			print 'Loading results..'
			line = f.next()
			for delim in try_delims:
				if len(line.split(delim)) > 1:
					break
			else:
				raise Exception("Appropriate delimiter wasn't found.")

			header = line.split(delim)
			for key in header:  #Handling old data formats..
				key = key.lower()
				if key == 'chromosome':
					key = 'chromosomes'
				if key == 'position':
					key = 'positions'
				if key == 'score':
					key = 'scores'
				if key == 'marf':
					key = 'mafs'
				if key == 'maf':
					key = 'macs'
				self.keys.append(key)

			for key in self.keys:
				self.snp_results[key] = []

			chrom = 0
			chrom_splits = []
			for i, line in enumerate(f):
				l = map(float, map(str.strip, line.split(delim)))
				for v, key in it.izip(l, self.keys):
					self.snp_results[key].append(v)
				if int(l[0]) != chrom:
					chrom_splits.append(i)
					chrom = int(l[0])

			for si in chrom_splits[1:]:
				self.chromosome_ends.append(self.snp_results['positions'][si - 1])
			self.chromosome_ends.append(self.snp_results['positions'][-1])

			self._rank_scores_()


	def _load_snps_(self, snps_data):
		"""
		This only loads SNPs...
		"""
		snps = []
		positions = []
		chromosomes = []
		for i, snpsd in enumerate(snps_data.snpsDataList):
			snps.extend(snpsd.snps)
			pl = snpsd.positions
			positions.extend(pl)
			chromosomes.extend([snps_data.chromosomes[i]] * len(pl))
			self.chromosome_ends.append(snpsd.positions[-1])
		self.snp_results['snps'] = snps
		self.snp_results['positions'] = positions
		self.snp_results['chromosomes'] = chromosomes
		self.keys.append('chromosomes')
		self.keys.append('positions')
		self.keys.append('snps')


	def _load_mafs_(self, snps_data):
		snps = self.snp_results['snps']
		if snps_data.data_format == 'float':
			self.snp_results['mafs'] = [1] * len(snps)
			self.snp_results['macs'] = [1] * len(snps)
		else:
			maf_d = snps_data.get_mafs()
			self.snp_results['mafs'] = maf_d['marfs']
			self.snp_results['macs'] = maf_d['mafs']
		self.keys.append('mafs')
		self.keys.append('macs')

	def _rank_scores_(self):
		"""
		Generates two data structures:
		self.orders (SNP indices ordered by rank)
		self.ranks (ranks of the SNPs)
		"""
		scores = self.snp_results['scores']
		rank_ls = zip(scores, range(len(scores)))
		rank_ls.sort(reverse=True)
		self.orders = []
		for j in range(len(rank_ls)):
			(s, i) = rank_ls[j]
			self.orders.append(i)
			rank_ls[j] = (i, j)

		rank_ls.sort()
		self.ranks = []
		for (i, j) in rank_ls:
			self.ranks.append(j + 1)


	def _sort_by_chr_pos_(self):
		res_ls = zip(self.snp_results['chromosomes'], self.snp_results['positions'],
			range(len(self.snp_results['positions'])))
		res_ls.sort()
		orders = list(zip(*res_ls))[2]

		snp_results = {}
		for key in self.keys:
			snp_results[key] = []
			for i in orders:
				snp_results[key].append(self.snp_results[key][i])
		self.snp_results = snp_results


	def candidate_gene_enrichments(self, cgl=None, cgl_file=None, pval_thresholds=[0.01], gene_radius=20000,
				methods=['chi_square'], num_perm=500, file_prefix=None,
				obs_genes_file=None, early_stop_threshold=25, all_genes=None, cand_gene_indices=None):
		"""
		Performs CGR analysis on this results object.
		
		cgl is a list of genes.
		cgl_file is a file with a list of genes.
		
		method: 
			chi_square - statistics
			multinomial - statistics
			gene_perm - permute genes.
			snps_perm - permute SNPs
		"""
		import pylab
		import analyze_gene_enrichment as genr

		if not 'chi_square' in methods:
			methods.append('chi_square')

		chrom_ends = self.get_chromosome_ends()
		print 'chrom_ends', chrom_ends

		#Parse cgl file
		if cgl_file:
			cgl, cg_tair_ids = load_cand_genes_file(cgl_file)
		for cg in cgl:
			print str(cg)


		if not all_genes:
			#Load genes from DB.
			print 'Fetching all genes'
			all_genes = get_gene_list(include_intron_exons=False, verbose=False)
			print 'Fetched %d genes.' % len(all_genes)
			num_genes = len(all_genes)

		if not cand_gene_indices:
			#Pre-process cgl.
			cand_gene_indices = []
			for i, g in enumerate(all_genes):
				for cg in cgl:
					if g.dbRef == cg.dbRef:
						#print g.dbRef, cg.dbRef
						cand_gene_indices.append(i)
						break
			num_cand_genes = len(cand_gene_indices)
			#print num_cand_genes, cand_gene_indices


		method_res_dict = {}
		for m in methods:
			method_res_dict[m] = {'statistics':[], 'pvals':[]}

		if obs_genes_file:
			obs_gene_str = ''

		pval_thresholds.sort()
		last_thres = 1.0
		log_enrichments = []
		for pval_threshold in reversed(pval_thresholds):
			print 'Using p-value threshold % f' % pval_threshold
			thres = pval_threshold / last_thres
			print 'Using corrected threshold % f' % thres
			last_thres = pval_threshold

			#Filter pvalue file
			self.filter_percentile(1 - thres)

			#pre-process pvalues
			regions = self.get_regions(gene_radius=gene_radius)

			#Calculate observed candidate gene enrichment. 
			obs_enrichments = genr.calc_enrichment(all_genes, cand_gene_indices, regions)
			r1 = obs_enrichments[0] / float(obs_enrichments[1])
			r2 = (num_cand_genes / float(num_genes))
			obs_stat = sp.log(r1 / r2)
			print 'Observed statistics % f' % obs_stat
			log_enrichments.append(obs_stat)

			#What cand. genes overlap with regions?
			obs_cg_indices = obs_enrichments[2]
			if obs_cg_indices:
				obs_gene_str += str(pval_threshold) + ','
				tair_ids = [all_genes[cgi].tairID for cgi in obs_cg_indices]
				obs_gene_str += ','.join(tair_ids)
				obs_gene_str += '\n'



			for method in methods:
				if method == 'chi_square':
					chi_sq_pval, chi_sq_stat = genr.get_chi_square_pval(obs_enrichments[0],
											obs_enrichments[1],
											num_cand_genes, num_genes)
					method_res_dict[method]['statistics'].append(chi_sq_stat)
					method_res_dict[method]['pvals'].append(chi_sq_pval)

				if method == 'multinomial':
					pass
				if method == 'gene_perm':
					p_val, perm_stats = genr.get_gene_perm_pval(obs_stat, regions, all_genes,
									cand_gene_indices, num_perm=num_perm,
									early_stop_threshold=early_stop_threshold)
					method_res_dict[method]['statistics'].append(perm_stats)
					method_res_dict[method]['pvals'].append(p_val)

				if method == 'snps_perm':
					p_val, perm_stats = genr.get_snps_perm_pval(obs_stat, regions, all_genes,
									cand_gene_indices, chrom_ends,
									num_perm=num_perm,
									early_stop_threshold=early_stop_threshold)

					method_res_dict[method]['statistics'].append(perm_stats)
					method_res_dict[method]['pvals'].append(p_val)

#						h_res = pylab.hist(perm_stats)
#						pylab.vlines(obs_stat, 0, max(h_res[0]), colors='r')
#						pylab.savefig(env.env['tmp_dir'] + 'test.pdf', format='pdf')


		if obs_genes_file:
			with open(obs_genes_file, 'w') as f:
				f.write(obs_gene_str)

		#Now the plotting of the results.
		method_name_dict = {'chi_square':'Chi-square test', 'gene_perm':'Candidate gene permutations',
					'snps_perm':'SNP positions permutation (chromosome rotation)' }


		pval_thresholds.reverse()
		pos_list = range(len(pval_thresholds))
		for m in methods:
			pylab.figure()
			pvals = []
			for p in method_res_dict[m]['pvals']:
				if p != 0:
					pvals.append(p)
				else:
					pvals.append(0.5 / num_perm)
			neg_log_pvals = map(lambda x:-math.log10(x), pvals)
			pylab.barh(pos_list, neg_log_pvals, align='center', color='g', alpha=0.6)
			pylab.ylim((-1, len(pval_thresholds)))
			pylab.yticks(pos_list, map(str, pval_thresholds))
			pylab.xlabel('Enrichment -log(p-value)')
			pylab.ylabel('p-value percentile threshold')
			ymin, ymax = pylab.ylim()
			xmin, xmax = pylab.xlim()
			pylab.axvline(-math.log10(0.05), ymin=ymin, ymax=ymax, color='r')
			pylab.xlim((0, max(-math.log10(0.05), max(neg_log_pvals)) * 1.05))
			pylab.title(method_name_dict[m])
			pylab.savefig(file_prefix + '_' + m + '.png', format='png')


		pylab.figure()
		pylab.barh(pos_list, log_enrichments, align='center', color='g', alpha=0.6)
		pylab.ylim((-1, len(pval_thresholds)))
		pylab.yticks(pos_list, map(str, pval_thresholds))
		pylab.xlabel('Enrichment (log[ratio])')
		pylab.ylabel('p-value percentile threshold')
		pylab.savefig(file_prefix + '_enrichment_ratio.png', format='png')

		return {'enr_stats':log_enrichments, 'method_res_dict':method_res_dict}




	def _plot_small_manhattan_(self, pdf_file=None, png_file=None, min_score=0, max_score=None,
				type="pvals", ylab="$-$log$_{10}(p-$value$)$", plot_bonferroni=False,
				cand_genes=None, threshold=0, highlight_markers=None, chromosome=None,
				tair_file=None, plot_genes=True):
		import matplotlib
		#matplotlib.use('Agg')
		import matplotlib.pyplot as plt

		tair_genes = None
		min_x = min(self.snp_results['positions'])
		max_x = max(self.snp_results['positions'])
		if plot_genes:
			print 'Retrieving genes from DB.'
			gene_buffer = int((max_x - min_x) / 16) #bp
			tair_genes = get_gene_list(min_x, max_x, chromosome)
			print "Found", len(tair_genes), "genes in region:", min_x, "to", max_x
			g = tair_genes[0]
			gene_line_list = [g]
			gene_line_map = {g:1}
			used_lines = set([1])
			free_lines = set()#[2, 3, 4, 5, 6, 7, 8, 9, 10])
			for next_g in tair_genes[1:]:
				s_pos = next_g.startPos
				new_gene_line_list = [next_g]
				for g in gene_line_list:
					if next_g.startPos - gene_buffer <= g.endPos:
						new_gene_line_list.append(g)
					else:

						gl = gene_line_map[g]
						used_lines.remove(gl)
						free_lines.add(gl)
				#pdb.set_trace()
				if len(free_lines) > 0:
					gl = min(free_lines)
					free_lines.remove(gl)
				elif len(used_lines) > 0:
					gl = max(used_lines) + 1
				else:
					gl = 1
				used_lines.add(gl)
				gene_line_map[next_g] = gl
				#pdb.set_trace()
				gene_line_list = new_gene_line_list
			max_lines = max([gene_line_map[g] for g in tair_genes])
			min_score = -0.2 - 0.4 * max_lines


		if tair_file:
			with open(tair_file, 'w') as f:
				for gene in tair_genes:
					f.write(str(gene) + "\n")



		displayed_unit = 1000.0 #kbs
		scoreRange = max_score - min_score
		if plot_genes:
			plt.figure(figsize=(10, 4 + 0.5 * max_lines))
		else:
			plt.figure(figsize=(10, 3.5))
		plt.axes([0.045, 0.12, 0.95, 0.7])
		starPoints = [[], [], []]
		color = 'b'
		color_map = {1:'b', 2:'g', 3:'r', 4:'c', 5:'m'}
		if chromosome:
			color = color_map[chromosome]

		chrom = self.snp_results['chromosomes'][0]
		positions = map(lambda x: x / displayed_unit, self.snp_results['positions'])
		scores = self.snp_results['scores']
		for s_i, (score, pos) in enumerate(it.izip(scores, positions)):
			if score > max_score:
				starPoints[0].append(pos)
				starPoints[1].append(max_score)
				starPoints[2].append(score)
				score = max_score
			scores[s_i] = score

		plt.plot(positions, scores, ".", markersize=4, alpha=0.7, color=color)


		if tair_genes:
			print "Drawing TAIR genes"
			for g in tair_genes:
				y_value = -0.1 - 0.4 * gene_line_map[g]
				plt.text(g.startPos / displayed_unit, y_value + 0.08, g.tairID, size=5)
				if len(g.exons) > 0:

					for i in g.introns:
						plt.plot([i.startPos / displayed_unit, i.endPos / displayed_unit],
							[y_value, y_value], color=(0.6, 0.6, 0.6), linewidth=1)
					for e in g.exons:
						plt.plot([e.startPos / displayed_unit, e.endPos / displayed_unit],
							[y_value, y_value], color=(0.3, 0.3, 0.3), linewidth=3)
				else:
					plt.plot([g.startPos / displayed_unit, g.endPos / displayed_unit],
						[y_value, y_value], color=(0.3, 0.3, 0.3), linewidth=3)



		if cand_genes:
			for cg in cand_genes:
				plt.axvspan(cg.startPos / displayed_unit, cg.endPos / displayed_unit, facecolor='k',
						alpha=0.5, linewidth=0)


		if highlight_markers:
			ys = []
			xs = []
			for c, p, score in highlight_markers:
				xs.append(p / displayed_unit)
				if score > max_score:
					plt.text(x, max_score * 1.1, str(round(score, 2)), rotation=45, size="small")
					ys.append(max_score)
				else:
					ys.append(score)
			plt.plot(xs, ys, ".", color="#ff9944", markersize=6, alpha=0.7)

		if len(starPoints[0]) > 0:
			plt.plot(starPoints[0], starPoints[1], ".", color="#ee9922", markersize=4)
			i = 0
			while i < len(starPoints[0]):
				max_point = i
				cur_pos = starPoints[0][i]
				while i < len(starPoints[0]) and abs(starPoints[0][i] - cur_pos) < 1000000:
					if starPoints[2][i] > starPoints[2][max_point]:
						max_point = i
					i += 1
				plt.text(starPoints[0][max_point] - 200000, (starPoints[1][max_point] - 1) * 1.15, str(round(starPoints[2][max_point], 2)), rotation=45, size="small")




		if plot_bonferroni:
			b_threshold = -math.log10(1.0 / (len(scores) * 20.0))
			if threshold :
				plt.plot([0, max(positions)], [b_threshold, b_threshold], ":")
				threshold = -math.log10(threshold)
				plt.plot([0, max(positions)], [threshold, threshold], color='#6495ed', linestyle='-.')
			#Bonferroni threshold
			else:
				plt.plot([0, max(positions)], [b_threshold, b_threshold], color='#000000', linestyle="-.")

		x_range = max(positions) - min(positions)
		plt.axis([min(positions) - 0.05 * x_range, max(positions) + 0.05 * x_range, min_score - 0.05 * scoreRange, max_score + 0.05 * scoreRange])
		if not ylab:
			if type == "pvals":
				plt.ylabel('$ - log(p - $value$)$')

			else:
				plt.ylabel('score')
		else:
			plt.ylabel(ylab)
		plt.xlabel("kilobases")
		plt.title('Chromsome %d' % chrom)

		if pdf_file:
			plt.savefig(pdf_file, format="pdf")
		if png_file:
			plt.savefig(png_file, format="png", dpi=300, bbox_inches='tight')
		if not (pdf_file or png_file):
			plt.show()

		plt.clf()
		plt.close()


	def get_rare_haplotype_list(self, sd):
		"""
		Assumes SNPs are defined..
		"""
		#Get non-included accessions..
		#Find the same SNPs as in this object.
		#Sort by haplotypes.





	def plot_manhattan(self, pdf_file=None, png_file=None, min_score=None, max_score=None, percentile=90,
			type="pvals", ylab="$-$log$_{10}(p-$value$)$", plot_bonferroni=False, b_threshold=None,
			cand_genes=None, threshold=0, highlight_markers=None, tair_file=None, plot_genes=True,
			plot_xaxis=True):

		"""
		Plots a 'Manhattan' style GWAs plot.
		"""

		import matplotlib
		#matplotlib.use('Agg')
		import matplotlib.pyplot as plt
		num_scores = len(self.snp_results['scores'])

		"Plotting a Manhattan-style plot with %i markers." % num_scores

		chromosome_ends = self.get_chromosome_ends()
		result = self.simple_clone()

		chrom_set = set(result.snp_results['chromosomes'])
		chromosomes = list(chrom_set)
		chromosomes.sort()
		if len(chrom_set) == 1:
			percentile = 0.0
		if percentile != 0.0:
			result.filter_percentile(percentile / 100.0)

		if highlight_markers:
			new_h_markers = []
			for c, p, pval in highlight_markers:
				new_h_markers.append((c, p, -math.log10(pval)))
			highlight_markers = new_h_markers

		if not max_score:
			max_score = max(result.snp_results['scores'])
			if highlight_markers:
				h_scores = [s for c, p, s in highlight_markers]
				max_score = max(max_score, max(h_scores))
		if not min_score:
			if type == "pvals":
				min_score = 0
			else:
				min_score = min(result.snp_results['scores'])


		if cand_genes:
			#processing candidate genes by chromosome
			chr_cand_genes = {}
			for chrom in chromosomes:
				chr_cand_genes[chrom] = []
			for cg in cand_genes:
				chr_cand_genes[cg.chromosome].append(cg)

		if len(chrom_set) == 1:
			chrom = chrom_set.pop()
			if cand_genes:
				cand_genes = chr_cand_genes[chrom]
			return result._plot_small_manhattan_(pdf_file=pdf_file, png_file=png_file, min_score=min_score,
						max_score=max_score, ylab=ylab, plot_bonferroni=plot_bonferroni,
						cand_genes=cand_genes, threshold=threshold,
						highlight_markers=highlight_markers, chromosome=chrom,
						tair_file=tair_file, plot_genes=plot_genes)


		scoreRange = max_score - min_score
		offset = 0
		chromosome_splits = result.get_chromosome_splits()

		ticksList1 = []
		ticksList2 = []
		textPos = []
		plt.figure(figsize=(12, 2.8))
		plt.axes([0.045, 0.15, 0.95, 0.71])
		starPoints = [[], [], []]
		chr_offsets = []
		for i, chromosome_end in enumerate(chromosome_ends):
			chr_offsets.append(offset)
			index1 = chromosome_splits[i]
			index2 = chromosome_splits[i + 1]
			scoreList = result.snp_results['scores'][index1:index2]
			posList = result.snp_results['positions'][index1:index2]
			chrom = chromosomes[i]
			newPosList = [offset + pos for pos in posList]

			for s_i, (score, pos) in enumerate(it.izip(scoreList, newPosList)):
				if score > max_score:
					starPoints[0].append(pos)
					starPoints[1].append(max_score)
					starPoints[2].append(score)
					score = max_score
				scoreList[s_i] = score

			plt.plot(newPosList, scoreList, ".", markersize=2, alpha=0.7)

			if cand_genes:
				for cg in chr_cand_genes[chrom]:
					plt.axvspan(offset + cg.startPos, offset + cg.endPos, facecolor='#FF9900', alpha=0.6)

			oldOffset = offset
#			textPos.append(offset + chromosome_end / 2 - 2000000)
			offset += chromosome_end
			if plot_xaxis:
				for j in range(oldOffset, offset, 2000000):
					ticksList1.append(j)
				for j in range(0, chromosome_end, 2000000):
					if j % 4000000 == 0 and j < chromosome_end - 2000000 :
						ticksList2.append(j / 1000000)
					else:
						ticksList2.append("")



		plt.plot(starPoints[0], starPoints[1], ".", color="#ee9922", markersize=4)
		if len(starPoints[0]) > 0:
			i = 0
			while i < len(starPoints[0]):
				max_point = i
				cur_pos = starPoints[0][i]
				while i < len(starPoints[0]) and abs(starPoints[0][i] - cur_pos) < 3000000:
					if starPoints[2][i] > starPoints[2][max_point]:
						max_point = i
					i += 1
				plt.text(starPoints[0][max_point] - 1000000, (starPoints[1][max_point] - 1) * 1.15, str(round(starPoints[2][max_point], 2)), rotation=45, size="small")


		if highlight_markers:
			ys = []
			xs = []
			for c, p, score in highlight_markers:
				x = chr_offsets[c - 1] + p
				xs.append(x)
				if score > max_score:
					plt.text(x, max_score * 1.1, str(round(score, 2)), rotation=45, size="small")
					ys.append(max_score)
				else:
					ys.append(score)
			plt.plot(xs, ys, ".", color="#ff9944", markersize=6, alpha=0.8)


		if plot_bonferroni:
			if not b_threshold:
				b_threshold = -math.log10(1.0 / (num_scores * 20.0))
			if threshold :
				plt.plot([0, sum(result.chromosome_ends)], [b_threshold, b_threshold], ":")
				threshold = -math.log10(threshold)
				plt.plot([0, sum(result.chromosome_ends)], [threshold, threshold], color='#6495ed', linestyle='-.')
			#Bonferroni threshold
			else:
				plt.plot([0, sum(result.chromosome_ends)], [b_threshold, b_threshold], color='#000000', linestyle="-.")

		x_range = sum(result.chromosome_ends)
		plt.axis([-x_range * 0.01, x_range * 1.01, min_score - 0.05 * scoreRange, max_score + 0.05 * scoreRange])
		if plot_xaxis:
			plt.xticks(ticksList1, ticksList2, fontsize='x-small')
		#pdb.set_trace()
		if not ylab:
			if type == "pvals":
				plt.ylabel('$ - log(p - $value$)$')

			else:
				plt.ylabel('score')
		else:
			plt.ylabel(ylab)
		plt.xlabel("Mb")

		if pdf_file:
			plt.savefig(pdf_file, format="pdf")
		if png_file:
			plt.savefig(png_file, format="png", dpi=300, bbox_inches='tight')
		if not (pdf_file or png_file):
			plt.show()

		plt.clf()
		plt.close()



	def get_chromosome_splits(self):
		"""
		Returns list of indices (and prev chromosome), for the when the chromosomes
		change in the scores, positions indices.
		"""
		last_chrom = 0
		chromosome_splits = []
		for i, chrom in enumerate(self.snp_results['chromosomes']):
			if last_chrom != chrom:
				while last_chrom < chrom:
					last_chrom += 1
					chromosome_splits.append(i)
		chromosome_splits.append(i)
		return chromosome_splits


	def neg_log_trans(self):
		"""
		Apply - log(x) to the pvalues (scores)
		"""
		import math
		f = lambda x: 323.0 if x == 0.0 else - math.log10(x)
		self.snp_results['scores'] = map(f, self.snp_results['scores'])


	def filter_percentile(self, percentile, reversed=False):
		"""
		Filter above (or below) percentile.
		"""
		sl = self.snp_results['scores'][:]
		sl.sort()
		score_cutoff = sl[int(len(sl) * percentile)]
		self.filter_attr("scores", score_cutoff, reversed=reversed)


	def filter_common_top_scores(self, result, top_fraction=0.1):
		import bisect
		self._rank_scores_()
		result._rank_scores_()
		keep_indices_1 = set()
		keep_indices_2 = set()
		chr_pos_list_1 = self.get_chr_pos_list()
		chr_pos_list_2 = result.get_chr_pos_list()
		for i in self.orders[:int(top_fraction * len(self.orders))]:
			j = bisect.bisect(chr_pos_list_2, chr_pos_list_1[i]) - 1
			if chr_pos_list_2[j] == chr_pos_list_1[i]:
				keep_indices_1.add(i)

		for i in result.orders[:int(top_fraction * len(result.orders))]:
			j = bisect.bisect(chr_pos_list_1, chr_pos_list_2[i]) - 1
			if chr_pos_list_1[j] == chr_pos_list_2[i]:
				keep_indices_2.add(i)

		keep_indices_list_1 = list(keep_indices_1)
		for i in keep_indices_2:
			keep_indices_1.add(bisect.bisect(chr_pos_list_1, chr_pos_list_2[i]) - 1)

		for i in keep_indices_list_1:
			keep_indices_2.add(bisect.bisect(chr_pos_list_2, chr_pos_list_1[i]) - 1)

		print len(keep_indices_1), len(keep_indices_2)
		self.filter_indices(list(keep_indices_1))
		result.filter_indices(list(keep_indices_2))
		return keep_indices_1, keep_indices_2


	def filter_indices(self, indices_to_keep):
		indices_to_keep.sort()
		snp_results = {}
		for info in self.snp_results:
			snp_results[info] = []
		count = len(self.scores)
		for i in indices_to_keep:
			for info in self.snp_results:
				if self.snp_results[info]:
					snp_results[info].append(self.snp_results[info][i])

		self.snp_results = snp_results
		print "%i scores were removed." % (count - len(self.scores))


	def filter_attr(self, attr_name, attr_threshold, reversed=False, verbose=False, return_clone=False):
		"""
		Filter out scores / pvalues etc. which have attr < attr_threshold.

		attr are e.g.
		'mafs', 'macs', 'scores', etc.
		"""
		if verbose:
			print "Filtering for attribute ' % s' with threshold: %g" % (attr_name, attr_threshold)
		attr = self.snp_results[attr_name]
		snp_results = {}
		for info in self.snp_results:
			snp_results[info] = []
		count = len(self.snp_results['scores'])
		for i in range(count):
			if reversed:
				if attr[i] <= attr_threshold:
					for info in self.snp_results:
						if len(self.snp_results[info]) > 0:
							snp_results[info].append(self.snp_results[info][i])
			else:
				if attr[i] >= attr_threshold:
					for info in self.snp_results:
						if len(self.snp_results[info]) > 0:
							snp_results[info].append(self.snp_results[info][i])

		num_scores = len(snp_results['scores'])
		if verbose:
			print "%i scores were removed out of %i." % (count - num_scores, count)

		if return_clone:
			return Result(snp_results=snp_results, accessions=self.accessions, name=self.name,
		     			result_type=self.result_type, phen_id=self.phen_id, keys=self.keys)
		else:
			self.snp_results = snp_results
			return num_scores



	def get_min_distances(self, chrom_pos_list):
		"""
		Return distance statistics on this result object in relation to the chromosome position list.
		"""
		cpl = self.get_chr_pos_list()
		cp_array = sp.array(chrom_pos_list)
		cpl_dists = map(sp.absolute, [sp.array(cpt) - cp_array for cpt in cpl])
		min_dists = []
		for cp_dists in cpl_dists:
			curr_min = -1 #different chromosome.
			for cp_dist in cp_dists:
				if cp_dist[0] == 0:
					if curr_min == -1:
						curr_min = cp_dist[1]
					else:
						curr_min = min(curr_min, cp_dist[1])
			min_dists.append(curr_min)
		return min_dists


	def get_distances(self, chrom_pos_list):
		"""
		Return a list of lists of distances.
		"""
		cpl = self.get_chr_pos_list()
		cp_array = sp.array(chrom_pos_list)
		cpl_dists = map(sp.absolute, [sp.array(cpt) - cp_array for cpt in cpl])
		return cpl_dists


	def get_farthest_w_stronger_association(self, chrom_pos_list, assume_ranked=False):
		"""
		Return distance to farthest locus with greater association than some of the given loci.
		(chr_dist,pos_dist)
		"""
		if not assume_ranked:
			self._rank_scores_()
		c_indices = self.get_indices(chrom_pos_list)
		c_ranks = [self.ranks[i] for i in c_indices]
		max_rank_i = sp.argmax(c_ranks)
		stronger_indices = self.orders[:max_rank_i - 1]
		cpl = self.get_chr_pos_list()
		stronger_cpl = sp.array([cpl[i] for i in stronger_indices])
		cp = chrom_pos_list[0]
		dist_l = sp.absolute(stronger_cpl - sp.array(cp)).tolist()
		for cp in chrom_pos_list[1:]:
			dist_list = sp.absolute(stronger_cpl - sp.array(cp)).tolist()
			dist_l = [min(d1, d2) for d1, d2 in it.izip(dist_l, dist_list)]
		max_dist = max(dist_l)
		return max_dist










#	def filter_non_segregating_snps(self, ecotype1, ecotype2, accessions=None):
#		"""
#		Filter out all SNPs which are not segregating in the two accessions.
#
#		Assumes the accessions map the results objects SNPs. (and that they are defined)
#		"""
#		newScores = []
#		newPositions = []
#		newChromosomes = []
#		newMafs = []
#		newMarfs = []
#		new_snps = []
#
#		if accessions:
#			ecotypes = accessions
#		else:
#			ecotypes = self.accessions
#
#		e_i1 = ecotypes.index(ecotype1)
#		e_i2 = ecotypes.index(ecotype2)
#
#		for i in range(len(self.snps)):
#			snp = self.snps[i]
#			if snp[e_i1] != snp[e_i2]:
#				newScores.append(self.scores[i])
#				newPositions.append(self.positions[i])
#				newChromosomes.append(self.chromosomes[i])
#				newMafs.append(self.mafs[i])
#				newMarfs.append(self.marfs[i])
#				new_snps.append(snp)
#
#		self.scores = newScores
#		self.positions = newPositions
#		self.chromosomes = newChromosomes
#		self.mafs = newMafs
#		self.marfs = newMarfs
#		self.snps = new_snps


	def get_top_snps_result(self, n):
		"""
		returns top n SNPs
		"""
		import copy
		result = copy.deepcopy(self) #Cloning
		result.filter_top_snps(n)
		return result


	def get_top_snps(self, n):
		"""
		returns top n SNPs
		"""
		self._rank_scores_() #Making sure the ranks are updated
		chr_pos_list = []
		for i in self.orders[:n]:
			chr_pos_list.append((self.snp_results['chromosomes'][i], self.snp_results['positions'][i]))
		return chr_pos_list


	def get_indices(self, chr_pos_list):
		"""
		Return the indices of the given locus.
		"""
		import bisect
		indices = []
		for chromosome, position in chr_pos_list:
			cpl = self.get_chr_pos_list()
			i = bisect.bisect(cpl, (chromosome, position))
			if cpl[i - 1] == (chromosome, position):
				indices.append(i - 1)
			else:
				indices.append(-1)
		return indices



	def get_top_genes(self, n, window_size=5000, conn=None):
		"""
		Returns a set of (chromosome, start_pos, end_pos), for genes found.
		"""
		self._rank_scores_() #Making sure the ranks are updated
		genes = set()
		snp_ix = []
		i = 0
		while len(genes) < n:
			snp_i = self.orders[i]
			p = self.snp_results['positions'][snp_i]
			c = self.snp_results['chromosomes'][snp_i]
			c_genes = get_gene_list(p - window_size, p + window_size, c, include_intron_exons=False, conn=conn)
			for g in c_genes:
				genes.add((g.chromosome, g.startPos, g.endPos, g.tairID))
			#print len(genes)
			i += 1
		print 'found % d genes' % len(genes)
		return genes


	def get_chr_pos_score_list(self):
		return zip(self.snp_results['chromosomes'], self.snp_results['positions'], self.snp_results['scores'])


	def get_chrom_score_pos_dict(self):
		cps_list = self.get_chr_pos_score_list()
		d = {}
		for chrom in [1, 2, 3, 4, 5]:
			d[chrom] = {'scores':[], 'positions':[]}
		for chrom, pos, score in cps_list:
			d[chrom]['scores'].append(score)
			d[chrom]['positions'].append(pos)
		return d



	def get_chr_pos_list(self):
		return zip(self.snp_results['chromosomes'], self.snp_results['positions'])


	def get_ranks_of_loci(self, chr_pos_list):
		"""
		Returns the ranks of loci.
		"""


	def analyse_loci(self,):
		"""
		Performs basic analysis on the given loci, such as 
			- ranks, 
			- distance from most significant
			- distances from all significant
			, of all SNPs with gre
		"""

	def get_top_regions(self, n, distance_threshold=25000):
		"""
		Returns a list of regions, defined by (chromosome, start_pos, end_pos).
		"""
		self._rank_scores_()
		chromosome_ends = self.get_chromosome_ends()

		i = 0 #how many SNPs are needed
		region_count = 0
		regions = [[], [], [], [], []]  #a list of (chromosome,start,stop)
		while sum(map(len, regions)) < n:
			snp_i = self.orders[i]
			snp_pos = self.snp_results['positions'][snp_i]
			chromosome = self.snp_results['chromosomes'][snp_i]
			new_snp_reg = (chromosome, max(snp_pos - distance_threshold, 0), min(snp_pos + distance_threshold, chromosome_ends[chromosome - 1]))
			covered = False
			for reg_i, reg in enumerate(regions[chromosome - 1]):
				if (new_snp_reg[1] < reg[2] and new_snp_reg[2] > reg[2]) or (new_snp_reg[1] < reg[1] and new_snp_reg[2] > reg[1]): #then merge					
					regions[chromosome - 1].pop(reg_i) #Removing the region which we're merging with.
					new_snp_reg = (reg[0], min(reg[1], new_snp_reg[1]), max(reg[2], new_snp_reg[2]))
				elif (new_snp_reg[1] >= reg[1] and new_snp_reg[2] <= reg[2]):
					covered = True  #It's already covered
					break
			if not covered:
				regions[chromosome - 1].append(new_snp_reg)
			i += 1
		print "It took %i SNPS to find %i regions." % (i, n)
		regions_list = regions[0] + regions[1] + regions[2] + regions[3] + regions[4]
		return regions_list


	def get_regions(self, gene_radius=20000):
		"""
		Returns a list of [[chr,start],[chr,end]] elements, 
		"""
		regions = []
		num_scores = len(self.snp_results['scores'])
		iter = enumerate(it.izip(self.snp_results['chromosomes'], self.snp_results['positions']))
		th = sp.array([0, gene_radius])
		(pos_i, t) = iter.next()
		chr_pos = sp.array(t)
		while pos_i < num_scores - 1:
			start_chr_pos = chr_pos
			end_chr_pos = start_chr_pos
			while pos_i < num_scores - 1 and sp.all((chr_pos - end_chr_pos) <= th):
				end_chr_pos = chr_pos
				(pos_i, t) = iter.next()
				chr_pos = sp.array(t)
			if tuple(end_chr_pos) < tuple(start_chr_pos):
				pdb.set_trace()
			if pos_i == num_scores - 1: #Last one..
				if sp.all((chr_pos - end_chr_pos) <= th):
					regions.append([start_chr_pos - th, chr_pos + th])
				else:
					regions.append([start_chr_pos - th, end_chr_pos + th])
					regions.append([chr_pos - th, chr_pos + th])

			else:
				regions.append([start_chr_pos - th, end_chr_pos + th])
		start_chr_pos = chr_pos
		end_chr_pos = start_chr_pos
		regions.append([start_chr_pos - th, end_chr_pos + th])
		print 'Found % d regions' % len(regions)
		return regions


	def count_nearby(self, chrom_pos_list, radius):
		"""
		Counts how many of the chrom_pos tuples are within a radius of a SNP in the result.
		"""
		cpl = self.get_chr_pos_list()
		count = 0
		for (chrom1, pos1) in chrom_pos_list:
			for (chrom2, pos2) in cpl:
				if chrom1 == chrom2 and abs(pos1 - pos2) <= radius:
					count += 1
					break
		return count


	def get_power_analysis(self, caus_chrom_pos_list, window_sizes=[0]):
		"""
		Calculate Power and FDR..
		"""
		for window_size in __window_sizes:
			cpl = self.get_chr_pos_list()
			num_caus_found = 0
			num_false_found = 0
			for (chrom1, pos1) in caus_chrom_pos_list:
				caus_found = False
				for (chrom2, pos2) in cpl:
					if chrom1 == chrom2 and abs(pos1 - pos2) <= radius:
						caus_found = True
					else:
						num_false_found += 1
				if caus_found:
					num_caus_found += 1
			tprs.append(float(num_caus_found) / len(caus_chrom_pos_list))
			fdrs.append(float(num_false_found) / len(cpl))

		return tprs, fdrs



	def get_region_result(self, chromosome, start_pos, end_pos, buffer=0):
		"""
		returns a result object with only the SNPs, etc. within the given boundary.
		"""
		snp_results = {}
		for k in self.keys:
			snp_results[k] = []
		i = 0
		start_pos -= buffer
		end_pos += buffer

		chromosomes = self.snp_results['chromosomes']
		positions = self.snp_results['positions']
		while i < len(chromosomes) and chromosomes[i] != chromosome:
			i += 1
		while i < len(chromosomes) and positions[i] < start_pos:
			i += 1
		if i == len(chromosomes):
			raise Exception("region results never found!")

		while i < len(chromosomes) and positions[i] <= end_pos and chromosomes[i] == chromosome:
			for k in self.keys:
				snp_results[k].append(self.snp_results[k][i])
			i += 1

		return Result(result_type=self.result_type, snp_results=snp_results, phen_id=self.phen_id,
			accessions=self.accessions, keys=self.keys[:])



	def get_top_region_results(self, n, distance_threshold=25000, buffer=25000):
		reg_results = []
		i = 0
		for reg in self.get_top_regions(n, distance_threshold):
			reg_results.append(self.get_region_result(*reg, buffer=buffer))
		print "Regions were retrieved."
		return reg_results




	def get_chromosome_ends(self):
		if not self.chromosome_ends:
			self._sort_by_chr_pos_()
			i = 0
			chromosome_ends = []
			chromosomes = self.snp_results['chromosomes']
			while i < len(chromosomes):
				curr_chr = chromosomes[i]
				while i < len(chromosomes) and chromosomes[i] == curr_chr:
					i += 1
				chromosome_ends.append(self.snp_results['positions'][i - 1])
			self.chromosome_ends = chromosome_ends
		return self.chromosome_ends



	def clone(self):
		import copy
		result = copy.deepcopy(self) #Cloning
		return result


	def simple_clone(self):
		snp_results = {}
		for k in self.keys:
			snp_results[k] = self.snp_results[k][:]
		if self.accessions:
			accessions = self.accessions[:]
		else:
			accessions = None
		result = Result(snp_results=snp_results, accessions=accessions, keys=self.keys[:])
		result.chromosome_ends = self.chromosome_ends[:]
		return result


#
#
#	def na_mafs(self, min_maf=10):
#		"""
#		NA scores/pvalues which have maf<minMaf.		
#		"""
#		for i in range(0, len(self.scores)):
#			if self.mafs[i] < min_maf:
#				self.scores[i] = "NA"


	def get_num_scores(self):
		return len(self.snp_results['scores'])


	def get_max_snp(self):
		max_val = max(self.snp_results['scores'])
		mi = self.snp_results['scores'].index(max_val)
		return (self.snp_results['snps'][mi], self.snp_results['scores'][mi],
			self.snp_results['chromosomes'][mi], self.snp_results['positions'][mi])


	def min_score(self):
		return min(self.snp_results['scores'])

	def max_score(self):
		return max(self.snp_results['scores'])


	def arg_min_attr(self, attr_name='scores'):
		"""
		Returns the index of the minimum attribute.
		"""
		return sp.argmin(self.snp_results[attr_name])

	def arg_max_attr(self, attr_name='scores'):
		"""
		Returns the index of the maximum attribute.
		"""
		return sp.argmax(self.snp_results[attr_name])


	def num_scores(self):
		return len(self.snp_results['scores'])

	def write_to_file(self, filename, additional_columns=None, auto_pickling_on=True, only_pickled=False):
		columns = ['chromosomes', 'positions', 'scores', 'mafs', 'macs']
		if additional_columns:
			for info in additional_columns:
				if info in self.snp_results:
					columns.append(info)
#		try:
		if not only_pickled:
			with open(filename, "w") as f:
				f.write(','.join(columns) + "\n")
				for i in range(len(self.snp_results[columns[0]])):
					l = [self.snp_results[c][i] for c in columns]
					l = map(str, l)
					f.write(",".join(l) + "\n")
		if auto_pickling_on or only_pickled:
			pickle_file = filename + '.pickled'
			d = {}
			for attr in self.pickle_attributes:
				if attr == 'snp_results':
					snp_results = getattr(self, attr)
					d = {}
					for k in snp_results:
						if k != 'snps':
							d[k] = snp_results[k]
					d[attr] = d
				else:
					d[attr] = getattr(self, attr)
			with open(pickle_file, 'wb') as f:
				cPickle.dump(d, f)

#		except Exception, err_str:
#			print 'Failed writing the resultfile:', err_str
#			print 'Make sure the given path is correct, and you have write rights.'


class SNP(object):
	"""
	A class to represent a SNP. 

	It's used only when analysing a SNP.
	"""

	def __init__(self, position, chromosome, accessions=None, alleles=None, snpsds=None, score=None, rank=None, regionRank=None):
		self.position = position
		self.chromosome = chromosome

		self.alleles = None
		self.accessions = None
		self.score = None
		self.rank = None
		self.regionRank = None

		if not alleles and snpsds:
			self._getAllele_(snpsds)
		else:
			self.alleles = alleles
		if accessions:
			self.accessions = accessions
		if score:
			self.score = score
		if rank:
			self.rank = rank
		if regionRank:
			self.regionRank = regionRank

	def _getAllele_(self, snpsds):
		chr = self.chromosome - 1
		snpsd = snpsds[chr]
		self.snpsdIndex = -1
		self.alleles = None
		for i in range(0, len(snpsd.positions)):
			if snpsd.position == snpsd.positions[i]:
				self.alleles = snpsd.snps[i]
				self.snpsdIndex = i
				break
		if not self.alleles:
			print "The corresponding allele was not found in the data."




class Gene(object):
	"""
	A class which encompasses basic information about a gene.
	"""
	def __init__(self, chromosome=None, startPos=None, endPos=None, name="", description=None, dbRef="", tairID=""):
		self.chromosome = chromosome
		self.startPos = startPos
		self.endPos = endPos
		self.exons = []
		self.introns = []
		self.tairID = tairID
		self.dbRef = dbRef
		self.name = name
		self.description = description
		self.functionDescriptions = []
		self.shortDescriptions = []
		self.direction = None
		self.highlight = False


	def __str__(self):
		if not self.description:
			return "Chromosome=" + str(self.chromosome) + ", position=(" + str(self.startPos) + "," + str(self.endPos) + "), tair ID=" + self.tairID + ", short descriptions=" + str(self.shortDescriptions) + ", function descriptions=" + str(self.functionDescriptions) + "."
		else:
			return "Chromosome=" + str(self.chromosome) + ", position=(" + str(self.startPos) + "," + str(self.endPos) + "), tdbRef=" + self.dbRef + ", name=" + str(self.name) + ", description=" + str(self.description) + "."


	def _update_introns_(self):
		"""
		updates the introns, given the exons.  It uses the Region object..
		"""
		introns = []
		for i in range(len(self.exons) - 1):
			e1 = self.exons[i]
			e2 = self.exons[i + 1]
			intron = Region(self.chromosome, e1.endPos, e2.startPos)
			introns.append(intron)
		self.introns = introns


def getCandidateGeneList(cgl_id, host="papaya.usc.edu", user="bvilhjal", passwd="bjazz32", db="stock_250k", tair9=True):
	import MySQLdb
	#Load cand. gene list.	
	print "Connecting to db, host=" + host
	if not user:
		import sys
		sys.stdout.write("Username: ")
		user = sys.stdin.readline().rstrip()
	if not passwd:
		import getpass
		passwd = getpass.getpass()
	try:
		conn = MySQLdb.connect (host=host, user=user, passwd=passwd, db=db)
	except MySQLdb.Error, e:
		print "Error %d: %s" % (e.args[0], e.args[1])
		sys.exit (1)
	cursor = conn.cursor ()
	#Retrieve the filenames
	print "Fetching data"

	#select c.locustag, b.start, b.stop, a.comment from genome.gene_commentary a, genome.entrezgene_mapping b, genome.gene c where b.start > 25000 and b.stop < 75000 and b.chromosome=1 and b.gene_id = c.gene_id and c.gene_id = a.gene_id and a.gene_commentary_type_id = 8
	#select distinct t8_fd.tair_id, t8.chromosome, t8.start, t8.end, t8_fd.type, t8_fd.short_description from T8_annotation_TH.t8_063008 t8, T8_annotation_TH.t8_func_desc t8_fd, stock_250k.candidate_gene_list cgl where t8.pub_locus+'.1' = t8_fd.tair_id and cgl.list_type_id=129  and cgl.original_name=t8.pub_locus and t8.chromosome =1 order by t8.chromosome, t8.start
	#select distinct gm.chromosome, gm.start, gm.stop, g.locustag from genome.entrezgene_mapping gm, genome.gene g, stock_250k.candidate_gene_list cgl where cgl.list_type_id=129 and gm.gene_id = g.gene_id and cgl.gene_id=g.gene_id order by gm.chromosome, gm.start, gm.stop

	numRows = int(cursor.execute("select distinct gm.chromosome, gm.start, gm.stop, g.locustag, g.gene_symbol, g.description, g.dbxrefs from genome.entrezgene_mapping gm, genome.gene g, stock_250k.candidate_gene_list cgl where cgl.list_type_id=" + str(cgl_id) + " and gm.gene_id = g.gene_id and cgl.gene_id=g.gene_id order by gm.chromosome, gm.start, gm.stop"))
	if tair9:
		t_map = tc.tair8_to_tair9_map()
	candGenes = []
	while(1):
		row = cursor.fetchone()
		if not row:
			break;
		if tair9:
			start_pos = t_map.get_tair9_pos(int(row[1]))
			end_pos = t_map.get_tair9_pos(int(row[2]))
		else:
			start_pos = int(row[1])
			end_pos = int(row[2])
		gene = Gene(int(row[0]), start_pos, end_pos, name=row[4], description=row[5], dbRef=row[6])
		candGenes.append(gene)
	cursor.close ()
	conn.close ()
	print "Candiate gene-lists fetched"
	return candGenes


def get_gene_list(start_pos=None, end_pos=None, chr=None, include_intron_exons=True, \
		verbose=True, conn=None, tair9=True):
	"""
	Fetch genes within a region or all genes from DB.
	"""
	import dbutils
	if not conn:
		new_conn = dbutils.connect_to_default_lookup('genome')
		cursor = new_conn.cursor()
	else:
		cursor = conn.cursor()
	#Retrieve the filenames
	#if verbose:
	#	print "Fetching data"  
	#print "Fetching data"  

	#select c.locustag, b.start, b.stop, a.comment from genome.gene_commentary a, genome.entrezgene_mapping b, genome.gene c where b.start > 25000 and b.stop < 75000 and b.chromosome=1 and b.gene_id = c.gene_id and c.gene_id = a.gene_id and a.gene_commentary_type_id = 8
	#select distinct t8_fd.tair_id, t8.chromosome, t8.start, t8.end, t8_fd.type, t8_fd.short_description from T8_annotation_TH.t8_063008 t8, T8_annotation_TH.t8_func_desc t8_fd, stock_250k.candidate_gene_list cgl where t8.pub_locus+'.1' = t8_fd.tair_id and cgl.list_type_id=129  and cgl.original_name=t8.pub_locus and t8.chromosome =1 order by t8.chromosome, t8.start
	#select distinct gm.chromosome, gm.start, gm.stop, g.locustag from genome.entrezgene_mapping gm, genome.gene g, stock_250k.candidate_gene_list cgl where cgl.list_type_id=129 and gm.gene_id = g.gene_id and cgl.gene_id=g.gene_id order by gm.chromosome, gm.start, gm.stop
	if chr and start_pos and end_pos:
		sql_statement = "SELECT DISTINCT gm.chromosome, gm.start, gm.stop, g.locustag, \
		g.gene_symbol, g.description, g.dbxrefs FROM genome.entrezgene_mapping gm, genome.gene g WHERE \
		gm.gene_id = g.gene_id AND gm.chromosome=" + str(chr) + " AND gm.stop>" + str(start_pos) + " AND \
		gm.start<" + str(end_pos) + " ORDER BY gm.chromosome, gm.start, gm.stop"
		numRows = int(cursor.execute(sql_statement))
	else:
		sql_statement = "select distinct gm.chromosome, gm.start, gm.stop, g.locustag, g.gene_symbol, \
			g.description, g.dbxrefs from genome.entrezgene_mapping gm, genome.gene g \
			where gm.gene_id = g.gene_id order by gm.chromosome, gm.start, gm.stop"
		numRows = int(cursor.execute(sql_statement))
	if numRows == 0:
		pass
		#print sql_statment

	if tair9:
		t_map = tc.tair8_to_tair9_map()
	genes = []
	while(1):
		row = cursor.fetchone()
		if not row:
			break;
		#try:
			#chr, start, stop, gene_symbol, description, dbref,  
		if row[1] and  row[2] and row[0] in ['1', '2', '3', '4', '5']:
			chrom = int(row[0])
			if tair9:
				start_pos = t_map.get_tair9_pos(chrom, int(row[1]))
				end_pos = t_map.get_tair9_pos(chrom, int(row[2]))
			else:
				start_pos = int(row[1])
				end_pos = int(row[2])

			gene = Gene(int(row[0]), start_pos, end_pos, name=row[4], description=row[5], dbRef=row[6], tairID=row[3])
			gene.tairID = row[6][5:]
			genes.append(gene)
#		except Exception, err_str:
#			#pass
#			if verbose:
#				print err_str, ':'
#				print row

	if include_intron_exons:
		for g in genes:
			sql_stat = "select distinct gs.start, gs.stop, g.dbxrefs from genome.entrezgene_mapping gm, genome.gene g, genome.gene_segment gs, genome.gene_commentary gc where g.dbxrefs='" + str(g.dbRef) + "' and gm.gene_id = g.gene_id and gs.gene_commentary_id=gc.id and gc.gene_id=gm.gene_id order by gs.start, gs.stop"
			numRows = int(cursor.execute(sql_stat))
			segments = []
			while(1):
				row = cursor.fetchone()
				if not row:
					break;
				try:
					if tair9:
						start_pos = t_map.get_tair9_pos(g.chromosome, int(row[0]))
						end_pos = t_map.get_tair9_pos(g.chromosome, int(row[1]))
					else:
						start_pos = int(row[0])
						end_pos = int(row[1])
					segments.append(Region(g.chromosome, start_pos, end_pos))
				except Exception, err_str:
					print err_str, ':'
					print row
			exons = []
			i = 1
			#print "len(segments):",len(segments)
			while i < len(segments):
				curr_exon = segments[i - 1]
				while i < len(segments) and curr_exon.overlapping(segments[i]):
					curr_exon.merge(segments[i])
					i += 1
				exons.append(curr_exon)
				i += 1
			g.exons = exons
			g._update_introns_()
			#print "len(exons):",len(exons)
			#for e in g.exons:
			#	print e.startPos, e.endPos
			#print "len(g.introns):",len(g.introns)
			#for e in g.introns:
			#	print e.startPos, e.endPos
	cursor.close()
	if not conn:
		new_conn.close ()
	if verbose:
		print "Gene-lists fetched"
	return genes


#"""
#select distinct gm.chromosome, gm.start, gm.stop, g.locustag, g.gene_symbol, g.description, g.dbxrefs, gs.start, gs.stop
#from genome.entrezgene_mapping gm, genome.gene g, gene_segment gs, gene_commentary gc
#where gm.gene_id = g.gene_id and gm.chromosome=5 and gm.stop>3173000 and gm.start<3181000 and gs.gene_commentary_id=gc.id and gc.gene_id=gm.gene_id 
#order by gm.chromosome, gm.start, gm.stop, gs.start, gs.stop
#"""


def load_cand_genes_file(file_name, format=1):
	"""
	Loads a candidate gene list from a csv file...
	"""
	f = open(file_name, "r")
	#print(f.readline())
	reader = csv.reader(f)
	tair_ids = []
	gene_names = []
	if format == 1:
		for row in reader:
		      tair_ids.append(row[0].upper())
		      if len(row) > 1:
		      	gene_names.append(row[1])
	f.close()
	return get_genes_w_tair_id(tair_ids), tair_ids




def get_genes_w_tair_id(tair_ids, tair9=True):
	conn = dbutils.connect_to_default_lookup("genome")
	cursor = conn.cursor()
	if tair9:
		t_map = tc.tair8_to_tair9_map()
	genes = []
	#print tair_ids
	tair_ids.sort()
	for tair_id in tair_ids:
		sql_statment = "select distinct gm.chromosome, gm.start, gm.stop, g.locustag, g.gene_symbol, g.description, g.dbxrefs from genome.entrezgene_mapping gm, genome.gene g where g.dbxrefs='TAIR:" + tair_id.upper() + "' and gm.gene_id = g.gene_id order by gm.chromosome, gm.start, gm.stop"
		#print sql_statment
		numRows = int(cursor.execute(sql_statment))
		if numRows > 1:
			print "Found 2 copies:", sql_statment
		while(1):
			row = cursor.fetchone()
			if not row:
				break;
			try:
				if row[1] and  row[2]:
					chrom = int(row[0])
					if tair9:
						start_pos = t_map.get_tair9_pos(chrom, int(row[1]))
						end_pos = t_map.get_tair9_pos(chrom, int(row[2]))
					else:
						start_pos = int(row[1])
						end_pos = int(row[2])
					#chr, start, stop, gene_symbol, description, dbref,  
					gene = Gene(chrom, start_pos, end_pos, name=row[4], description=row[5], dbRef=row[6], tairID=row[3])
					gene.tairID = row[6][5:]
					genes.append(gene)
			except Exception, err_str:
				pass
				print err_str, ':'
				print row

	cursor.close()
	conn.close()
	return genes


def get_result_filename(cm_id, pm_id, am_id):
	"""
	Retrieve the filename with the results.
	
	call method ID
	pehnotype method ID
	analysis method ID
	"""
	conn = dbutils.connect_to_default_lookup('stock_250k')
	cursor = conn.cursor ()

	#Retrieve the filenames
	print "Fetching data"
	numRows = int(cursor.execute("select rm.filename from stock_250k.results_method rm \
				where rm.call_method_id=%d and rm.phenotype_method_id=%d and analysis_method_id=%d " \
				% (cm_id, pm_id, am_id)))
	filenames = []
	while(1):
		row = cursor.fetchone()
		if not row:
			break;
		filenames.append(row[0])
	cursor.close ()
	conn.close ()
	return filenames




def load_result_from_db(pid, aid, cmid=54, host='gmi-ara-devel-be', conn=None):
	"""
	Imports a result object from the filesystem/DB.
	"""
	import dbutils
	if conn:
		cursor = conn.cursor()
	else:
		new_conn = dbutils.connect_to_db(host, 'stock_250k')
		cursor = new_conn.cursor()
	sql_statement = "SELECT short_name, filename, id FROM stock_250k.results_method \
			 WHERE call_method_id=%d and phenotype_method_id=%d and analysis_method_id=%d"\
			 % (cmid, pid, aid)
	#print sql_statement
	numRows = int(cursor.execute(sql_statement))
	row = cursor.fetchone()
	r, r_id = None, None
	if row:
		fname = row[1]
		r_id = int(row[2])
		print "File for %s, found at:%s" % (row[0], fname)
		r = Result(fname)

	else:
		print "Result not found with pid=%d, aid=%d, cmid=%d" % (pid, aid, cmid)

	cursor.close ()
	if not conn:
		new_conn.close ()

	return r, r_id


#def loadResults(phenotypeIndices, resultTypes=None, phed=None, snpsds=None, filterPercentile=None, filterCutoffs=None, phenotypeFile="/Network/Data/250k/dataFreeze_080608/phenotypes_all_raw_111008.tsv", secondRun=False):
#
#	if not phed:
#		phed = phenotypeData.readPhenotypeFile(phenotypeFile, delimiter='\t')
#
#	if not resultTypes:
#		if secondRun:
#			resultTypes = _getStandardSecondRunResultTypes_()
#		else:
#			resultTypes = _getStandardResultTypes4_()
#
#	results_map = {}
#	for i in phenotypeIndices:
#
#		results = []
#		for j in range(0, len(resultTypes)):
#			resultType = resultTypes[j]
#			phenName = phed.getPhenotypeName(i)
#			if phenName:
#				resultFile = resultType.getFileName(phed, i, secondRun=(secondRun and j % 2 == 1))  #Modify back to get old results 120708
#				try:
#					print "Loading result file", resultFile
#					if snpsds:
#						result = SNPResult(resultFile, snpsds=snpsds, name=str(resultType) + "_" + phenName, resultType=resultType, phenotypeID=i)
#					else:
#						result = Result(resultFile, name=str(resultType) + "_" + phenName, resultType=resultType, phenotypeID=i)
#					if resultType.logTransform:
#						print "Log transformed the p-values"
#						result.negLogTransform()
#
#					result.filterMARF(minMaf=resultType.mafCutoff)
#					#result.filterMAF(minMaf=resultType.mafCutoff)
#					if filterPercentile:
#						result.filterPercentile(filterPercentile)
#					elif filterCutoffs:
#						result.filterScoreCutoff(filterCutoffs[j])
#
#					results.append(result)
#				except Exception, e:
#					print e.message
#					print "Couldn't load", resultFile
#
#		results_map[i] = results
#		gc.collect()  #Calling garbage collector, in an attempt to clean up memory..
#	return results_map


#def _add_results_to_db_():
#	"""
#	TEST
#	"""
#	host = 'gmi-ara-devel-be'
#	for pid in range(1, 2):
#		for aid in range(20):
#			for cmid in range(60):
#				try:
#					r, r_id = load_result_from_db(pid, aid, cmid, host=host)
#					if not r:
#						raise Exception
#					r.insert_into_db(r_id)
#		                except Exception, err_str:
#		                	print "File not found, pid=%d, aid=%d, cmid=%d, err_str=%s" % (pid, aid, cmid, err_str)





if __name__ == "__main__":
	#load_cand_genes_file("/Users/bjarnivilhjalmsson/Projects/Ales_Pecinka/UV_cand_genes_021710.csv")
	#load_result_from_db(513,1)
	_add_results_to_db_()
