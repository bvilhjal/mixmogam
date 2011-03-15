"""
Usage: gwa.py [OPTIONS] 

Option:

	-h					Show this help
	-d ...					Filter SNPs randomly (speed-up for debugging)

	-o ...					ID string, used for output files.
	-i ...					The phenotype IDs, to be run. 

	-t ...					What data set is used.  Default is 75.
	-f ...					Load a specific data file, e.g. for heteroplasmy.
	-r ...					Phenotype file, if left out then phenotypes are retireved from the DB 
						(transformed values).
	-k ...					Specify the file containing the kinship matrix.  
						(otherwise default file is used or it's generated.)

	-a ...					Apply specific methods, otherwise all available are applied:
						lm, emma, emmax, kw, ft, emmax_anova, lm_anova, emmax_step etc.
	-b ...				 	Apply a transformation to the data, default is none, other possibilities are 
						log, sqrt, exp, sqr, most_normal (picks a most Gaussian looking transformation).
	-c ...					Should phenotype outliers be removed.  0 (no fence) is the default, 
						else the outlier fence is given in IQRs. (An int is required).
												 
	-u					Use existing results when possible, to speed up.  
						(Otherwise existing files are overwritten)

	-e					Generate analysis plots, histograms, QQ-plots, test different 
						transformations, etc.
	-n					Adds the result(s) to the DB.		
	--comment=...				Comment for DB. (Only applicable if result is added to DB.)	
	--no_phenotype_ids			Phenotypes don't have DB id as an prefix in their names.  
						(The phenotype index should be used instead.)
	--suggest_new_genotypes=...		Suggest new genotypes base on the top X snps, where X is the given argument

	--region_plots=...			Include region plots for the top N (given num) peaks.
	--cand_genes_file=...			A file with a list of candidate genes.  (added to plots)	

	-m ...					MAC threshold which is used for LM, EMMA, EMMAX, etc.  Default is 15.
	--data_format=...			What type of data should it aim for, binary (default), int, float, etc.
	--emmax_perm=...			Number of permutations, used for estimating a significant threshold.
	--with_replicates			Run EMMAX with replicates (if any, otherwise it uses the mean)
	--with_betas				Output betas (effect sizes), this is a tad slower
	--num_steps=...				Max number of steps, for EMMAX stepwise

	--local_gwas=chrom,start,stop		Run local GWAs only..		
	
	--use_imputed_full_data			Use the imputed full data, 1.7 million SNPs (memory heavy).
	
	#ONLY APPLICABLE FOR CLUSTER RUNS
	-p ...					Run mapping methods on the cluster with standard parameters.  The argument is used for runid 
						as well as output files.  Note that if this option is used then no output file should be specified.
						If phenotype index is missing, then all phenotypes are used.
	-q ...					Request memory (on cluster), otherwise use defaults 4GB for Kruskal-Wallis, 12GB for Emma.
	-l ...			 		Request time limit (on cluster), otherwise use defaults
	--proc_per_node=...			Request number of processors per node, default is 8. (Works only for EMMA.)
	--only_add_2_db				Does not submit jobs, but only adds available result files to DB. (hack for usc hpc)
	
	
	
						
	#NOT YET IMPLEMENTED
	--cofactor_chr_pos=...			A list of SNP (chromosome,positions) to be added as co-factors in the analysis.
	--cofactor_phen_id=...			A list of SNP positions to be added as co-factors in the analysis.
	--cofactor_file=...			A file specifying the cofactor.
	--cofactor_no_interact			Exclude interaction terms for cofactors.

Examples:
~/gwas_data$ python gwa.py -o test -i 1,5 -a kw,emmax -r ~/Projects/Data/phenotypes/phen_raw_092910.tsv 
Description:
  Applies various GWA methods to to phenotypes.

  Methods include: Kruskal-Wallis, Fisher's Exact, EMMA, EMMAX, etc.
  
  If PHENOTYPE_DATA_FILE is left out, then papaya DB is searched.  
  A PHENOTYPE_ID is required.
    
"""
import sys, getopt, traceback
import os

import phenotypeData
import dataParsers
import snpsdata
import gwaResults
import util
import warnings
import time
import cPickle

import scipy as sp
import linear_models as lm
from numpy import *
from env import *
import copy
import pdb

#For annoying linux computers, which don't have a display..
import matplotlib
matplotlib.use("Agg")


transformation_method_dict = {
			'sqrt':0,
			'none':1,
			'log':2,
			'sqr':3,
			'exp':4,
			}


analysis_methods_dict = {"kw":1,
			 "ft":2,
			 "emma":4,
			 'lm':16,
			 "emmax":32,
			 'emmax_anova':None,
			 'emmax_step':None,
			 'lm_anova':None,
			 }


def _parse_pids_(pid_arg_str):
	t_pids = pid_arg_str.split(',')
	pids = []
	for s in t_pids:
		if '-' in s:
			pid_range = map(int, s.split('-'))
		        for pid in range(pid_range[0], pid_range[1] + 1):
		        	pids.append(pid)
		else:
			pids.append(int(s))
	return pids


def parse_parameters():
	'Parse the parameters into a dict, etc.'
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)

	long_options_list = ["comment=", 'no_phenotype_ids', 'region_plots=', 'cand_genes_file=', 'proc_per_node=',
			'only_add_2_db', 'data_format=', 'emmax_perm=', 'with_replicates', 'with_betas', 'num_steps=',
			'local_gwas=', 'use_imputed_full_data']
	try:
		opts, args = getopt.getopt(sys.argv[1:], "o:i:p:a:b:c:d:ef:t:r:k:nm:q:l:hu", long_options_list)

	except:
		traceback.print_exc()
		print __doc__
		sys.exit(2)


	p_dict = {'run_id':'donald_duck', 'parallel':None, 'add_to_db':False, 'comment':'', 'mem_req':'1800mb',
		'call_method_id':75, 'walltime_req':'12:00:00', 'proc_per_node':8,
		'specific_methods':['kw', 'ft', 'lm', 'emma', 'emmax'], 'specific_transformations':['none'],
		'remove_outliers':0, 'kinship_file':None, 'analysis_plots':False, 'use_existing_results':False,
		'region_plots':0, 'cand_genes_file':None, 'debug_filter':1, 'phen_file':None,
		'no_phenotype_ids':False, 'only_add_2_db':False, 'mac_threshold':15, 'data_file':None,
		'data_format':'binary', 'emmax_perm':None, 'with_replicates':False, 'with_betas':False,
		'num_steps':50, 'local_gwas':None, 'use_imputed_full_data':False}


	for opt, arg in opts:
		if opt in ("-h"):
			print __doc__
			return
		elif opt in ('-i'): p_dict['pids'] = _parse_pids_(arg)
		elif opt in ('-o'): p_dict['run_id'] = arg
		elif opt in ('-p'): p_dict['parallel'] = arg
		elif opt in ('-t'): p_dict['call_method_id'] = int(arg)
		elif opt in ('-n'): p_dict['add_to_db'] = True
		elif opt in ('-m'): p_dict['mac_threshold'] = int(arg)
		elif opt in ('-f'): p_dict['data_file'] = arg
		elif opt in ('-r'): p_dict['phen_file'] = arg
		elif opt in ('-k'): p_dict['kinship_file'] = arg
		elif opt in ('-a'): p_dict['specific_methods'] = arg.split(',')
		elif opt in ("-b"): p_dict['specific_transformations'] = arg.split(',')
		elif opt in ("-c"): p_dict['remove_outliers'] = int(arg)
		elif opt in ("-d"): p_dict['debug_filter'] = float(arg)
		elif opt in ('-u'): p_dict['use_existing_results'] = True
		elif opt in ('-e'): p_dict['analysis_plots'] = True
		elif opt in ('-q'): p_dict['mem_req'] = arg
		elif opt in ('-l'): p_dict['walltime_req'] = arg
		elif opt in ("--comment"): p_dict['comment'] = arg
		elif opt in ("--no_phenotype_ids"): p_dict['no_phenotype_ids'] = True
		elif opt in ("--region_plots"): p_dict['region_plots'] = int(arg)
		elif opt in ("--cand_genes_file"): p_dict['cand_genes_file'] = arg
		elif opt in ("--proc_per_node"): p_dict['proc_per_node'] = int(arg)
		elif opt in ("--only_add_2_db"): p_dict['only_add_2_db'] = True
		elif opt in ("--data_format"): p_dict['data_format'] = arg
		elif opt in ("--emmax_perm"): p_dict['emmax_perm'] = int(arg)
		elif opt in ("--with_replicates"): p_dict['with_replicates'] = True
		elif opt in ("--with_betas"): p_dict['with_betas'] = True
		elif opt in ("--num_steps"): p_dict['num_steps'] = int(arg)
		elif opt in ("--local_gwas"): p_dict['local_gwas'] = map(int, arg.split(','))
		elif opt in ("--use_imputed_full_data"): p_dict['use_imputed_full_data'] = True
		else:
			print "Unkown option:", opt
			print __doc__
			sys.exit(2)
	return p_dict, args



def _prepare_transformation_(phed, pid, trans_type, remove_outliers):
	num_outliers_removed = 0
	phed.transform(pid, trans_type)
	if remove_outliers:
		print 'Removing outliers above IQR fence of', remove_outliers, '..'
		num_outliers_removed = phed.na_outliers(pid, iqr_threshold=remove_outliers)
	return num_outliers_removed



def prepare_data(sd, phed, pid, trans_type, remove_outliers, with_replicates=False):
	"""
	Coordinates phenotype and snps data for different mapping methods.
	"""
	if not with_replicates:
		print 'Converting replicates of phenotypes to averages'
		phed.convert_to_averages([pid])
	num_outliers_removed = _prepare_transformation_(phed, pid, trans_type, remove_outliers)
	sd.coordinate_w_phenotype_data(phed, pid)
	return num_outliers_removed



def get_perm_pvals(snps, phen_vals, mapping_method='kw', num_perm=100, snps_filter=0.05):
	import random
	if snps_filter < 1.0:
		snps = random.sample(snps, int(snps_filter * len(snps)))
	pvals = []
	if mapping_method == 'kw':
		for i in range(num_perm):
			random.shuffle(phen_vals)
			kw_res = util.kruskal_wallis(snps, phen_vals, verbose=False)
			pvals.extend(kw_res['ps'])
	elif mapping_method == 'ft':
		for i in range(num_perm):
			random.shuffle(phen_vals)
			pvals.extend(run_fet(snps, phen_vals))
	return pvals



def _get_file_prefix_(id, p_i, phenotype_name, mapping_method=None, trans_method=None, remove_outliers=None,
		with_replicates=False):
	prefix = env['results_dir'] + id + "_pid" + str(p_i) + "_" + phenotype_name
	if mapping_method:
		prefix += "_" + mapping_method
	if trans_method:
		prefix += "_" + trans_method
	if remove_outliers:
		prefix += '_no' + str(remove_outliers)
	if with_replicates:
		prefix += '_with_replicates'
	return prefix



def _get_genotype_data_(p_dict):
	if p_dict['use_imputed_full_data']:
		if p_dict['local_gwas']:
			chrom = p_dict['local_gwas'][0]
			sd = dataParsers.load_1001_full_snps(debug_filter=p_dict['debug_filter'], chromosomes=[chrom])
		else:
			sd = dataParsers.load_1001_full_snps(debug_filter=p_dict['debug_filter'])
	else:
		sd = dataParsers.parse_snp_data(p_dict['data_file'] , format=p_dict['data_format'], filter=p_dict['debug_filter'])
	return sd



def run_parallel(p_i, phed, p_dict, mapping_method="analysis", trans_method='none'):
	"""
	If no mapping_method, then analysis run is set up.
	"""

	if not p_i in phed.phen_dict:
		print "Phenotype ID not found:%i" % p_i
		return
	if phed.is_binary(p_i):
		if trans_method != 'none':
			return
		elif mapping_method in ["kw"]:
			return
		elif p_dict['analysis_plots']:
			specific_methods = ["ft", 'lm', "emma", 'emmax']
	else:
		if mapping_method in ["ft"]:
			return
		elif p_dict['analysis_plots']:
			specific_methods = ["kw", 'lm', "emma", 'emmax']

	phenotype_name = phed.get_name(p_i)
	#Cluster specific parameters
	print "Setting up a gwa run for phenotype:%s, pid=%d, using method:%s, with transformation as:%s"\
		% (phenotype_name, p_i, mapping_method, trans_method)
	run_id = p_dict['parallel']

	#Add results to DB (a hack for papaya, to add results from main node to DB).
	if p_dict['only_add_2_db']:
		file_name_prefix = _get_file_prefix_(p_dict['parallel'], p_i, phenotype_name, mapping_method, \
						trans_method, p_dict['remove_outliers'], p_dict['with_replicates'])
		pval_file = file_name_prefix + ".pvals"
		score_file = file_name_prefix + ".scores"
		sys.stdout.write("Looking for files %s or %s." % (pval_file, score_file))
		result_file = None
		if os.path.isfile(pval_file):
			result_file = pval_file
		elif os.path.isfile(score_file):
			result_file = score_file
		if result_file:
			sys.stdout.write("..... found!\n")
			if p_dict['no_phenotype_ids']:
				db_pid = phed.get_db_pid(p_i)
			else:
				db_pid = p_i

			import results_2_db as rdb
			short_name = "cm" + str(p_dict['call_method_id']) + "_pid" + str(db_pid) + "_" + phenotype_name \
				+ "_" + mapping_method + "_" + trans_method
			if p_dict['remove_outliers']:
				short_name += "_no"
			tm_id = transformation_method_dict[trans_method]
			rdb.add_results_to_db(result_file, short_name, p_dict['call_method_id'], db_pid, \
					analysis_methods_dict[mapping_method], tm_id, \
					remove_outliers=p_dict['remove_outliers'])
			return
		else:
			sys.stdout.write("Result files not found!\n")
			sys.stdout.write("Setting up the run.\n")
			sys.stdout.flush()



	shstr = "#!/bin/csh\n"
	shstr += "#PBS -l walltime=%s \n" % p_dict['wall_time_req']
	if mapping_method in ['emma']:
		shstr += "#PBS -l nodes=1:ppn=%d \n" % p_dict['proc_per_node']
	shstr += "#PBS -l mem=%s \n" % p_dict['mem_req']
	shstr += "#PBS -q cmb\n"

	job_id = '%d_%s_%s_%s' % (p_i, mapping_method, trans_method, p_dict['parallel'])
	shstr += "#PBS -N p%s \n" % job_id
	shstr += "(python %s gwa.py -o %s " % (env['script_dir'], run_id)
	shstr += "-i %d -a %s -b %s -c %d -d %n -t %d -m %d " % (p_i, p_dict['mapping_method'], p_dict['trans_method'],
					p_dict['remove_outliers'], p_dict['debug_filter'],
					p_dict['call_method_id'], p_dict['mac_threshold'])

	shstr += "--region_plots=%d --proc_per_node=%d --data_format=%s " % \
		(p_dict['region_plots'], p_dict['proc_per_node'], p_dict['data_format'])


	if p_dict['use_existing_results']: shstr += "--use_existing_results  "
	if p_dict['analysis_plots']: shstr += "-e "
	if p_dict['phen_file']: shstr += "-r %s " % p_dict['phen_file']
	if p_dict['kinship_file']: shstr += "-k %s " % p_dict['kinship_file']
	if p_dict['cand_genes_file']: shstr += "--cand_genes_file=%s " % p_dict['cand_genes_file']
	if p_dict['no_phenotype_ids']: shstr += "--no_phenotype_ids "
	if p_dict['add_to_db']: shstr += "-n "
	if p_dict['comment']: shstr += "--comment=%s " % comment

	shstr += "> " + run_id + "_job" + ".out) >& " + run_id + "_job" + ".err\n"
	#print '\n',shstr,'\n'
	script_file_name = p_dict['parallel'] + ".sh"
	f = open(script_file_name, 'w')
	f.write(shstr)
	f.close()

	#Execute qsub script
	os.system("qsub " + script_file_name)



def analysis_plots(phed, p_dict):
	print "\nAnalysing GWAs results jointly... QQ plots etc."

	#Genotype and phenotype data is only used for permutations.
	sd = _get_genotype_data_(p_dict)

	#try:
	print "Plotting accession phenotype map"
	for p_i in p_dict['pids']:
		sys.stdout.flush()
		file_prefix = _get_file_prefix_(p_dict['run_id'], p_i, phed.get_name(p_i))
		accession_map_pdf_file = file_prefix + "_acc_phen_map.pdf"
		phed.plot_accession_map(p_i, pdf_file=accession_map_pdf_file)
	#except Exception, err_str:
	#	print 'Skipping accession - phenotype map... basemap package is probably not installed:', err_str

	#Load gwas results
	results = {}
	perm_pvals = None
	method_types = []
	#Iterate over all possible combination of results
	for p_i in p_dict['pids']:
		phenotype_name = phed.get_name(p_i)
		phen_is_binary = phed.is_binary(p_i)
		print "Plotting analysis plots for phenotype:%s, phenotype_id:%s" % (phenotype_name, p_i)
		for trans_method in p_dict['specific_transformations']:
			prepare_data(sd, phed, p_i, trans_method, p_dict['remove_outliers'], p_dict['with_replicates'])
			for mapping_method in p_dict['specific_methods']:
				file_prefix = _get_file_prefix_(p_dict['run_id'], p_i, phed.get_name(p_i),
							mapping_method, trans_method, p_dict['remove_outliers'])
				for with_reps in ['', '_with_replicates']:
					file_prefix = file_prefix + with_reps
					mapping_label = mapping_method + '_' + trans_method + with_reps
					snps = sd.getSnps()
					phen_vals = phed.get_values(p_i)
					res_name = "%s_%s_%s" % (phenotype_name, mapping_method, trans_method)
					try:
						file_name = file_prefix + ".pvals"
						sys.stdout.write("Looking for file %s." % (file_name))
						if os.path.isfile(file_name):
							sys.stdout.write("..... found!\n")
							res = gwaResults.Result(result_file=file_name, name=res_name, snps=snps)
							pvals = True
							if mapping_method in ['lm', 'emma', 'emmax']:
								res.filter_attr("macs", p_dict['mac_threshold'])
							results[mapping_label] = res
							method_types.append(mapping_method)

						else:
							sys.stdout.write("..... not found.\n")
							sys.stdout.flush()
							file_name = file_prefix + ".scores"
							sys.stdout.write("Looking for file %s." % (file_name))
							if os.path.isfile(file_name):
								sys.stdout.write("..... found!\n")
								res = gwaResults.Result(result_file=file_name, name=res_name,
										snps=snps)
								pvals = False
								results[mapping_label] = res
								method_types.append(mapping_method)
							else:
								sys.stdout.write("..... not found.\n")
						sys.stdout.flush()


						#Permutation tests for KW and FT..
						if mapping_method in ['kw', 'ft']:
							print "Retrieving permutations for %s" % (mapping_method)
							sys.stdout.flush()
							perm_pvals = get_perm_pvals(snps, phen_vals, mapping_method)
					except Exception, err_str:
						print err_str
						print "Failed loading file %s" % (file_name)
						sys.stdout.flush()
		if len(results.keys()) > 0:
			print "Drawing QQ plots for methods: %s" % (str(results.keys()))
			for k in results:
				print k, results[k].name
			sys.stdout.flush()

                        #Plotting the QQ plots.
			num_dots = 1000
			max_log_val = 6
			qq_file_prefix = _get_file_prefix_(p_dict['run_id'], p_i, phed.get_name(p_i))
			gwaResults.qq_plots(results, num_dots, max_log_val, qq_file_prefix, method_types=method_types,
					mapping_labels=results.keys(), phen_name=phenotype_name, perm_pvalues=perm_pvals,
					is_binary=phen_is_binary)



def map_phenotype(p_i, phed, mapping_method, trans_method, p_dict):
	phed = copy.deepcopy(phed)
	phenotype_name = phed.get_name(p_i)
	phen_is_binary = phed.is_binary(p_i)
	file_prefix = _get_file_prefix_(p_dict['run_id'], p_i, phed.get_name(p_i),
				mapping_method, trans_method, p_dict['remove_outliers'], p_dict['with_replicates'])
	result_name = "%s_%s_%s" % (phenotype_name, mapping_method, trans_method)
	emmax_perm_threshold = None
	k = None

	res = None
	#Check whether result already exists.
	if p_dict['use_existing_results']:
		if p_dict['region_plots']:
			sd = _get_genotype_data_(p_dict)
			num_outliers = prepare_data(sd, phed, p_i, trans_method, p_dict['remove_outliers'], p_dict['with_replicates'])
			if p_dict['remove_outliers']:
				assert num_outliers != 0, "No outliers were removed, so it makes no sense to go on and perform GWA."

			snps = sd.getSnps()
		else:
			snps = None

		print "\nChecking for existing results."
		result_file = file_prefix + ".pvals"
		if os.path.isfile(result_file):
			res = gwaResults.Result(result_file=result_file, name=result_name, snps=snps)
			pvals = True
		else:
			result_file = file_prefix + ".scores"
			if os.path.isfile(result_file):
				res = gwaResults.Result(result_file=result_file, name=result_name, snps=snps)
				pvals = False
		if res:
			print "Found existing results.. (%s)" % (result_file)

		sys.stdout.flush()


	if not res: #If results weren't found in a file... then do GWA.
		#Loading data
		sd = _get_genotype_data_(p_dict)
		#Do we need to calculate the K-matrix?
		if mapping_method in ['emma', 'emmax', 'emmax_anova', 'emmax_step']:
			#Load genotype file (in binary format)
			sys.stdout.write("Retrieving the Kinship matrix K.\n")
			sys.stdout.flush()
			if p_dict['use_imputed_full_data']:
				k_file = env['data_1001_dir'] + 'kinship_matrix.pickled'
			else:
				k_file = env['data_dir'] + "kinship_matrix_cm" + str(p_dict['call_method_id']) + ".pickled"
			kinship_file = p_dict['kinship_file']
			if not kinship_file and os.path.isfile(k_file): #Check if corresponding call_method_file is available
				kinship_file = k_file
			if kinship_file:   #Kinship file was somehow supplied..
				num_outliers = prepare_data(sd, phed, p_i, trans_method, p_dict['remove_outliers'],
							p_dict['with_replicates'])
				print 'Loading supplied kinship'
				k = lm.load_kinship_from_file(kinship_file, sd.accessions)
			else:
				print "No kinship file was found.  Generating kinship file:", k_file
				k_accessions = sd.accessions[:]
				k = sd.get_ibs_kinship_matrix(p_dict['debug_filter'])
				f = open(k_file, 'w')
				cPickle.dump([k, k_accessions], f)
				f.close()
				num_outliers = prepare_data(sd, phed, p_i, trans_method, p_dict['remove_outliers'],
							p_dict['with_replicates'])
				k = lm.filter_k_for_accessions(k, k_accessions, sd.accessions)
			sys.stdout.flush()
			sys.stdout.write("Done!\n")
		else:
			num_outliers = prepare_data(sd, phed, p_i, trans_method, p_dict['remove_outliers'],
						p_dict['with_replicates'])

		if p_dict['remove_outliers']:
			assert num_outliers != 0, "No outliers were removed, so it makes no sense to go on and perform GWA."
		phen_vals = phed.get_values(p_i)

		if p_dict['local_gwas']: #Filter SNPs, etc..
			sd = snpsdata.SNPsDataSet([sd.get_region_snpsd(*p_dict['local_gwas'])],
						[p_dict['local_gwas'][0]])
		snps = sd.getSnps()

		#Loading candidate genes
		cand_genes = None
		if p_dict['cand_genes_file']:
			cand_genes, tair_ids = gwaResults.load_cand_genes_file(p_dict['cand_genes_file'])
		else:
			cand_genes = None
			tair_ids = None

		sys.stdout.write("Finished loading and handling data!\n")

		print "Plotting a histogram"
		p_her = None
		if k is not None:
			p_her = phed.get_pseudo_heritability(p_i, k)
		hist_file_prefix = _get_file_prefix_(p_dict['run_id'], p_i, phenotype_name, trans_method,
						p_dict['remove_outliers'], p_dict['with_replicates'])
		hist_png_file = hist_file_prefix + "_hist.png"
		phed.plot_histogram(p_i, png_file=hist_png_file, p_her=p_her)


		print "Applying %s to data." % (mapping_method)
		sys.stdout.flush()
		kwargs = {}
		additional_columns = []
		if "kw" == mapping_method:

			if phen_is_binary:
				warnings.warn("Warning, applying KW to a binary phenotype")

			kw_res = util.kruskal_wallis(snps, phen_vals)
			pvals = kw_res['ps']
			kwargs['statistics'] = kw_res['ds']
			additional_columns.append('statistics')


		elif "ft" == mapping_method:
			pvals, or_est = run_fet(snps, phen_vals)
			kwargs['odds_ratio_est'] = or_est
			additional_columns.append('odds_ratio_est')

		else:  #Parametric tests below:		

			if mapping_method in ['emma']:
				res = lm.emma(snps, phen_vals, k)
			elif mapping_method in ['emmax']:
				if p_dict['emmax_perm']:
					perm_sd = _get_genotype_data_(p_dict)
					num_outliers = prepare_data(perm_sd, phed, p_i, 'none', 0, p_dict['with_replicates'])
					perm_sd.filter_mac_snps(p_dict['mac_threshold'])
					t_snps = perm_sd.getSnps()
					t_phen_vals = phed.get_values(p_i)
					res = lm.emmax_perm_test(t_snps, t_phen_vals, k, p_dict['emmax_perm'])
					emmax_perm_threshold = res['threshold_05'][0]
					import pylab
					hist_res = pylab.hist(-sp.log10(res['min_ps']))
					threshold = -sp.log10(emmax_perm_threshold)
					b_threshold = -sp.log10(1.0 / (len(t_snps) * 20.0))
					pylab.vlines(threshold, 0, max(hist_res[0]), color='g')
					pylab.vlines(b_threshold, 0, max(hist_res[0]), color='r')
					pylab.savefig(env['tmp_dir'] + 'test.png', format='png')
				if p_dict['with_replicates']:
					#Get values, with ecotypes, construct Z and do GWAM
					phen_vals = phed.get_values(p_i)
					Z = phed.get_incidence_matrix(p_i)
					res = lm.emmax(snps, phen_vals, k, Z=Z, with_betas=p_dict['with_betas'])
				else:
					res = lm.emmax(snps, phen_vals, k, with_betas=p_dict['with_betas'])

			elif mapping_method in ['emmax_step']:
				sd.filter_mac_snps(p_dict['mac_threshold'])
				local = False
				if p_dict['local_gwas']:
					local = True
					file_prefix += '_' + '_'.join(map(str, p_dict['local_gwas']))
				res = lm.emmax_step_wise(phen_vals, k, sd=sd, num_steps=p_dict['num_steps'],
							file_prefix=file_prefix, local=local, cand_gene_list=cand_genes)
				print 'Step-wise EMMAX finished!'
				return
			elif mapping_method in ['lm']:
				res = lm.linear_model(snps, phen_vals)
			elif mapping_method in ['emmax_anova']:
				res = lm.emmax_anova(snps, phen_vals, k)
			elif mapping_method in ['lm_anova']:
				res = lm.anova(snps, phen_vals)
			else:
				print "Mapping method", mapping_method, 'was not found.'
				return

			if mapping_method in ['lm', 'emma', 'emmax']:
				kwargs['genotype_var_perc'] = res['var_perc']
				additional_columns.append('genotype_var_perc')
				if p_dict['with_betas'] or mapping_method not in ['emmax', 'lm']:
					betas = map(list, zip(*res['betas']))
					kwargs['beta0'] = betas[0]
					additional_columns.append('beta0')
					if len(betas) > 1:
						kwargs['beta1'] = betas[1]
						additional_columns.append('beta1')
				pvals = res['ps']
				sys.stdout.write("Done!\n")
				sys.stdout.flush()

			if mapping_method in ['lm_anova', 'emmax_anova']:
				kwargs['genotype_var_perc'] = res['var_perc']
				pvals = res['ps']
				sys.stdout.write("Done!\n")
				sys.stdout.flush()


#		print 'Calculating SNP-phenotype correlations.'
#		kwargs['correlations'] = calc_correlations(snps, phen_vals)
#		additional_columns.append('correlations')
		print 'Writing result to file.'
		res = gwaResults.Result(scores=pvals.tolist(), snps_data=sd, name=result_name, **kwargs)
		if mapping_method in ["kw", "ft", "emma", 'lm', "emmax", 'emmax_anova', 'lm_anova']:
		 	result_file = file_prefix + ".pvals"
		else:
		 	result_file = file_prefix + ".scores"
		res.write_to_file(result_file, additional_columns)

	#add results to DB..

	if p_dict['add_to_db']:
		print 'Adding results to DB.'
		if p_dict['no_phenotype_ids']:
			db_pid = phed.get_db_pid(p_i)
		else:
			db_pid = p_i

		import results_2_db as rdb

		short_name = 'cm%d_pid%d_%s_%s_%s_%d_%s' % (p_dict['call_method_id'], db_pid, phenotype_name,
							mapping_method, trans_method, p_dict['remove_outliers'],
							str(p_dict['with_replicates']))
		tm_id = transformation_method_dict[trans_method]
		rdb.add_results_to_db(result_file, short_name, p_dict['call_method_id'], db_pid,
					analysis_methods_dict[mapping_method],
					tm_id, remove_outliers=p_dict['remove_outliers'])



	if p_dict['data_format'] != 'float':
		#Load candidate genes from a file, if it is given

		print "Generating a GW plot."
		sys.stdout.flush()
		if p_dict['local_gwas']:
			png_file = file_prefix + "_gwa_plot_" + '_'.join(map(str, p_dict['local_gwas'])) + "_.png"
		else:
			png_file = file_prefix + "_gwa_plot.png"
		#png_file_max30 = file_prefix+"_gwa_plot_max30.png"
		if mapping_method in ["kw", "ft", "emma", 'lm', "emmax", 'emmax_anova', 'lm_anova']:
			res.neg_log_trans()
			if mapping_method in ["kw", "ft"]:# or p_dict['data_format'] != 'binary':
				#res.plot_manhattan(png_file=png_file_max30,percentile=90,type="pvals",ylab="$-$log$_{10}(p)$", 
				#	       plot_bonferroni=True,cand_genes=cand_genes,max_score=30)
				res.plot_manhattan(png_file=png_file, percentile=90, type="pvals",
						ylab="$-$log$_{10}(p)$", plot_bonferroni=True,
						cand_genes=cand_genes, threshold=emmax_perm_threshold)
			else:
				if res.filter_attr("macs", p_dict['mac_threshold']) > 0:
					#res.plot_manhattan(png_file=png_file_max30,percentile=90,type="pvals",ylab="$-$log$_{10}(p)$", 
					#	       plot_bonferroni=True,cand_genes=cand_genes,max_score=30)				
					res.plot_manhattan(png_file=png_file, percentile=90, type="pvals",
							ylab="$-$log$_{10}(p)$", plot_bonferroni=True,
							cand_genes=cand_genes, threshold=emmax_perm_threshold)
		else:
			pass




		if p_dict['region_plots']:
			import regionPlotter as rp
			regions_results = res.get_top_region_results(p_dict['region_plots'], distance_threshold=60000, buffer=60000)
			#plotter = rp.RegionPlotter()
			print "Starting region plots..."
			for reg_res in regions_results:
				chromosome = reg_res.snp_results['chromosomes'][0]
				caption = phenotype_name + "_c" + str(chromosome) + "_" + mapping_method
				res_positions = reg_res.snp_results['positions']

				png_file = '%s_reg_plot_c%d_s%d_e%d.png' % \
						(file_prefix , chromosome, res_positions[0], res_positions[-1])
				tair_file = '%s_reg_plot_c%d_s%d_e%d_tair_info.txt' % \
						(file_prefix , chromosome, res_positions[0], res_positions[-1])
#				plotter.plot_small_result([reg_res], png_file=png_file, highlight_gene_ids=tair_ids,
#							  caption=caption, tair_file=tair_file)
				reg_res.plot_manhattan(png_file=png_file, percentile=0, tair_file=tair_file)

				#Plot Box-plot
				box_plot_file = '%s_reg_plot_c%d_s%d_e%d_box_plot.png' % \
						(file_prefix, chromosome, res_positions[0], res_positions[-1])
				(marker, score, chromosome, pos) = reg_res.get_max_snp()
				m_accessions = sd.accessions
				phed.revert_to_raw_values(p_i)
				phed.plot_marker_box_plot(p_i, marker, m_accessions, plot_file=box_plot_file,
							m_chromosome=chromosome, m_position=pos, m_score=score)
				acc_hist_plot_file = '%s_reg_plot_c%d_s%d_e%d_acc_hist_plot.png' % \
						(file_prefix, chromosome, res_positions[0], res_positions[-1])
				phed.plot_marker_accessions_hist(p_i, marker, m_accessions, plot_file=acc_hist_plot_file,
							m_chromosome=chromosome, m_position=pos, m_score=score)




def _run_():

	p_dict, args = parse_parameters()
	print "GWA runs are being set up with the following parameters:"
	for k, v in p_dict.iteritems(): print k + ': ' + str(v)
	print ''

	#Load phenotype file
	if p_dict['phen_file']:
		print 'Loading phenotypes from file.'
		phed = phenotypeData.parse_phenotype_file(p_dict['phen_file'], with_db_ids=(not p_dict['no_phenotype_ids']))  #load phenotype file
	else:
		print 'Retrieving the phenotypes from the DB.'
		phed = phenotypeData.get_all_phenotypes_from_db()

	#If on the cluster, then set up runs..
	if p_dict['parallel']:
		if len(p_dict['pids']) == 0:  #phenotype index arguement is missing, hence all phenotypes are run/analyzed.
			if not p_dict['phen_file']:
				raise Exception('Phenotype file or phenotype ID is missing.')
			p_dict['pids'] = phed.phen_dict.keys()
		else:
			raise Exception('Too many arguments..')

		if analysis_plots:  #Running on the cluster..
			for p_i in p_dict['pids']:
				run_parallel(p_i, phed, p_dict)
		else:
			for mapping_method in p_dict['specific_methods']:
				for trans_method in p_dict['specific_transformations']:
					for p_i in pids:
						run_parallel(p_i, phed, p_dict, mapping_method, trans_method)
		return #Exiting the program...


	#SNPs data file name
	if not p_dict['data_file']:
		p_dict['data_file'] = '%s250K_t%d.csv' % (env['data_dir'], p_dict['call_method_id'])


	#Plot analysis plots...
	if p_dict['analysis_plots']:
		analysis_plots(phed, p_dict)
	else:
		#If not analysis plots... then GWAS
		for p_i in p_dict['pids']:
			if p_i in phed.phen_dict:
				print '-' * 120, '\n'
				phenotype_name = phed.get_name(p_i)
				print "Performing GWAS for phenotype: %s, phenotype_id: %s" % (phenotype_name, p_i)
				for trans_method in p_dict['specific_transformations']:
					print 'Phenotype transformation:', trans_method

					for mapping_method in p_dict['specific_methods']:
						#DO ANALYSIS
						print 'Mapping method:', mapping_method
						map_phenotype(p_i, phed, mapping_method, trans_method, p_dict)




#def _run_emma_mp_(in_que, out_que, confirm_que, num_splits=10):
#	from rpy import r
#	args = in_que.get(block=True)
#	#r_source(script_dir+"emma_fast.R")
#	#r_emma_REML_t = robjects.r['emma.REML.t']
#	r.source(env['script_dir'] + "emma_fast.R")
#	phenArray = array([args[2]])
#	snps = args[1]
#	print "Found data nr. %i, starting run.  pid=%d" % (args[0], os.getpid())
#	sys.stdout.flush()
#	confirm_que.put([args[0], "start confirmation!"], block=False)
#	results = []
#	for i in range(num_splits):
#		snp_slice = snps[(len(snps) * i) / num_splits:(len(snps) * (i + 1)) / num_splits]
#		snpsArray = array(snp_slice)
#	        results.append(r.emma_REML_t(phenArray, snpsArray, args[3], [args[4]]))
#		#results.append(r_emma_REML_t(phenArray,snpsArray,args[3],[args[4]]))
#		if args[0] == 0:
#			print "About %d.0%% is done now." % (100 * (i + 1.0) / num_splits)
#			sys.stdout.flush()
#	new_res = {}
#	for k in results[0]:
#		new_res[k] = []
#	for nr in results:
#		for k in nr:
#			new_res[k].extend(list(nr[k]))
#	out_que.put([args[0], new_res])
#
	#def run_emma_parallel(snps, phen_vals, k, num_proc=8):
	#	"""
	#	Run EMMA on multiple processors, using multiprocessing
	#	"""
	#	processes = []
	#	in_que = mp.Queue()
	#	out_que = mp.Queue()
	#	confirm_que = mp.Queue()
	#	import scipy
	#	phen_var = scipy.var(phen_vals)
	#
	#	#Populating the in_que
	#	for i in range(num_proc):
	#		snp_slice = snps[(len(snps) * i) / num_proc:(len(snps) * (i + 1)) / num_proc]
	#		print len(snp_slice)
	#		in_que.put([i, snp_slice, phen_vals, k, phen_var], block=False)
	#		p = mp.Process(target=_run_emma_mp_, args=(in_que, out_que, confirm_que))
	#		p.start()
	#		done = False
	#		fail_count = 0
	#		while not done:
	#			try:
	#				confirm = confirm_que.get(block=True, timeout=90)
	#				done = True
	#			except:
	#				p.terminate()
	#				fail_count += 1
	#				if fail_count > 10:
	#					print "Giving up!!"
	#					for p in processes:
	#						p.terminate()
	#					raise Exception("Failed to start all processes.")
	#				else:
	#					print "Trying again with a new process!"
	#					in_que.put([i, snp_slice, phen_vals, k, phen_var], block=False)
	#					p = mp.Process(target=_run_emma_mp_, args=(in_que, out_que, confirm_que))
	#					p.start()
	#
	#
	#		print "ID=%i, Recieved %s" % (confirm[0], confirm[1])
	#		sys.stdout.flush()
	#		processes.append(p)
	#
	#	for p in processes:
	#		print p, p.is_alive()
	#
	#	results = []
	#	for i in range(num_proc):
	#		#import pdb;pdb.set_trace()
	#		if i > 0:
	#			try:
	#				res = out_que.get(block=True, timeout=10000) #waits about 3 hours, if needed. 
	#				results.append(res)
	#			except Exception, err_str:
	#				print "The parent process didn't receive the results, after waiting over almost 3 hours."
	#
	#		else:
	#			res = out_que.get()
	#			results.append(res)
	#
	#	for p in processes:
	#		p.terminate()
	#
	#	results.sort()
	#	new_res = {}
	#	for k in results[0][1]:
	#		new_res[k] = []
	#	for nr in results:
	#		for k in nr[1]:
	#			new_res[k].extend(list(nr[1][k]))
	#
	#	#import pdb;pdb.set_trace()
	#	return new_res
	#	#print "len(q):",len(q)



#def runEmma(snps, phenValues, k):
#	from rpy import r
#	#Assume that the accessions are ordered.
#	r.source("emma.R")
#	#r_emma_REML_t = robjects.r['emma.REML.t']
#
#	phenArray = array([phenValues])
#	snpsArray = array(snps)
#	res = r.emma_REML_t(phenArray, snpsArray, k)
#	return res


def run_fet(snps, phenotypeValues, verbose=False):
	from rpy import r
	"""
	Fisher's exact test.
	"""
	print "Running Fisher's exact test on", len(snps), "snps, and", len(phenotypeValues), "phenotype values."
	pvals = []
	or_est = []  #Odds ratio estimates
	for snp in snps:
		res = r.fisher_test(phenotypeValues, snp)
		pvals.append(res["p.value"])
		or_est.append(res['estimate'])
	return pvals, or_est


def calc_correlations(snps, phen_vals):
	import scipy as sp
	corrs = sp.zeros(len(snps))
	for i, snp in enumerate(snps):
		corrs[i] = sp.corrcoef(snp, phen_vals)[1, 0]
	return corrs



if __name__ == '__main__':
	_run_()
	print "Done!"
