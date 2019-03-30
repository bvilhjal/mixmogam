""" 
This library offers functions to parse different types of SNPs data from multiple formats into a snpsData objects.

Author: Bjarni J. Vilhjalmsson
Email: bjarni.vilhjalmsson@gmail.com
"""
import time, random
import cPickle
from snpsdata import SNPsData, RawSnpsData, SNPsDataSet
import scipy as sp
import h5py
import numpy as np
import os
import sys

# this should be fixed

# Standard missing value is N (Used to be NA)
missing_val = 'N'
# Standard nt decoding, is using the IUPAC alphabet
nt_decoder = {'A':'A',
          'C':'C',
          'G':'G',
          'T':'T',
          'AG':'R',
          'AC':'M',
          'GT':'K',
          'CT':'Y',
          'AT':'W',
          'CG':'S',
          'Y':'Y',
          'R':'R',
          'W':'W',
          'S':'S',
          'K':'K',
          'M':'M',
          'D':'D',
          'H':'H',
          'V':'V',
          'B':'B',
          'X':'X',  # Unknown base(s)
          'N':'X',  # Unknown base(s)
          '-':'-',  # Indel 
          '|':missing_val}

# An int decoder is useful for processing the data efficiently
nt_int_decoder = {'A':1,
          'C':2,
          'G':3,
          'T':4,
          'AG':5,
          'AC':6,
          'GT':7,
          'CT':8,
          'AT':9,
          'CG':10,
          'Y':11,
          'R':12,
          'W':13,
          'S':14,
          'K':15,
          'M':16,
          'D':17,
          'H':18,
          'V':19,
          'B':20,
          'X':21,  # Unknown base(s)
          'N':21,  # Unknown base(s)
          '-':22,  # Indel 
          '|':0}


def parse_raw_snps_data(datafile, target_format='nucleotides', deliminator=",", missing_val='N', return_chromosomes=False,
        use_decoder=True, debug_filter=1.0, marker_type=None, verbose=True):
    """
    Parses nucleotide SNPs data files into a RawSnpsData.
    """

    sys.stderr.write("Loading file: %s ... \n" % datafile)
    decoder = nt_decoder
    decoder[missing_val] = missing_val
    #code_fun = sp.vectorize(lambda x: decoder[x], otypes=['a1']) #FINISH!!!
    accessions = []

    # Reading column data
    with open(datafile, 'r') as f:
        array_ids = None
        first_line = map(str.strip, f.next().split(deliminator))
        second_line = map(str.strip, f.next().split(deliminator))
        if first_line[0] not in ["Chromosome", 'chr'] and second_line in ["Chromosome", 'chr']:
            array_ids = first_line[2:]
            accessions = second_line[2:]
        elif first_line[0] in ["Chromosome", 'chr']:
            accessions = first_line[2:]
        else:
            raise Exception('Unknown file format')

        num_accessions = len(accessions)
        print "Found", num_accessions, "arrays/strains."

        pos_snps_dict = {}
        num_snps = 0
        num_snps_w_3_alleles = 0
        d = {'snps':[], 'positions':[]}
        for snp_i, line in enumerate(f):
            if verbose and snp_i % 100000 == 0:
                print '%d SNPs have been read.' % snp_i
                print '%d SNPs with issues (three alleles etc.) have been ignored' % num_snps_w_3_alleles
            if random.random() >= debug_filter:
                continue
            num_snps += 1
            l = line.strip()
            l = l.split(deliminator)
            chrom = int(l[0])
            pos = int(l[1])
            if use_decoder:
                snp = sp.empty(num_accessions, 'a1')
#                s = l[2:]
#                if '0' in s or '' in s:
#                    num_snps_w_3_alleles += 1 #A hack
#                    continue
                for i in xrange(num_accessions):
                    snp[i] = decoder[l[2 + i]]
            else:
                snp = sp.empty(num_accessions, 'a')
                for i in xrange(num_accessions):
                    snp[i] = l[2 + i]

            try:
                d = pos_snps_dict[chrom]
            except KeyError:
                d = {'snps':[], 'positions':[]}
                pos_snps_dict[chrom] = d
                if verbose and snp_i:
                    print '%d SNPs have been loaded.' % num_snps

            d['snps'].append(snp)
            d['positions'].append(pos)

    snps_data_list = []
    chromosomes = sorted(pos_snps_dict.keys())
    for chrom in chromosomes:
        snps_data_list.append(RawSnpsData(accessions=accessions, positions=pos_snps_dict[chrom]['positions'], \
                snps=pos_snps_dict[chrom]['snps'], chromosome=chrom, arrayIds=array_ids))

    if target_format == 'binary':
        print "Converting raw SNPs data to binary SNPs."
        for i in range(len(chromosomes)):
            snps_data_list[i] = snps_data_list[i].getSnpsData()
    if return_chromosomes:
        return (snps_data_list, chromosomes)
    else:
        return snps_data_list


def parse_numerical_snp_data(data_file, delimiter=",", missing_val='NA', filter=1,
                             filter_accessions=None, use_pickle=True, dtype='int8',
                             data_format='binary'):
    """
    A sped-up version, to load a int (e.g. binary) file.
    
    If pickle is used then more memory is required, but it speeds up.
    """
    pickle_file = data_file + '.pickled'
    filter_accessions_ = None
    if use_pickle:
        if os.path.isfile(pickle_file):
            t_ = time.time()
            f = open(pickle_file, 'rb')
            sd = cPickle.load(f)
            f.close()
            if filter_accessions and set(filter_accessions) != set(sd.accessions):
                print 'Filtering accessions.'
                sd.filter_accessions(filter_accessions)
            if filter < 1:
                sd.sample_snps(filter)
                print 'Filtering random %d%%' % int(filter * 100)
            print 'Loading genotype data took %.2f s...' % (time.time() - t_)
            return sd
        else:
            print 'Pickle file not found: %s' % pickle_file
            filter_accessions_ = filter_accessions
            filter_accessions = None
    sys.stderr.write("Loading the SNPs data file: %s \n" % data_file)
    accessions = []
    snpsd_ls = []
    chromosomes = []

    # Reading accession data
    f = open(data_file, 'rb')
    line = f.next().split(delimiter)
    for acc in line[2:]:
        accessions.append(acc.strip())
    print "Found", len(accessions), "arrays/strains."

    if filter_accessions:
        indices_to_include = []
        new_accessions = []
        for i, acc in enumerate(accessions):
            if acc in filter_accessions:
                indices_to_include.append(i)
                new_accessions.append(acc)
        accessions = new_accessions
        print "Loading only", len(accessions), "arrays/strains."

    positions = []
    snps = []
    line = f.next().split(delimiter)
    assert len(line)==len(accessions)+2
    old_chromosome = int(line[0])
    positions.append(int(line[1]))
    snps.append(np.array(line[2:], dtype=dtype))
    i = 0

    if filter < 1.0:
        if filter_accessions:
            for i, line in enumerate(f):
                if random.random() > filter: continue
                line_list = line.split(delimiter)
                chromosome = int(line_list[0])
                if chromosome != old_chromosome:  # Then save snps_data
                    snpsd = SNPsData(snps, positions, accessions=accessions, chromosome=old_chromosome)
                    snpsd_ls.append(snpsd)
                    chromosomes.append(old_chromosome)
                    positions = []
                    snps = []
                    old_chromosome = chromosome
                    sys.stderr.write("Loaded %s SNPs.\n" % i)
                positions.append(int(line_list[1]))
                l_list = line_list[2:]
                # snp = [int(l_list[j]) for j in indices_to_include]
                snp = np.array([l_list[j] for j in indices_to_include], dtype=dtype)
                snps.append(snp)
        else:
            for i, line in enumerate(f):
                if random.random() > filter: continue
                line_list = line.split(delimiter)
                chromosome = int(line_list[0])
                if chromosome != old_chromosome:  # Then save snps_data
                    snpsd = SNPsData(snps, positions, accessions=accessions, chromosome=old_chromosome)
                    snpsd_ls.append(snpsd)
                    chromosomes.append(old_chromosome)
                    positions = []
                    snps = []
                    old_chromosome = chromosome
                    sys.stderr.write("Loaded %s SNPs.\n" % i)
                positions.append(int(line_list[1]))
                # snps.append(map(int,line_list[2:]))
                snps.append(np.array(line_list[2:], dtype=dtype))
    else:
        if filter_accessions:
            for i, line in enumerate(f):
                line_list = line.split(delimiter)
                chromosome = int(line_list[0])
                if chromosome != old_chromosome:  # Then save snps_data
                    snpsd = SNPsData(snps, positions, accessions=accessions, chromosome=old_chromosome)
                    snpsd_ls.append(snpsd)
                    chromosomes.append(old_chromosome)
                    positions = []
                    snps = []
                    old_chromosome = chromosome
                    sys.stderr.write("Loaded %s SNPs.\n" % i)
                positions.append(int(line_list[1]))
                l_list = line_list[2:]
                # snp = [int(l_list[j]) for j in indices_to_include]
                snp = np.array([l_list[j] for j in indices_to_include], dtype=dtype)
                snps.append(snp)
        else:
            for i, line in enumerate(f):
                line_list = line.split(delimiter)
                chromosome = int(line_list[0])
                if chromosome != old_chromosome:  # Then save snps_data
                    snpsd = SNPsData(snps, positions, accessions=accessions, chromosome=old_chromosome)
                    snpsd_ls.append(snpsd)
                    chromosomes.append(old_chromosome)
                    positions = []
                    snps = []
                    old_chromosome = chromosome
                    sys.stderr.write("Loaded %s SNPs.\n" % i)
                positions.append(int(line_list[1]))
                # snps.append(map(int,line_list[2:]))
                snps.append(np.array(line_list[2:], dtype=dtype))

    f.close()
    snpsd = SNPsData(snps, positions, accessions=accessions, chromosome=old_chromosome)
    snpsd_ls.append(snpsd)
    chromosomes.append(old_chromosome)
    sys.stderr.write("Loaded %s SNPs.\n" % i)
    sd = SNPsDataSet(snpsd_ls, chromosomes, data_format=data_format)
    if use_pickle:
        print 'Saving a pickled version of genotypes.'
        f = open(pickle_file, 'wb')
        cPickle.dump(sd, f, protocol=-1)
        f.close()
    if filter_accessions_:
        sd.filter_accessions(filter_accessions_)
    return sd


def parse_snp_data(data_file, delimiter=",", missingVal='NA', data_format='binary', filter=1,
                   useDecoder=True, look_for_binary=True, filter_accessions=None,
                   use_pickle=True):
    """
    Load snps data..
    """
    if data_format == 'binary' and look_for_binary:
        print 'Looking for binary SNPs'
        if os.path.isfile(data_file) or os.path.isfile(data_file + '.pickled'):
            sd = parse_numerical_snp_data(data_file, delimiter=delimiter, missing_val=missingVal,
                        filter=filter, filter_accessions=filter_accessions,
                        use_pickle=use_pickle, dtype='int8', data_format=data_format)
        else:
            raise Exception('Genotype file not found: %s'%data_file)
    elif data_format in ['int', 'diploid_int']:
        sd = parse_numerical_snp_data(data_file, delimiter=delimiter, missing_val=missingVal,
                    filter=filter, filter_accessions=filter_accessions,
                    use_pickle=use_pickle, dtype='int8', data_format=data_format)
    elif data_format == 'float':
        sd = parse_numerical_snp_data(data_file, delimiter=delimiter, missing_val=missingVal,
                    filter=filter, filter_accessions=filter_accessions,
                    use_pickle=use_pickle, dtype='float32', data_format=data_format)
    elif data_format == 'nucleotides':
        print 'Looking for nucleotide SNPs'
        (snpsds, chromosomes) = parse_raw_snps_data(data_file, deliminator=delimiter, missing_val=missingVal,
                                debug_filter=filter, return_chromosomes=True)
        sd = SNPsDataSet(snpsds, chromosomes, data_format=data_format)
    else:
        print "Unknown file format!"
        raise Exception()
    return sd


def parse_chr_pos_list(datafile, delimiter=","):
    """
    Return a chr_pos list without loading all data...
    """

    sys.stderr.write("Loading file: %s ... \n" % datafile)

    chr_pos_list = []

    # Reading column data
    f = open(datafile, 'r')
    lines = f.readlines()
    f.close()

    line = lines[1].split(delimiter)
    withArrayIds = line[0] == "Chromosome"
    i = 1
    if withArrayIds:
        i += 1
    while i < len(lines):
        line = lines[i].split(delimiter)
        chr_pos_list.append((int(line[0]), int(line[1])))
        i += 1
    sys.stderr.write("Chromosomes and positions read\n")
    return chr_pos_list


def parse_plink_ped_file(file_prefix, only_binary_snps=True, debug=False):
    """
    Requires a .ped and a .map file.
    
    - Converts (on-the-fly) to a integer format. 
    - Ignores missing data?  Or imputes missing data.
    
    """
    assert only_binary_snps, 'Can only deal with binary SNPs for now.'
    map_filename = file_prefix + '.map'
    map_pickled_filename = map_filename + '.pickled'
    ped_filename = file_prefix + '.ped'
    ped_pickled_filename = ped_filename + '.pickled'
    if os.path.isfile(map_pickled_filename):
        print 'Loading pickled map file'
        chrom_pos_dict = cPickle.load(open(map_pickled_filename))
    else:
        chrom_pos_dict = {}
        with open(map_filename) as f:
            cur_chrom = -1
            for line in f:
                l = map(str.strip, line.split())
                chrom = int(l[0])
                if chrom != cur_chrom:
                    chrom_pos_dict[chrom] = {'positions':[]}
                    cur_chrom = chrom
                chrom_pos_dict[chrom]['positions'].append(int(l[3]))
        print 'The map file was loaded:'
        print 'Pickling..'
        cPickle.dump(chrom_pos_dict, open(map_pickled_filename, 'wb'), protocol=-1)
    num_markers = 0
    for chrom in sorted(chrom_pos_dict):
        n = len(chrom_pos_dict[chrom]['positions'])
        print 'Chromosome %d has %d markers.' % (chrom, n)
        num_markers += n
    print 'In total there are %d markers.' % num_markers

    nt_pair_map = {('0', '0'):0}
    for i, nt1 in enumerate(['A', 'C', 'G', 'T']):
        for j, nt2 in enumerate(['A', 'C', 'G', 'T']):
            nt_pair_map[(nt1, nt2)] = i * 4 + j + 1

    if os.path.isfile(ped_pickled_filename):
        print 'Loading pickled ped file'
        individ_dict = cPickle.load(open(ped_pickled_filename))
        print 'Pickled file was loaded.'
    else:
        individ_dict = {}
        with open(ped_filename) as f:
            for line_i, line in enumerate(f):
                if line_i % 10 == 0:
                    print line_i
                if debug and line_i > 500:
                    break
                l = map(str.strip, line.split('\t'))
                ind_id = int(l[1])
#                fam_id = int(l[0])
                assert not ind_id in individ_dict, 'The indivual %d is already in dictionary??' % ind_id
                individ_dict[ind_id] = {'fam_id':int(l[0]), 'pat_id':int(l[2]), 'mat_id':int(l[3]),
                            'sex':int(l[4]), 'phen_val':int(l[5])}
                # print individ_dict[ind_id]
                nt_pairs = sp.zeros(num_markers, dtype='int8')
                # missing_count = 0
                missing_indices = []
                for snp_i, genotype in enumerate(l[6:]):
                    nt_pair = genotype.split()
                    if nt_pair == ['0', '0']:
                        # missing_count += 1
                        missing_indices.append(snp_i)
                    nt_pairs[snp_i] = nt_pair_map[tuple(nt_pair)]
                # print '%d missing values were found' % missing_count
                individ_dict[ind_id]['snps'] = nt_pairs
                individ_dict[ind_id]['missing_indices'] = missing_indices
        print 'The ped file was loaded:'
        print 'Pickling..'
        cPickle.dump(individ_dict, open(ped_pickled_filename, 'wb'), protocol=-1)

    missing_indices_set = set()
    missing_nums = []
    num_retained = 0
    for ind_id in individ_dict:
        num_missing = len(individ_dict[ind_id]['missing_indices'])
        if num_missing < 400:
            num_retained += 1
            for i in individ_dict[ind_id]['missing_indices']:
                missing_indices_set.add(i)
        missing_nums.append(num_missing)
    print len(missing_indices_set)
    print num_retained

    snp_mat = sp.zeros((len(individ_dict), num_markers), dtype='int8')
    for i, ind_id in enumerate(individ_dict):
                snp_mat[i] = individ_dict[ind_id]['snps']

    print 'Finished construcing matrix.'

    num_weird_snps = 0
    snps = []
    for snp_i in range(num_markers):
        if snp_i % 1000 == 0:
            print snp_i, num_weird_snps
            snp = snp_mat[:, snp_i]
        alleles = sp.unique(snp)
        if 0 in alleles:
            if 3 <= len(alleles) <= 4:
                snps.append(snp)
            else:
                num_weird_snps += 1
        else:
            if 2 <= len(alleles) <= 3:
                snps.append(snp)
            else:
                num_weird_snps += 1

    print 'Number of weird SNPs is %d, out of %d' % (num_weird_snps, num_markers)
#    accessions = individ_dict.keys()


def parse_plink_tped_file(file_prefix, imputation_type='simple', return_kinship=False):
    """
    Requires a .tped file in 12 format.
    
    - Converts (on-the-fly) to a integer format. 
    - Imputes missing data.
    """
    tped_filename = file_prefix + '.tped'
    tped_pickled_filename = tped_filename + '.imputed.pickled'
    tfam_filename = file_prefix + '.tfam'
    tfam_pickled_filename = tfam_filename + '.pickled'

    if os.path.isfile(tfam_pickled_filename):
        print 'Loading pickled tfam file'
        individs, sex_list = cPickle.load(open(tfam_pickled_filename))
        print 'Pickled tfam file was loaded.'
    else:
        individs = []
        sex_list = []
        with open(tfam_filename) as f:
            for line in f:
                l = map(str.strip, line.split())
                individs.append(l[1])
                sex_list.append(int(l[4]))
        cPickle.dump((individs, sex_list), open(tfam_pickled_filename, 'wb'), protocol=-1)
    num_individs = len(individs)

#    k_mat = sp.zeros((num_individs, num_individs))
    if os.path.isfile(tped_pickled_filename):
        print 'Loading pickled tped file'
        chrom_pos_snp_dict = cPickle.load(open(tped_pickled_filename))
        print 'Pickled tped file was loaded.'
    else:
        chrom_pos_snp_dict = {}
        with open(tped_filename) as f:
            cur_chrom = -1
            for line_i, line in enumerate(f):
                if line_i % 1000 == 0:
                    print line_i
                l = map(str.strip, line.split())
                chrom = int(l[0])
                if chrom != cur_chrom:
                    chrom_pos_snp_dict[chrom] = {'positions':[], 'snps':[]}
                    cur_chrom = chrom
                chrom_pos_snp_dict[chrom]['positions'].append(int(l[3]))
                snp = sp.zeros(num_individs, dtype='int8')
                j = 0
                w_missing = False
                for i in range(4, 2 * num_individs + 4, 2):
                    nt1 = int(l[i])
                    nt2 = int(l[i + 1])
                    if nt1 == 0  or nt2 == 0:
                        snp[j] = 3
                        w_missing = True
                    elif nt1 == 2 and nt2 == 2:
                        snp[j] = 2
                    elif nt1 != 1  or nt2 != 1:
                        snp[j] = 1
#                    #Calculating K
#                    for ind_i in range(j):
#                        if snp[j] != 3 and snp[ind_i] != 3:
#                            k_mat[ind_i, j] = int(snp[j] == snp[ind_i]) + 0.5 * int(sp.absolute(snp[j] - snp[ind_i]) == 1)
#                            k_mat[ind_i, j] += 1
                    j += 1
#                print k_mat

                bin_counts = sp.bincount(snp)
                if w_missing:

                    if imputation_type == 'simple':
                        mean = (bin_counts[1] + 2 * bin_counts[2]) / (bin_counts[0] + bin_counts[1] + bin_counts[2])
                        snp[snp == 3] = round(mean)
                    if imputation_type == 'simple2':
                        snp[snp == 3] = sp.argmax(bin_counts[:-1])

                chrom_pos_snp_dict[chrom]['snps'].append(snp)
        cPickle.dump(chrom_pos_snp_dict, open(tped_pickled_filename, 'wb'), protocol=-1)

    chromosomes = sorted(chrom_pos_snp_dict.keys())
    snpsds = []
    for chrom in chromosomes:
        snps = chrom_pos_snp_dict[chrom]['snps']
        positions = chrom_pos_snp_dict[chrom]['positions']
        snpsds.append(SNPsData(snps, positions, accessions=individs, chromosome=chrom))
    sd = SNPsDataSet(snpsds, chromosomes, data_format='diploid_int')
    print 'SNPsDataSet constructed!'

    if return_kinship:
        print 'Loading the kinship matrix'
        ibs_filename = file_prefix + '.mibs'
        ibs_pickled_filename = ibs_filename + '.pickled'
        if os.path.isfile(ibs_pickled_filename):
            print 'Loading pickled IBS kinship file'
            l = cPickle.load(open(ibs_pickled_filename))
            K = l[0]
            print 'Pickled IBS kinship was loaded.'
        else:
            print 'Loading K...'
            K = sp.zeros((num_individs, num_individs), dtype='double')
            with open(ibs_filename) as f:
                for i, line in enumerate(f):
                    K[i] = map(float, line.split())
            cPickle.dump([K, individs], open(ibs_pickled_filename, 'wb'), protocol=-1)
            print 'K was loaded.'
        return sd, K
    return sd

# def load_snps_call_method(call_method_id=75, data_format='binary', debug_filter=1.0, min_mac=5):
#    file_prefix = '%s%d/' % (env['cm_dir'], call_method_id)
#    data_file = file_prefix + 'all_chromosomes_%s.csv' % data_format
#    if os.path.isfile(data_file):
#        print 'Found data in one file: %s' % data_file
#        return parse_snp_data(data_file , format=data_format, filter=debug_filter)
#    else:
#        data_file = file_prefix + 'chr_1_%s_mac%d.csv' % (data_format, min_mac)
#        data_file_mac0 = file_prefix + 'chr_1_%s_mac%d.csv' % (data_format, 0)
#        if os.path.isfile(data_file) or os.path.isfile(data_file_mac0):
#            print 'Found data in 5 files, in data format: %s' % data_format
#            return load_full_sequence_data(file_prefix, data_format=data_format,
#                        min_mac=min_mac, debug_filter=debug_filter)
#        else:
#            print 'Files in %s format were not found.' % data_format
#
#
#    print 'Looking for raw nucleotide files..'
#    nt_data_file = file_prefix + 'all_chromosomes_nucleotides.csv'
#    if os.path.isfile(nt_data_file):
#        print 'Found data file, loading and then attempting to convert to %s format' % data_format
#        (snpsds, chromosomes) = parse_raw_snps_data(nt_data_file, deliminator=',', missing_val='',
#                            debug_filter=debug_filter, return_chromosomes=True)
#        sd = SNPsDataSet(snpsds, chromosomes, data_format='nucleotides')
#        sd.convert_data_format(data_format)
#        data_file = file_prefix + 'all_chromosomes_%s.csv' % data_format
#        sd.writeToFile(data_file)
#        return load_snps_call_method(call_method_id=call_method_id, data_format=data_format,
#                    debug_filter=debug_filter, min_mac=min_mac)
#
#    else:
#        nt_data_file = file_prefix + 'chr_1_nucleotides_mac0.csv'
#        if os.path.isfile(nt_data_file):
#            print 'Found data in one file, now attempting to convert to %s format' % data_format
#            for chrom in [1, 2, 3, 4, 5]:
#                nt_data_file = file_prefix + 'chr_%d_nucleotides_mac0.csv' % chrom
#                (snpsds, chromosomes) = parse_raw_snps_data(nt_data_file, deliminator=',',
#                                    missing_val='N', debug_filter=debug_filter,
#                                    return_chromosomes=True)
#                sd = SNPsDataSet(snpsds, chromosomes, data_format='nucleotides')
#                sd.convert_data_format(data_format)
#                data_file = file_prefix + 'chr_%d_%s_mac0.csv' % (chrom, data_format)
#                sd.writeToFile(data_file)
#
#            return load_snps_call_method(call_method_id=call_method_id, data_format=data_format,
#                        debug_filter=debug_filter, min_mac=min_mac)
#
#
#    print 'There were problems with loading the genotype data..'
#    print 'Genotype files were not found in the given directory: %s' % file_prefix
#    raise NotImplementedError
#
#
# def load_hdf5_snps_call_method(call_method_id=75, data_format='binary', debug_filter=1.0):
#    file_prefix = '%s%d/' % (env['cm_dir'], call_method_id)
#    data_file = file_prefix + 'all_chromosomes_%s.hdf5' % data_format
#    if os.path.isfile(data_file):
#        sd = snps_data_set(data_file)
#    else:
#        #If not found, then load old data and create the HDF5 file.
#        sd = load_snps_call_method(call_method_id=call_method_id, data_format=data_format, debug_filter=1.0, min_mac=0)
#        sd = snps_data_set(data_file, sd=sd)
#    sd.filter_random_snps(debug_filter)
#    return sd


def load_full_sequence_data(file_prefix, data_format='diploid_int', min_mac=5, chromosomes=[1, 2, 3, 4, 5],
                debug_filter=1.0):
    print "Loading sequence data in data format: %s, and with min MAC: %d" % (data_format, min_mac)
    file_name = file_prefix + 'chr_%d_%s_mac%d.csv' % (1, data_format, min_mac)
    pickled_file_name = file_name + '.pickled'
    print 'Looking for file: %s' % file_name
    if min_mac > 0 and not (os.path.isfile(file_name) or os.path.isfile(pickled_file_name)):
        file_name = file_prefix + 'chr_%d_%s_mac%d.csv' % (1, data_format, 0)
        print 'File not found, instead trying: %s' % file_name
        if os.path.isfile(file_name):
            file_mac = 0
        else:
            raise Exception('Data file not found')
    else:
        file_mac = min_mac

    t = time.time()
    snpsds = []
    num_snps = 0
    for chrom in chromosomes:
        file_name = file_prefix + 'chr_%d_%s_mac%d.csv' % (chrom, data_format, file_mac)
        pickled_file_name = file_name + '.pickled'
        print 'Attempting to load pickled file: %s ' % pickled_file_name
        if os.path.isfile(pickled_file_name):
            sd = cPickle.load(open(pickled_file_name))
        else:
            sd = parse_snp_data(file_name, format=data_format, filter=debug_filter, use_pickle=True)
        if min_mac != file_mac:
            sd.filter_mac_snps(min_mac)
            file_name = file_prefix + 'chr_%d_%s_mac%d.csv' % (chrom, data_format, min_mac)
            pickled_file_name = file_name + '.pickled'
            cPickle.dump(sd, open(pickled_file_name, 'wb'), protocol=-1)
        print "Done."

        if debug_filter < 1.0:
            sd.sample_snps(debug_filter)
        num_snps += len(sd.snpsDataList[0].snps)

        snpsds.append(sd.snpsDataList[0])
    t = time.time() - t
    print 'Loaded % d SNPs.' % num_snps
    print 'It took % d minutes and % 0.2f seconds to load the SNPs' % (t / 60, t % 60)
    sd = SNPsDataSet(snpsds, chromosomes, data_format=data_format)
    print 'Loaded % d SNPs in total.' % sd.num_snps()
    return sd


def parse_tair_gff_file(chrom=None, reg_start_pos=None, reg_end_pos=None, only_genes=False,
                        gff_file='at_data/TAIR10_GFF3_genes_transposons.gff.tsv',
                        func_desc_file='at_data/TAIR10_functional_descriptions.tsv'):
    """
    Loads the TAIR GFF file.
    """
    pickled_filename = gff_file + '.pickled'
    if os.path.isfile(pickled_filename) and chrom == None:
        gene_dict = cPickle.load(open(pickled_filename))
    else:
        print 'Plowing through the tsv files.'
        with open(gff_file) as f:
            lines = f.readlines()
        gene_dict = {}
        for line in lines:
            l = line.split()
            annotation_type = l[2]
            info_list = l[8].split(';')
            chromosome = l[0][3]
            if chrom != None and str(chrom) != chromosome:
                continue
            start_pos = int(l[3])
            end_pos = int(l[4])
            strand = l[6]
            frame = l[7]
            if annotation_type in ['gene', 'transposable_element_gene', 'pseudogene']:
                tair_id = info_list[0][3:].strip()
                note = info_list[1][5:].strip()
                name = info_list[2][5:].strip()
                assert tair_id not in gene_dict, 'Gene is already in the dictionary'
                gene_dict[tair_id] = {'gene_type':annotation_type, 'name':name, 'note':note,
                            'strand':strand, 'start_pos':start_pos, 'end_pos':end_pos,
                            'chromosome':chromosome}

            elif annotation_type in ['mRNA', 'tRNA', 'rRNA', 'miRNA', 'ncRNA', 'snoRNA', 'snRNA',
                        'pseudogenic_transcript']:
                tair_isoform_id = info_list[0][3:].strip()
                tair_id = tair_isoform_id.split('.')[0]
                name = info_list[2][5:].strip()
                assert tair_id  in gene_dict, "Gene isn't in the  dictionary"
                assert tair_isoform_id not in gene_dict[tair_id], "%s already found?" % annotation_type
                gene_dict[tair_id][tair_isoform_id] = {'RNA_type':annotation_type, 'name':name,
                                    'strand':strand, 'start_pos':start_pos,
                                    'end_pos':end_pos, 'exons':[], 'cds':[]}

            elif annotation_type == 'protein':
                tair_isoform_id = info_list[1][5:].strip()
                tair_id = tair_isoform_id.split('.')[0]
                name = info_list[1][5:].strip()
                gene_dict[tair_id][tair_isoform_id]['protein'] = {'name':name, 'start_pos':start_pos,
                                        'strand':strand, 'end_pos':end_pos}

            elif annotation_type in ['five_prime_UTR', 'three_prime_UTR']:
                tair_isoform_id = info_list[0][7:].strip()
                tair_id = tair_isoform_id.split('.')[0]
                gene_dict[tair_id][tair_isoform_id][annotation_type] = {'start_pos':start_pos,
                                            'end_pos':end_pos}

            elif annotation_type == 'CDS':
                tair_isoform_id = info_list[0].split(',')[0][7:].strip()
                tair_id = tair_isoform_id.split('.')[0]
                gene_dict[tair_id][tair_isoform_id]['cds'].append({'start_pos':start_pos,
                                            'end_pos':end_pos,
                                            'frame':int(frame)})

            elif annotation_type in ['exon', 'pseudogenic_exon']:
                tair_isoform_id = info_list[0][7:].strip()
                tair_id = tair_isoform_id.split('.')[0]
                d = {'start_pos':start_pos, 'end_pos':end_pos}
                gene_dict[tair_id][tair_isoform_id]['exons'].append(d)

            elif annotation_type == 'transposable_element':
                tair_id = info_list[0][3:].strip()
                name = info_list[1][5:].strip()
                alias = info_list[2][6:].strip()
                gene_dict[tair_id] = {'gene_type':annotation_type, 'name':name, 'strand':strand,
                            'start_pos':start_pos, 'end_pos':end_pos, 'alias':alias,
                            'chromosome':chromosome, 'fragments':[]}

            elif annotation_type == 'transposon_fragment':
                tair_id = info_list[0][7:].strip()
                assert tair_id in gene_dict, 'We have a fragment but where is the transposon?'
                gene_dict[tair_id]['fragments'].append({'start_pos':start_pos, 'end_pos':end_pos})

        print 'Done loading the GFF file'
        print 'Now adding functional description'
        with open(func_desc_file) as f:
            f.next()
            for line in f:
                l = map(str.strip, line.split('\t'))
                tair_isoform_id = l[0]
                tair_id = tair_isoform_id.split('.')[0]
                gene_type = l[1]
                short_description = l[2]
                curator_summary = l[3]
                computational_description = l[4]
                if tair_id in gene_dict:
                    gene_dict[tair_id][tair_isoform_id]['functional_description'] = \
                        {'type':gene_type, 'short_description':short_description,
                        'curator_summary':curator_summary,
                        'computational_description':computational_description}
        if chrom is None:
            cPickle.dump(gene_dict, open(pickled_filename, 'wb'), protocol=-1)
        elif reg_start_pos != None and reg_end_pos != None:
            tair_ids = gene_dict.keys()
            for tair_id in tair_ids:
                if not (gene_dict[tair_id]['start_pos'] <= reg_end_pos and \
                        gene_dict[tair_id]['end_pos'] >= reg_start_pos):
                    del gene_dict[tair_id]

    if only_genes:
        tids = gene_dict.keys()
        for tid in tids:
            if gene_dict[tid]['gene_type'] != 'gene':
                del gene_dict[tid]

    return gene_dict    


def _parse_map_file_():
    """
    Returns a dictionary for the NFBC data map file.
    """
    chrom_d = {}
    for chrom in range(1, 24):
        chrom_d[chrom] = {'positions':[], 'var_ids':[]}
    var_d = {}
    with open('/Users/bjarni.vilhjalmsson/Projects/Data/NFBC_results/NFBC_20091001.map') as f:
        for line in f:
            l = line.split()
            chrom = int(l[0])
            var_id = l[1]
            position = int(l[3])
            chrom_d[chrom]['positions'].append(position)
            chrom_d[chrom]['var_ids'].append(var_id)
            var_d[var_id] = {'chrom':chrom, 'position':position}
    return chrom_d, var_d


if __name__ == "__main__":
    pass
