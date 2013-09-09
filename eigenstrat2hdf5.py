"""
A basic parser for eigenstrat formated files to a more convenient HDF5 format.
"""

def load_eigenstrat_genotypes(in_file_prefix='eigenstrat_file_prefix',
                              out_file_prefix='hdf5_file_prefix',
                              impute_type='mode',
                              filter_monomorphic_snps=True,
                              missing_val_thr=0.1):
    """
    Parses eigenstrat formated genotype files to a HDF5 format.  It requires the h5py and scipy package.
    
    Ideally the genotypes are imputed apriory, otherwise a rough imputation 
    (the most common genotype) is used for missing genotypes.
    
    Notes: 
        Assumes the files are in diploid format!
    
    """
    import h5py
    import scipy as sp
    import os
    import sys
    
    data_file_prefix = '%s_mv%0.2f_imp_%s.' % (out_file_prefix, missing_val_thr, impute_type)    
     
    genotype_data = {}
    
    # Setting the HDF5 file up
    h5py_file_name = data_file_prefix + 'h5py'
    if os.path.isfile(h5py_file_name):
        print 'Overwriting: %s' % h5py_file_name
        os.remove(h5py_file_name)
    h5py_file = h5py.File(h5py_file_name)
    genotype_data['h5py_file'] = h5py_file_name
        
    
    # Fill out individuals data, if available
    i_filename = '%sind' % (in_file_prefix)
    if os.path.isfile(i_filename):
        iids = []
        phens = []
        genders = []
        with open(i_filename) as f:
            for line in f:
                l = (line.strip()).split()
                iids.append(l[0])
                genders.append(l[1])
                phens.append(l[2])
        ind_group = h5py_file.create_group('indivs')
        ind_group.create_dataset('indiv_ids', data=iids)
        ind_group.create_dataset('sex', data=genders)
        ind_group.create_dataset('phenotype', data=phens)
    else:
        print 'Individual information file not found: %s' % i_filename
        
    tot_num_snps = 0
    tot_num_duplicated_snps_removed = 0
    tot_num_missing_val_snps_removed = 0
    tot_num_monomorphic_snps_removed = 0
    
    
    # Open the genotype files.
    s_filename = '%ssnp' % (in_file_prefix) 
    g_filename = '%sgeno' % (in_file_prefix)
    print 'Starting to parse files:\n\t %s \n\t %s' % (s_filename, g_filename)
    sf = open(s_filename) 
    gf = open(g_filename) 
    

    # Figure out sample size, number of SNPs, etc. 
    # Initialize HDF5 file.

    # Setting up containers.
    curr_chrom = 1
    curr_hdf5_group = h5py_file.create_group('chrom_%d' % curr_chrom)
    snps_mat = []
    positions = []
    sids = []
    nts_list = []
    nt_counts_list = []
    missing_counts = []
    freqs = []
    num_missing_removed = 0
    num_monomorphic_removed = 0
    num_duplicated_snps_removed = 0

    print 'Starting to parse SNP files'
    for s_line in sf:
        g_line = gf.next()
        sl = s_line.split()
        pos = int(sl[3])
        chrom = int(sl[1])
        sid = sl[0]

        if chrom != curr_chrom:
            # Report statistics and store stuff
            print 'Finished with Chromosome %d' % curr_chrom
            print 'Number of SNPs removed due to too many missing values: %d' % num_missing_removed
            print 'Number of duplicated SNPs removed: %d' % num_duplicated_snps_removed
            print 'Number of monomorphic SNPs removed: %d' % num_monomorphic_removed
            print 'Number of SNPs retained: %d' % len(positions)
            snps = sp.array(snps_mat, dtype='int8')
            curr_hdf5_group.create_dataset('raw_snps', compression='lzf', data=snps)
            h5py_file.flush()
            print 'Raw SNPs stored'
            snps = snps.T
            snps = (snps - sp.mean(snps, 0)) / sp.std(snps, 0)
            curr_hdf5_group.create_dataset('snps', compression='lzf', data=snps.T)
            h5py_file.flush()
            print 'Normalized SNPs stored'
            del snps
            del snps_mat
            curr_hdf5_group.create_dataset('positions', compression='lzf', data=positions)
            curr_hdf5_group.create_dataset('nts', compression='lzf', data=nts_list)
            curr_hdf5_group.create_dataset('nt_counts', compression='lzf', data=sp.array(nt_counts_list))
            curr_hdf5_group.create_dataset('missing_counts', compression='lzf', data=missing_counts)
            curr_hdf5_group.create_dataset('freqs', compression='lzf', data=freqs)
            curr_hdf5_group.create_dataset('snp_ids', compression='lzf', data=sids)        
            h5py_file.flush()
            sys.stdout.flush()

            # Reset containers
            curr_chrom = chrom
            curr_hdf5_group = h5py_file.create_group('chrom_%d' % curr_chrom)
            snps_mat = []
            positions = []
            sids = []
            nts_list = []
            nt_counts_list = []
            missing_counts = []
            freqs = []
            num_missing_removed = 0
            num_monomorphic_removed = 0
            num_duplicated_snps_removed = 0
            
        
        # Debug filter
                    
        nt = (sl[4], sl[5])

        snp = sp.array(map(int, g_line.strip()), dtype='int8')
        num_indiv = len(snp)
        bin_counts = sp.bincount(snp)
#        print bin_counts
        missing_count = bin_counts[-1]

        # Filtering SNPs with too many missing values
        if missing_count > missing_val_thr * 2 * num_indiv:
            num_missing_removed += 1
            tot_num_missing_val_snps_removed += 1
            continue

        nt_counts = list(bin_counts[:3])
        # Imputing the SNPs roughly by replacing missing values with the mode value.
        if impute_type == 'mode':
            v = sp.argmax(nt_counts)
            snp[snp == 9] = v
        else:
            raise Exception('Imputation type is unknown')

        bin_counts = sp.bincount(snp)
        nt_counts = list(bin_counts[:3])
        # Removing monomorphic SNPs
        if max(nt_counts) == sum(nt_counts):
            num_monomorphic_removed += 1
            tot_num_monomorphic_snps_removed += 1
            continue
        if len(nt_counts) == 2:
            nt_counts.append(0)
            
#        assert len(nt_counts) == 3, 'ARrrg'                    

        # Is this position already there?
        if len(positions) > 0 and pos == positions[-1]:
            num_duplicated_snps_removed += 1
            tot_num_duplicated_snps_removed += 1
            continue
        
        freq = sp.mean(snp) / 2.0            
        snps_mat.append(snp)
        positions.append(pos)
        sids.append(sid)
        nts_list.append(nt)
        nt_counts_list.append(nt_counts)
        missing_counts.append(missing_count)
        freqs.append(freq)

        tot_num_snps += 1
        


    # Report statistics and store stuff
    print 'Number of SNPs removed due to too many missing values: %d' % num_missing_removed
    print 'Number of duplicated SNPs removed: %d' % num_duplicated_snps_removed
    print 'Number of monomorphic SNPs removed: %d' % num_monomorphic_removed
    print 'Number of SNPs retained: %d' % len(positions)
    snps = sp.array(snps_mat, dtype='int8')
    curr_hdf5_group.create_dataset('raw_snps', compression='lzf', data=snps)
    h5py_file.flush()
    print 'Raw SNPs stored'
    snps = snps.T
    snps = (snps - sp.mean(snps, 0)) / sp.std(snps, 0)
    curr_hdf5_group.create_dataset('snps', compression='lzf', data=snps.T)
    h5py_file.flush()
    print 'Normalized SNPs stored'
    del snps
    del snps_mat
    curr_hdf5_group.create_dataset('positions', compression='lzf', data=positions)
    curr_hdf5_group.create_dataset('nts', compression='lzf', data=nts_list)
    curr_hdf5_group.create_dataset('nt_counts', compression='lzf', data=sp.array(nt_counts_list))
    curr_hdf5_group.create_dataset('missing_counts', compression='lzf', data=missing_counts)
    curr_hdf5_group.create_dataset('freqs', compression='lzf', data=freqs)
    curr_hdf5_group.create_dataset('snp_ids', compression='lzf', data=sids)        
    
                
    gf.close()
    sf.close()
    
    print 'Genotypes for %d individuals were parsed.' % num_indiv
    print 'Total number of SNPs parsed successfully was: %d' % tot_num_snps
    print 'Total number of SNPs removed due to too many missing values: %d' % tot_num_missing_val_snps_removed
    print 'Total number of SNPs removed due to monomorphicity: %d' % tot_num_monomorphic_snps_removed
    print 'Total number of duplicated SNPs removed: %d' % tot_num_duplicated_snps_removed
    h5py_file.close()
    sys.stdout.flush()
    
    print 'Done parsing genotypes.'
  
  
