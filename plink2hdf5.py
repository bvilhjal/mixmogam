"""
A basic parser for tped plink formated files to a more convenient HDF5 format.
"""

def convert_plink_files(in_file_prefix='plink_file_prefix',
                         out_file_prefix='hdf5_file_prefix',
                         impute_type='mode',
                         filter_monomorphic_snps=True,
                         missing_val_thr=0.1):
    """
    Parses plink tped format to a HDF5 format.  It requires the h5py and scipy package.
    
    Ideally the genotypes are imputed apriory, otherwise a rough imputation 
    (the most common genotype) is used for missing genotypes.

    Notes: 
        Assumes the files are in diploid format!
    
    """
    import time
    import h5py
    import scipy as sp
    
    
    print 'Starting to parse genotypes'
    genotype_data = {}
    h5py_file = h5py.File(out_file_prefix + '.hdf5')
    genotype_data['hdf5p_file'] = h5py_file
            
            
    tot_num_snps = 0
    tot_num_missing_val_snps_removed = 0
    tot_num_ambiguous_loc_removed = 0
    curr_chrom = 1
    print 'Working on chromosome %d' % curr_chrom
    
    g_filename = '%s.tped' % (in_file_prefix) 
    s_filename = '%s.bim' % (in_file_prefix)
    i_filename = '%s.tfam' % (in_file_prefix)  

        
    
    indiv_ids = []
    phenotypes = [] 
    sex = []
    print 'Parsing individuals file: %s' % i_filename
    with open(i_filename) as f:
        for line in f:
            l = line.split()
            iid = l[0]
            indiv_ids.append(iid)
            sex.append(int(l[4]))
            phenotypes.append(float(l[5]))
    tot_num_indiv = len(indiv_ids)
            
    num_indiv = len(indiv_ids)
    print 'Found %d Individuals' % (num_indiv)

    print 'Parsing nucleotide map'
    nt_map = {}
    chromsomoes = []
    curr_chrom = 0
    with open(s_filename) as f:
        for line in f:
            l = line.split()
            chrom = l[0]
            if chrom != curr_chrom:
                chromsomoes.append(chrom)
            nt_map[l[1]] = (l[4], l[5]) 
    assert len(chromsomoes) == len(set(chromsomoes)), 'Chromosomes need to be in order.'
    
        
    position = -1
    # Initializing containers.
    snps_mat = [] 
    positions = []
    sids = []
    nts_list = []
    nt_counts_list = []
    missing_counts = []
    freqs = []
    num_missing_removed = 0
    num_monomorphic_removed = 0
    num_ambiguous_loc_removed = 0
    t0 = time.time()

    print 'Starting to parse SNP files'
    gf = open(g_filename)
    for g_line in gf:
#        if random.random() > 0.01:
#            continue
        gl = g_line.split()
        chrom = int(gl[0])
        if chrom != curr_chrom:
            
            # Store everything and reset.
            print 'Number of SNPs removed due to too many missing values: %d' % num_missing_removed
            print 'Number of SNPs removed due to ambiguous location: %d' % num_ambiguous_loc_removed
            print 'Number of monomorphic SNPs removed: %d' % num_monomorphic_removed
            print 'Number of SNPs retained: %d' % len(positions)
            print 'Number of individuals: %d' % num_indiv
            snps = sp.array(snps_mat, dtype='int8')
            h5py_chrom_group = h5py_file.create_group('chrom_%d' % curr_chrom)
            h5py_chrom_group.create_dataset('indiv_ids', data=indiv_ids)
            h5py_chrom_group.create_dataset('sex', data=sex)
            h5py_chrom_group.create_dataset('phenotypes', data=phenotypes)
            h5py_chrom_group.create_dataset('raw_snps', compression='lzf', data=snps)
            h5py_chrom_group.create_dataset('positions', compression='lzf', data=positions)
            h5py_chrom_group.create_dataset('nts', compression='lzf', data=nts_list)
            h5py_chrom_group.create_dataset('nt_counts', compression='lzf', data=nt_counts_list)
            h5py_chrom_group.create_dataset('missing_counts', compression='lzf', data=missing_counts)
            h5py_chrom_group.create_dataset('freqs', compression='lzf', data=freqs)
            h5py_chrom_group.create_dataset('snp_ids', compression='lzf', data=sids)        
            tot_num_snps += len(positions)
            tot_num_missing_val_snps_removed += num_missing_removed
            tot_num_ambiguous_loc_removed += num_ambiguous_loc_removed
            h5py_file.flush()         
            t1 = time.time()
            t = t1 - t0
            print 'It took %d minutes and %0.2f seconds to parse Chromosome %d.' % (t / 60, t % 60, curr_chrom)
            t0 = time.time()

            

            # Reset containers
            snps_mat = [] 
            positions = []
            sids = []
            nts_list = []
            nt_counts_list = []
            missing_counts = []
            freqs = []
            num_missing_removed = 0
            num_monomorphic_removed = 0
            num_ambiguous_loc_removed = 0
               
            curr_chrom = chrom

        sid = gl[1]
        prev_position = position
        position = int(gl[3])

        # Skipping unmappable locations
        if position == prev_position:
            num_ambiguous_loc_removed += 1
            continue
        if position == 0:
            num_ambiguous_loc_removed += 1
            continue

        nt = nt_map[sid]
                
        snp0 = sp.array(map(int, (g_line.strip()).split()[4:]), 'int8')
        a = sp.arange(tot_num_indiv * 2)
        even_map = a % 2 == 0
        odd_map = a % 2 == 1
        snp = snp0[even_map] + snp0[odd_map] - 2
        snp[snp < 0] = 9
                   
        bin_counts = sp.bincount(snp)
        

        if len(bin_counts) > 3:
            missing_count = bin_counts[-1]
            # Filtering SNPs with too many missing values
            if missing_count > missing_val_thr * 2 * num_indiv:
                num_missing_removed += 1
                continue
            elif impute_type == 'mode':
                nt_counts = bin_counts[:3]                    
                v = sp.argmax(nt_counts)
                snp[snp == 9] = v
                bin_counts = sp.bincount(snp)
            else:
                raise Exception('Imputation type is unknown')
        else:
            missing_count = 0

        assert len(bin_counts) < 4, 'Issues with nucleotides.'
        nt_counts = bin_counts[:3]                    
        if len(nt_counts) == 2:
            nt_counts = sp.array([nt_counts[0], nt_counts[1], 0])
        elif len(nt_counts) == 1:
            nt_counts = sp.array([nt_counts[0], 0, 0])
            

        # Removing monomorphic SNPs
        if filter_monomorphic_snps:
            if max(nt_counts) == sum(nt_counts):
                num_monomorphic_removed += 1
                continue
        
        freq = sp.mean(snp) / 2.0            
        snps_mat.append(snp)
        positions.append(position)
        sids.append(sid)
        nts_list.append(nt)
        nt_counts_list.append(nt_counts)
        missing_counts.append(missing_count)
        freqs.append(freq) 

    # Store everything and reset.
    print 'Number of SNPs removed due to too many missing values: %d' % num_missing_removed
    print 'Number of SNPs removed due to ambiguous location: %d' % num_ambiguous_loc_removed
    print 'Number of monomorphic SNPs removed: %d' % num_monomorphic_removed
    print 'Number of SNPs retained: %d' % len(positions)
    print 'Number of individuals: %d' % num_indiv
    snps = sp.array(snps_mat, dtype='int8')
    h5py_chrom_group = h5py_file.create_group('chrom_%d' % chrom)
    h5py_chrom_group.create_dataset('indiv_ids', data=indiv_ids)
    h5py_chrom_group.create_dataset('sex', data=sex)
    h5py_chrom_group.create_dataset('phenotypes', data=phenotypes)
    h5py_chrom_group.create_dataset('raw_snps', compression='lzf', data=snps)
    h5py_chrom_group.create_dataset('positions', compression='lzf', data=positions)
    h5py_chrom_group.create_dataset('nts', compression='lzf', data=nts_list)
    h5py_chrom_group.create_dataset('nt_counts', compression='lzf', data=nt_counts_list)
    h5py_chrom_group.create_dataset('missing_counts', compression='lzf', data=missing_counts)
    h5py_chrom_group.create_dataset('freqs', compression='lzf', data=freqs)
    h5py_chrom_group.create_dataset('snp_ids', compression='lzf', data=sids)        
    tot_num_snps += len(positions)
    tot_num_missing_val_snps_removed += num_missing_removed
    tot_num_ambiguous_loc_removed += num_ambiguous_loc_removed
    h5py_file.flush()         
    t1 = time.time()
    t = t1 - t0
    print 'It took %d minutes and %0.2f seconds to parse chromosome %d.' % (t / 60, t % 60, chrom)

    
    gf.close()
    
    print 'Total number of SNPs parsed successfully was: %d' % tot_num_snps
    print 'Total number of SNPs removed due to too many missing values: %d' % tot_num_missing_val_snps_removed
    print 'Total number of SNPs removed due to ambiguous locations: %d' % tot_num_ambiguous_loc_removed
    h5py_file.close()
    
    print 'Done parsing genotypes.'
    
    
