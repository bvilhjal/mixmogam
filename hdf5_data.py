"""
Container for functions that work with HDF5 genotype/phenotype datasets.
(Usually useful for analysing human data.)
"""

import h5py
import scipy as sp
import sys
import linear_models as lm
import time


def calculate_ibd_kinship(hdf5_filename='/home/bv25/data/Ls154/Ls154_12.hdf5',
                          chunk_size=1000, overwrite=False):
    """
    Calculates a kinship matrix and stores it in the HDF5 file.
    """
    
    h5f = h5py.File(hdf5_filename)
    n_indivs = len(h5f['indiv_data']['indiv_ids'][...])
    if overwrite or not 'kinship' in h5f.keys():
        if 'kinship' in h5f.keys():
            print 'Overwriting kinship.'
            del h5f['kinship']
        print 'Calculating kinship.'
        k_mat = sp.zeros((n_indivs, n_indivs), dtype='single')
        
        gg = h5f['genot_data']
        chromosomes = gg.keys()
        n_snps = 0
        for chrom in chromosomes:
            print 'Working on Chromosome %s' % chrom
            cg = gg[chrom]
            num_snps = len(cg['raw_snps'])
            normalized = 'snps' in cg.keys()                
            
            for chunk_i, i in enumerate(range(0, num_snps, chunk_size)):
    #            if chunk_i % 2 != 0:
    #                continue
                end_i = min(i + chunk_size, num_snps)
                if normalized:
                    x = cg['snps'][i:end_i]
                else:
                    x = cg['raw_snps'][i:end_i]
                    x = x.T
                    x = (x - sp.mean(x, 0)) / sp.std(x, 0)
                    x = x.T           
                n_snps += len(x)
                k_mat += sp.dot(x.T, x)
                del x
                sys.stdout.write('\b\b\b\b\b\b\b%0.2f%%' % (100.0 * (min(1, ((chunk_i + 1.0) * chunk_size) / num_snps))))
                sys.stdout.flush()
            sys.stdout.write('\b\b\b\b\b\b\b100.00%\n')
        k_mat = k_mat / float(n_snps)
        c = sp.sum((sp.eye(len(k_mat)) - (1.0 / len(k_mat)) * sp.ones(k_mat.shape)) * sp.array(k_mat))
        scalar = (len(k_mat) - 1) / c
        print 'Kinship scaled by: %0.4f' % scalar
        k = scalar * k_mat

        h5f.create_dataset('kinship', data=k)
    else:
        print 'kinship already there.'
    


def run_emmax(hdf5_filename='/home/bv25/data/Ls154/Ls154_12.hdf5',
              out_file='/home/bv25/data/Ls154/Ls154_results.hdf5'):
    """
    Apply the EMMAX algorithm to hdf5 formated genotype/phenotype data 
    """
    
    ih5f = h5py.File(hdf5_filename)
    n_indivs = len(ih5f['indiv_data']['indiv_ids'][...])
    assert 'kinship' in ih5f.keys(), 'Kinship is missing.  Please calculate that first!'
    k = ih5f['kinship']
    gg = ih5f['genot_data']
    ig = ih5f['indiv_data']
    
    # Get the phenotypes
    phenotypes = ig['phenotypes'][...]
    num_snps = ih5f['num_snps'][...]

    # Initialize the mixed model
    lmm = lm.LinearMixedModel(phenotypes)
    lmm.add_random_effect(k)
    # Calculate pseudo-heritability, etc.
    print 'Calculating the eigenvalues of K'
    s0 = time.time()
    eig_L = lmm._get_eigen_L_()
    print 'Done.'
    print 'Took %0.2f seconds' % (time.time() - s0)
    print "Calculating the eigenvalues of S(K+I)S where S = I-X(X'X)^-1X'"
    s0 = time.time()
    eig_R = lmm._get_eigen_R_(X=lmm.X)
    print 'Done'
    print 'Took %0.2f seconds' % (time.time() - s0)

    print 'Getting variance estimates'
    s0 = time.time()
    res = lmm.get_estimates(eig_L, method='REML', eig_R=eig_R)  # Get the variance estimates..
    print 'Done.'
    print 'Took %0.2f seconds' % (time.time() - s0)
    print 'pseudo_heritability:', res['pseudo_heritability']

    # Initialize results file
    oh5f = h5py.File(out_file)
    
    # Store phenotype_data
    oh5f.create_dataset('pseudo_heritability', data=sp.array(res['pseudo_heritability']))
    oh5f.create_dataset('ve', data=sp.array(res['ve']))
    oh5f.create_dataset('vg', data=sp.array(res['vg']))
    oh5f.create_dataset('max_ll', data=sp.array(res['max_ll']))
    
    # Construct results data containers
    chrom_res_group = oh5f.create_group('chrom_results')
    
    for chrom in gg.keys():
        crg = chrom_res_group.create_group(chrom)
        # Get the SNPs
        print 'Working on Chromosome: %s' % chrom
        snps = gg[chrom]['raw_snps'][...]
        
        # Now run EMMAX
        
        print "Running EMMAX"
        s1 = time.time()
        r = lmm._emmax_f_test_(snps, res['H_sqrt_inv'], with_betas=False, emma_num=0, eig_L=eig_L)
        secs = time.time() - s1
        if secs > 60:
            mins = int(secs) / 60
            secs = secs % 60
            print 'Took %d mins and %0.1f seconds.' % (mins, secs)
        else:
            print 'Took %0.1f seconds.' % (secs)
        crg.create_dataset('ps', data=r['ps'])
        crg.create_dataset('positions', data=gg[chrom]['positions'])
        crg.flush()

    ih5f.close()   
    oh5f.close()
        
        


def qq_plot(hdf5_results_file='/home/bv25/data/Ls154/Ls154_results.hdf5'):
    h5f = h5py.File(hdf5_results_file)
    for chrom in gg.keys():
    



def manhattan_plot(hdf5_results_file='/home/bv25/data/Ls154/Ls154_results.hdf5'):
        
    
