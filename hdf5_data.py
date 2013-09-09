"""
Container for functions that work with HDF5 genotype/phenotype datasets.
(Usually useful for analysing human data.)
"""

import h5py
import scipy as sp
import sys

def calculate_ibd_kinship(hdf5_filename='/home/bv25/data/Ls154/Ls154_12.hdf5',
                          chunk_size=1000, overwrite=False):
    """
    Calculates a kinship matrix and stores it in the HDF5 file.
    """
    
    h5f = h5py.File(hdf5_filename)
    c22g = h5f['chrom_22']
    n_indivs = len(c22g['raw_snps'][0])
    if overwrite or not 'kinship' in h5f.keys():
        if 'kinship' in h5f.keys():
            print 'Overwriting kinship.'
        print 'Calculating kinship.'
        k_mat = sp.zeros((n_indivs, n_indivs), dtype='single')
        
        n_snps = 0
        for chrom in range(1, 6):
            chrom_str = 'chrom_%d' % chrom
            print 'Working on chromosome %d' % chrom
            cg = h5f[chrom_str]
            num_snps = len(cg['snps'])
            assert len(cg['snps']) == len(cg['raw_snps']), 'WTF?'
            for chunk_i, i in enumerate(range(0, num_snps, chunk_size)):
    #            if chunk_i % 2 != 0:
    #                continue
                end_i = min(i + chunk_size, num_snps)
                x = cg['snps'][i:end_i]
                n_snps += len(x)
                k_mat += sp.dot(x.T, x)
                del x
                sys.stdout.write('\b\b\b\b\b\b\b%0.2f%%' % (100.0 * (min(1, ((chunk_i + 1.0) * chunk_size) / num_snps))))
                sys.stdout.flush()
            sys.stdout.write('\b\b\b\b\b\b\b100.00%\n')
        k_mat = k_mat / float(n_snps)
        h5f.create_dataset('kinship', data=k_mat)
    else:
        print 'kinship already there.'
    
