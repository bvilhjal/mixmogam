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
    n_indivs = len(h5f['indiv_data']['indiv_ids'][...])
    if overwrite or not 'kinship' in h5f.keys():
        if 'kinship' in h5f.keys():
            print 'Overwriting kinship.'
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
    
