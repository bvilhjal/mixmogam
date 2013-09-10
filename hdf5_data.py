"""
Container for functions that work with HDF5 genotype/phenotype datasets.
(Usually useful for analysing human data.)
"""

import h5py
import scipy as sp
import sys
import linear_models as lm
import analyze_gwas_results as agr
import time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


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
    oh5f.create_dataset('num_snps', data=ih5f['num_snps'])
    
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
        oh5f.flush()

    ih5f.close()   
    oh5f.close()
        
        


def qq_plot(hdf5_results_file='/home/bv25/data/Ls154/Ls154_results.hdf5',
            png_file_prefix='/home/bv25/data/Ls154/Ls154_results'):
    """
    Plot QQ-plot for a single HDF5 result
    """
    h5f = h5py.File(hdf5_results_file)
    chrom_res_group = h5f['chrom_results']
    pvals = sp.empty(h5f['num_snps'][...])
    i = 0
    for chrom in chrom_res_group.keys():
        crg = chrom_res_group[chrom]
        n = len(crg['ps'])
        pvals[i:i + n] = crg['ps'][...]
        i += n
    
    quantiles = agr.get_quantiles(pvals)
    log_quantiles = agr.get_log_quantiles(pvals, max_val=7)
    qq_plot_png_filename = png_file_prefix + '_qq.png'
    qq_log_plot_png_filename = png_file_prefix + '_qq_log.png'
    agr.simple_qqplot([quantiles], png_file=qq_plot_png_filename)
    agr.simple_log_qqplot([log_quantiles], png_file=qq_log_plot_png_filename, max_val=7)
    


def manhattan_plot(hdf5_results_file='/home/bv25/data/Ls154/Ls154_results.hdf5', png_file='/home/bv25/data/Ls154/Ls154_results_manhattan.png',
                   max_log_pval=None, filter_pval=0.02, ylab="$-$log$_{10}(p-$value$)$", plot_bonferroni=False,
                   b_threshold=None, markersize=3, chrom_col_map=None):
    """
    Plot a Manhattan plot for a single HDF5 result
    """

    chrom_res_dict = {}
    h5f = h5py.File(hdf5_results_file)
    chrom_res_group = h5f['chrom_results']
    num_snps = h5f['num_snps'][...]
    for chrom in chrom_res_group.keys():
        crg = chrom_res_group[chrom]
        ps = crg['ps'][...]
        positions = crg['positions'][...]
        ps_filter = ps < filter_pval
        chrom_end = positions[-1]
        chrom_res_dict[chrom] = {'log_ps':-sp.log10(ps[ps_filter]), 'positions': positions[ps_filter], 'chrom_end':chrom_end}
        
    h5f.close()

    chromosomes = chrom_res_dict.keys()
    chromosomes.sort()

    if not max_log_pval:
        max_log_pvals = []
        for chrom in chromosomes:
            max_log_pvals.append(chrom_res_dict[chrom]['log_ps'].max())
        max_log_pval = max(max_log_pvals)


    offset = 0
    tick_positions = []
    tick_strings = []
    plt.figure(figsize=(11, 2.8))
    plt.axes([0.045, 0.15, 0.95, 0.71])
    chr_offsets = []
    for chrom in chromosomes:
        chr_offsets.append(offset)
        log_ps = chrom_res_dict[chrom]['log_ps']
        positions = chrom_res_dict[chrom]['positions'] / 1000000.0
        plot_positions = offset + positions 

        pval_truc_filter = log_ps > max_log_pval
        if sum(pval_truc_filter) > 0:
            print '%d -log p-values were truncated at %0.f.' % max_log_pval
            log_ps[pval_truc_filter] = max_log_pval

        if not chrom_col_map:
            plt.plot(plot_positions, log_ps, ".", markersize=markersize, alpha=0.7, mew=0)
        else:
            color = chrom_col_map[chrom]
            plt.plot(plot_positions, log_ps, ".", markersize=markersize, alpha=0.7, color=color, mew=0)

        chrom_end = chrom_res_dict[chrom]['chrom_end'] / 1000000.0
        for j in sp.arange(offset, offset + chrom_end, 10):
            tick_positions.append(j)
            if  j % 20 == 0 and j < chrom_end - 10 :
                tick_strings.append(j)
            else:
                tick_strings.append("")

        offset += chrom_end + 1  # one Mb buffer

    if plot_bonferroni:
        if not b_threshold:
            b_threshold = -sp.log10(1.0 / (num_snps * 20.0))
        else:
            plt.plot([0, offset], [b_threshold, b_threshold], color='#000000', linestyle="-.")

    x_range = offset
    plt.axis([-x_range * 0.01, x_range * 1.01, -0.05 * max_log_pval, 1.05 * max_log_pval])
    plt.xticks(tick_positions, tick_strings, fontsize='x-small')
    plt.ylabel('$ - log(p - $value$)$')

#    plt.xlabel("Chromosome")
#        else:
    plt.xlabel("Mb")


    if png_file:
        plt.savefig(png_file, format="png", dpi=300, bbox_inches='tight')

    plt.clf()
    plt.close()

    
