"""
Examples for how to perform GWAS using mixed models, and stepwise mixed models.

Author: Bjarni J. Vilhjalmsson
Email: bjarni.vilhjalmsson@gmail.com

"""


def load_a_thaliana_genotypes():
    """
    Loads A. thaliana genotypes (Horton et al., 2012) and returns a snps_data object
    """
    import dataParsers as dp
    sd = dp.parse_snp_data('at_data/all_chromosomes_binary.csv')
    return sd


def load_a_thaliana_phenotypes():
    """
    Loads A. thaliana phenotypes (Atwell et al., 2010) and returns a phenotype_data 
    object containing 107 different phenotypes.
    """
    import phenotypeData as pd
    phend = pd.parse_phenotype_file('at_data/199_phenotypes.csv')
    return phend


def linear_regression_gwas(phenotype_id=5, pvalue_file='lr_results.pvals',
                           manhattan_plot_file='lr_manhattan.png',
                           qq_plot_file_prefix='lr_qq'):
    """
    Perform linear regression GWAS for flowering time (phenotype_id=5 in the phenotype file) 
    in plants grown under 10C conditions. 
    """
    import linear_models as lm
    import gwaResults as gr
    # Load genotypes
    sd = load_a_thaliana_genotypes()

    # Load phenotypes
    phend = load_a_thaliana_phenotypes()

    # Coordinate phenotype of interest and genotypes.  This filters the genotypes and
    # phenotypes, leaving only accessions (individuals) which overlap between both,
    # and SNPs that are polymorphic in the resulting subset.
    sd.coordinate_w_phenotype_data(phend, phenotype_id)

    # Perform linear regression GWAS
    lr_results = lm.linear_model(sd.get_snps(), phend.get_values(phenotype_id))

    # Construct a results object
    res = gr.Result(scores=lr_results['ps'], snps_data=sd)

    # Save p-values to file
    res.write_to_file(pvalue_file)

    # Plot Manhattan plot
    res.plot_manhattan(png_file=manhattan_plot_file, percentile=90, plot_bonferroni=True,
                       neg_log_transform=True)
    # Plot a QQ-plot
    res.plot_qq(qq_plot_file_prefix)


def mixed_model_gwas(phenotype_id=5, pvalue_file='mm_results.pvals',
                     manhattan_plot_file='mm_manhattan.png',
                     qq_plot_file_prefix='mm_qq'):
    """
    Perform mixed model (EMMAX) GWAS for flowering time (phenotype_id=5 in the phenotype file) 
    in plants grown under 10C conditions. 
    """
    import linear_models as lm
    import kinship
    import gwaResults as gr
    # Load genotypes
    sd = load_a_thaliana_genotypes()

    # Load phenotypes
    phend = load_a_thaliana_phenotypes()

    # Coordinate phenotype of interest and genotypes.  This filters the genotypes and
    # phenotypes, leaving only accessions (individuals) which overlap between both,
    # and SNPs that are polymorphic in the resulting subset.
    sd.coordinate_w_phenotype_data(phend, phenotype_id)

    # Calculate kinship (IBS)
    K = kinship.calc_ibs_kinship(sd.get_snps())

    # Perform mixed model GWAS
    mm_results = lm.emmax(sd.get_snps(), phend.get_values(phenotype_id), K)

    # Construct a results object
    res = gr.Result(scores=mm_results['ps'], snps_data=sd)

    # Save p-values to file
    res.write_to_file(pvalue_file)

    # Plot Manhattan plot
    res.plot_manhattan(png_file=manhattan_plot_file, percentile=90, plot_bonferroni=True,
                       neg_log_transform=True)
    # Plot a QQ-plot
    res.plot_qq(qq_plot_file_prefix)


def multiple_loci_mixed_model_gwas(phenotype_id=5, pvalue_file_prefix='mlmm_results',
                                   result_files_prefix='mlmm_manhattan', max_num_steps=10, snp_priors=None):
    """
    Perform multiple loci mixed model GWAS for flowering time (phenotype_id=5 in the phenotype file) 
    in plants grown under 10C conditions. 
    """
    import linear_models as lm
    import kinship
    # Load genotypes
    sd = load_a_thaliana_genotypes()

    # Load phenotypes
    phend = load_a_thaliana_phenotypes()

    # Coordinate phenotype of interest and genotypes.  This filters the genotypes and
    # phenotypes, leaving only accessions (individuals) which overlap between both,
    # and SNPs that are polymorphic in the resulting subset.
    sd.coordinate_w_phenotype_data(phend, phenotype_id)

    # Calculate kinship (IBS)
    K = kinship.calc_ibs_kinship(sd.get_snps())

    # Perform multiple loci mixed model GWAS
    mlmm_results = lm.mlmm(phend.get_values(phenotype_id), K, sd=sd,
                           num_steps=max_num_steps, file_prefix=result_files_prefix,
                           save_pvals=True, pval_file_prefix=result_files_prefix, snp_priors=snp_priors)


def perform_cegs_gwas(kinship_type='ibd', phen_type='medians'):
    """
    Perform a simple MLM GWAS for the 8 traits
    """
    import hdf5_data
    import kinship
    import linear_models as lm
    import time
    import scipy as sp
    from matplotlib import pyplot as plt
    import analyze_gwas_results as agr
    phen_dict = hdf5_data.parse_cegs_drosophila_phenotypes()

    phenotypes = ['Protein', 'Sugar', 'Triglyceride', 'weight']
    envs = ['mated', 'virgin']
    for phenotype in phenotypes:
        for env in envs:
            print phenotype, env
            s1 = time.time()
            d = hdf5_data.coordinate_cegs_genotype_phenotype(
                phen_dict, phenotype, env)
            print 'Calculating kinship'
            if kinship_type == 'ibs':
                K = kinship.calc_ibs_kinship(d['snps'])
            elif kinship_type == 'ibd':
                K = kinship.calc_ibd_kinship(d['snps'])
            else:
                raise NotImplementedError

            if phen_type == 'means':
                lmm = lm.LinearMixedModel(d['Y_means'])
            elif phen_type == 'medians':
                lmm = lm.LinearMixedModel(d['Y_medians'])
            else:
                raise NotImplementedError
            lmm.add_random_effect(K)

            print "Running EMMAX"
            res = lmm.emmax_f_test(d['snps'], emma_num=1000)
            print 'Mean p-value:', sp.mean(res['ps'])

            secs = time.time() - s1
            if secs > 60:
                mins = int(secs) / 60
                secs = secs - mins * 60
                print 'Took %d mins and %f seconds.' % (mins, secs)
            else:
                print 'Took %f seconds.' % (secs)

            # Now generating QQ-plots
            label_str = '%s_%s_%s_%s' % (
                kinship_type, phenotype, env, phen_type)
            agr.plot_simple_qqplots_pvals('/Users/bjarnivilhjalmsson/data/tmp/cegs_qq_%s' % (label_str),
                                          [res['ps']], result_labels=[
                                              label_str], line_colors=['green'],
                                          num_dots=1000, title=None, max_neg_log_val=6)

            # Perform multiple loci mixed model GWAS
            chromosomes = d['positions'][:, 0]
            positions = sp.array(d['positions'][:, 1], 'int32')
            x_positions = []
            y_log_pvals = []
            colors = []
            x_shift = 0
            for i, chrom in enumerate(sp.unique(chromosomes)):
                if chrom in ['2L', '2LHet', '3L', '3LHet', '4', 'X', 'XHet']:
                    colors.append('c')
                else:  # chrom in ['2R', '2RHet', '3R', '3RHet', 'U', 'Uextra']
                    # Toss U and Hets
                    colors.append('m')
                chrom_filter = sp.in1d(chromosomes, chrom)
                positions_slice = positions[chrom_filter]
                x_positions.append(positions_slice + x_shift)
                x_shift += positions_slice.max()
                log_ps_slice = -sp.log10(res['ps'][chrom_filter])
                y_log_pvals.append(log_ps_slice)

            m = len(positions)
            log_bonf = -sp.log10(1 / (20.0 * m))
            print m, log_bonf

            # Plot manhattan plots?
            plt.figure(figsize=(12, 4))
            plt.axes([0.03, 0.1, 0.95, 0.8])
            for i, chrom in enumerate(sp.unique(chromosomes)):
                plt.plot(x_positions[i], y_log_pvals[i],
                         c=colors[i], ls='', marker='.')
            xmin, xmax = plt.xlim()
            plt.hlines(log_bonf, xmin, xmax, colors='k',
                       linestyles='--', alpha=0.5)
            plt.title('%s, %s' % (phenotype, env))
            plt.savefig('/Users/bjarnivilhjalmsson/data/tmp/cegs_gwas_%s_%s_%s_%s.png' %
                        (kinship_type, phenotype, env, phen_type))


def leave_k_out_blup(num_repeats=20, num_cvs=5, genotype_file='/Users/bjarnivilhjalmsson/data/cegs_lehmann/', k_thres=0.5):
    """

    """
    import h5py
    import hdf5_data
    import kinship
    import linear_models as lm
    import time
    import scipy as sp
    from matplotlib import pyplot as plt
    import analyze_gwas_results as agr
    phen_dict = hdf5_data.parse_cegs_drosophila_phenotypes()

    phenotypes = ['Protein', 'Sugar', 'Triglyceride', 'weight']
    envs = ['mated', 'virgin']
    rep_dict = {}
    for rep_i in range(num_repeats):
        res_dict = {}
        for phenotype in phenotypes:
            env_dict = {}
            for env in envs:
                print phenotype, env
                s1 = time.time()
                # Load data..
                d = hdf5_data.coordinate_cegs_genotype_phenotype(
                    phen_dict, phenotype, env, k_thres=k_thres)
                Y_means = d['Y_means']
                snps = d['snps']
                assert sp.all(sp.negative(sp.isnan(snps))), 'WTF?'
                K = kinship.calc_ibd_kinship(snps)
                print '\nKinship calculated'
                assert sp.all(sp.negative(sp.isnan(K))), 'WTF?'
                n = len(Y_means)
                # partition genotypes in k parts.
                gt_ids = d['gt_ids']
                num_ids = len(gt_ids)
                chunk_size = num_ids / num_cvs

                # Create k CV sets of prediction and validation data

                cv_chunk_size = int((n / num_cvs) + 1)
                ordering = sp.random.permutation(n)

                a = sp.arange(n)
                osb_ys = []
                pred_ys = []
                p_herits = []
                for cv_i, i in enumerate(range(0, n, cv_chunk_size)):
                    cv_str = 'cv_%d' % cv_i
                    # print 'Working on CV %d' % cv_i
                    end_i = min(n, i + cv_chunk_size)
                    validation_filter = sp.in1d(a, ordering[i:end_i])
                    training_filter = sp.negative(validation_filter)

                    train_snps = snps[:, training_filter]
                    val_snps = snps[:, validation_filter]

                    train_Y = Y_means[training_filter]
                    val_Y = Y_means[validation_filter]

                    #Calc. kinship
                    K_train = K[training_filter, :][:, training_filter]
                    K_cross = K[validation_filter, :][:, training_filter]
                    # Do gBLUP
                    lmm = lm.LinearMixedModel(train_Y)
                    lmm.add_random_effect(K_train)
                    r1 = lmm.get_REML()

                    # Now the BLUP.
                    y_mean = sp.mean(lmm.Y)
                    Y = lmm.Y - y_mean
                    p_herit = r1['pseudo_heritability']
                    p_herits.append(p_herit)
                    #delta = (1 - p_herit) / p_herit
            #        if K_inverse == None:
            #            K_inverse = K.I
            #        M = (sp.eye(K.shape[0]) + delta * K_inverse)
            #        u_blup = M.I * Y
                    M = sp.mat(p_herit * sp.mat(K_train) +
                               (1 - p_herit) * sp.eye(K_train.shape[0]))
                    u_mean_pred = sp.array(K_cross * (M.I * Y)).flatten()
                    osb_ys.extend(val_Y)
                    pred_ys.extend(u_mean_pred)
                corr = sp.corrcoef(osb_ys, pred_ys)[1, 0]
                print 'Correlation:', corr
                r2 = corr**2
                print 'R2:', r2
                mean_herit = sp.mean(p_herits)
                print 'Avg. heritability:', mean_herit
                env_dict[env] = {'R2': r2, 'obs_y': osb_ys,
                                 'pred_y': pred_ys, 'corr': corr, 'avg_herit': mean_herit}

            res_dict[phenotype] = env_dict
        rep_dict[rep_i] = res_dict
    res_hdf5_file = '/Users/bjarnivilhjalmsson/data/tmp/leave_%d_BLUP_results_kthres_%0.1f.hdf5' % (
        num_cvs, k_thres)
    h5f = h5py.File(res_hdf5_file)
    for rep_i in range(num_repeats):
        res_dict = rep_dict[rep_i]
        rep_g = h5f.create_group('repl_%d' % rep_i)
        for phenotype in phenotypes:
            phen_g = rep_g.create_group(phenotype)
            for env in envs:
                d = res_dict[phenotype][env]
                env_g = phen_g.create_group(env)
                env_g.create_dataset('R2',  data=[d['R2']])
                env_g.create_dataset('corr',  data=[d['corr']])
                env_g.create_dataset('obs_y',  data=d['obs_y'])
                env_g.create_dataset('pred_y',  data=d['pred_y'])
                env_g.create_dataset('avg_herit',  data=[d['avg_herit']])
    h5f.close()


def _test_GxE_mixed_model_gwas(num_indivs=1000, num_snps=10000, num_trait_pairs=10,
                               plot_prefix='/Users/bjarnivilhjalmsson/tmp/test'):
    """
    Test for the multiple environment mixed model

    Simulates correlated trait pairs with exponentially distributed effects. 
    """

    import simulations
    import kinship
    import scipy as sp
    import linear_models as lm
    import gwaResults as gr
    num_trait_pairs = 10
    num_indivs = 200
    num_snps = 10000
    # Number of causal SNPs per trait (in total there may be up to twice that,
    # depending on genetic correlation)
    num_causals = 10

    # Simulating (unlinked) genotypes and phenotype pairs w. random positive
    # correlation
    d = simulations.get_simulated_data(num_indivs=num_indivs, num_snps=num_snps,
                                       num_trait_pairs=num_trait_pairs, num_causals=num_causals)

    for i in range(num_trait_pairs):
        # The two different phenotypes.
        phen1 = d['trait_pairs'][i][0]
        phen2 = d['trait_pairs'][i][1]
        # Stacking up the two phenotypes into one vector.
        Y = sp.hstack([phen1, phen2])

        # The higher genetic correlation, the better the model fit (since we
        # assume genetic correlation is 1).
        print 'The genetic correlation between the two traits is %0.4f' % d['rho_est_list'][i][0, 1]

        # The genotypes
        sd = d['sd']
        snps = sd.get_snps()
        # Doubling the genotype data as well.
        snps = sp.hstack([snps, snps])

        # Calculating the kinship using the duplicated genotypes
        K = kinship.calc_ibd_kinship(snps)
        print ''

        # Calculating the environment vector
        E = sp.zeros((2 * num_indivs, 1))
        E[num_indivs:, 0] = 1

        print 'Here are the dimensions:'
        print 'Y.shape: ', Y.shape
        print 'snps.shape: ', snps.shape
        print 'E.shape: ', E.shape
        print 'K.shape: ', K.shape

        mm_results = lm.emmax_w_two_env(snps, Y, K, E)
        gtres = mm_results["gt_res"]
        gtgres = mm_results["gt_g_res"]
        gres = mm_results["g_res"]

        # Figuring out which loci are causal
        highlight_loci = sp.array(sd.get_chr_pos_list())[
            d['causal_indices_list'][i]]
        highlight_loci = highlight_loci.tolist()
        highlight_loci.sort()

        # Plotting stuff
        res = gr.Result(scores=gtres['ps'], snps_data=sd)
        res.plot_manhattan(png_file='%s_%d_gtres_manhattan.png' % (plot_prefix, i),
                           percentile=50, highlight_loci=highlight_loci,
                           plot_bonferroni=True,
                           neg_log_transform=True)
        res.plot_qq('%s_%d_gtres_qq.png' % (plot_prefix, i))
        res = gr.Result(scores=gtgres['ps'], snps_data=sd)
        res.plot_manhattan(png_file='%s_%d_gtgres_manhattan.png' % (plot_prefix, i),
                           percentile=50, highlight_loci=highlight_loci,
                           plot_bonferroni=True,
                           neg_log_transform=True)
        res.plot_qq('%s_%d_gtgres_qq.png' % (plot_prefix, i))
        res = gr.Result(scores=gres['ps'], snps_data=sd)
        res.plot_manhattan(png_file='%s_%d_gres_manhattan.png' % (plot_prefix, i),
                           percentile=50, highlight_loci=highlight_loci,
                           plot_bonferroni=True,
                           neg_log_transform=True)
        res.plot_qq('%s_%d_gres_qq.png' % (plot_prefix, i))


def lotus_data_analysis(phenotype_id=1,
                        result_files_prefix='/Users/bjarnivilhjalmsson/Dropbox/Cloud_folder/tmp/lmm_results',
                        manhattan_plot_file='/Users/bjarnivilhjalmsson/Dropbox/Cloud_folder/tmp/lmm_manhattan.png',
                        qq_plot_file_prefix='/Users/bjarnivilhjalmsson/Dropbox/Cloud_folder/tmp/lmm_qq'):
    """
    Lotus GWAS (data from Stig U Andersen)
    """
    import linear_models as lm
    import kinship
    import gwaResults as gr
    import dataParsers as dp
    import phenotypeData as pd

    # Load genotypes
    print 'Parsing genotypes'
    sd = dp.parse_snp_data(
        '/Users/bjarnivilhjalmsson/Dropbox/Lotus_GWAS/20140603_NonRep.run2.vcf.matrix.ordered.csv')

    # Load phenotypes
    print 'Parsing phenotypes'
    phend = pd.parse_phenotype_file(
        '/Users/bjarnivilhjalmsson/Dropbox/Lotus_GWAS/141007_FT_portal_upd.csv')

    print 'Box-cox'
    phend.box_cox_transform(1)

    # Coordinate phenotype of interest and genotypes.  This filters the genotypes and
    # phenotypes, leaving only accessions (individuals) which overlap between both,
    # and SNPs that are polymorphic in the resulting subset.
    print 'Coordinating data'
    sd.coordinate_w_phenotype_data(phend, phenotype_id)

    # Calculate kinship (IBS/IBD)
#     print 'Calculating kinship'
#     K = kinship.calc_ibd_kinship(sd.get_snps())
#     print K

    # Perform mixed model GWAS
    print 'Performing mixed model GWAS'
#     mm_results = lm.emmax(sd.get_snps(), phend.get_values(phenotype_id), K)

#     mlmm_results = lm.mlmm(phend.get_values(phenotype_id), K, sd=sd,
#                          num_steps=10, file_prefix=result_files_prefix,
# save_pvals=True, pval_file_prefix=result_files_prefix)

    lg_results = lm.local_vs_global_mm_scan(phend.get_values(phenotype_id), sd,
                                            file_prefix='/Users/bjarnivilhjalmsson/Dropbox/Cloud_folder/tmp/lotus_FT_loc_glob_0.1Mb',
                                            window_size=100000, jump_size=50000, kinship_method='ibd', global_k=None)

#     # Construct a results object
    print 'Processing results'
#     res = gr.Result(scores=mm_results['ps'], snps_data=sd)

    # Save p-values to file
#     res.write_to_file(pvalue_file)

    # Plot Manhattan plot
#     res.plot_manhattan(png_file=manhattan_plot_file, percentile=90, plot_bonferroni=True,
#                         neg_log_transform=True)
    # Plot a QQ-plot
#     res.plot_qq(qq_plot_file_prefix)

    # Local-global scan


def lotus_mixed_model_gwas(phenotype_id='OW_2014_15', phen_file='/home/bjarni/LotusGenome/cks/Lotus31012019/20181113_136LjAccessionData.csv',
                           gt_file='/home/bjarni/LotusGenome/cks/Lotus31012019/all_chromosomes_binary.csv',
                           pvalue_file='mm_results.pvals', manhattan_plot_file='mm_manhattan.png', qq_plot_file_prefix='mm_qq'):
    """
    Perform mixed model (EMMAX) GWAS for Lotus data
    """
    import linear_models as lm
    import kinship
    import gwaResults as gr
    import dataParsers as dp
    # Load genotypes
    sd = dp.parse_snp_data(gt_file)

    # Load phenotypes
    import phenotypeData as pd
    phend = pd.parse_phenotype_file(phen_file)

    # Coordinate phenotype of interest and genotypes.  This filters the genotypes and
    # phenotypes, leaving only accessions (individuals) which overlap between both,
    # and SNPs that are polymorphic in the resulting subset.
    sd.coordinate_w_phenotype_data(phend, phenotype_id)

    # Calculate kinship (IBS)
    K = kinship.calc_ibs_kinship(sd.get_snps())

    # Perform mixed model GWAS
    mm_results = lm.emmax(sd.get_snps(), phend.get_values(phenotype_id), K)

    # Construct a results object
    res = gr.Result(scores=mm_results['ps'], snps_data=sd)

    # Save p-values to file
    res.write_to_file(pvalue_file)

    # Plot Manhattan plot
    res.plot_manhattan(png_file=manhattan_plot_file, percentile=90, plot_bonferroni=True,
                       neg_log_transform=True)
    # Plot a QQ-plot
    res.plot_qq(qq_plot_file_prefix)


if __name__ == '__main__':
    #    lotus_data_analysis()
    #     _test_GxE_mixed_model_gwas()
    lotus_mixed_model_gwas()
    #    linear_regression_gwas()
#    multiple_loci_mixed_model_gwas()
    pass
