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
                     result_files_prefix='mlmm_manhattan', max_num_steps=10):
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
                         save_pvals=True, pval_file_prefix=result_files_prefix,)
            


def _test_GxE_mixed_model_gwas(num_indivs=1000, num_snps=10000, num_trait_pairs=10,
                               plot_prefix='/Users/bjarnivilhjalmsson/tmp/test'):
    """
    Test for the multiple environment mixed model
    """
    
    import simulations
    import kinship
    import scipy as sp
    import linear_models as lm
    import gwaResults as gr
    num_trait_pairs = 10
    num_indivs = 200
    num_snps = 100000
    # Simulating (unlinked) genotypes and phenotype pairs w. random positive correlation
    d = simulations.get_simulated_data(num_indivs=num_indivs, num_snps=num_snps, num_trait_pairs=num_trait_pairs)
    
    for i in range(num_trait_pairs):
        # The two different phenotypes.
        phen1 = d['trait_pairs'][0][0]
        phen2 = d['trait_pairs'][0][1]
        # Stacking up the two phenotypes into one vector.
        Y = sp.hstack([phen1, phen2])
        
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
        res = gr.Result(scores=gtres['ps'], snps_data=sd)
        res.plot_manhattan(png_file='%s_%d_gtres_manhattan.png' % (plot_prefix , i),
                           percentile=50,
                           plot_bonferroni=True,
                           neg_log_transform=True)
        res.plot_qq('%s_%d_gtres_qq.png' % (plot_prefix , i))
        res = gr.Result(scores=gtgres['ps'], snps_data=sd)
        res.plot_manhattan(png_file='%s_%d_gtgres_manhattan.png' % (plot_prefix , i),
                           percentile=50,
                           plot_bonferroni=True,
                           neg_log_transform=True)
        res.plot_qq('%s_%d_gtgres_qq.png' % (plot_prefix , i))
        res = gr.Result(scores=gres['ps'], snps_data=sd)
        res.plot_manhattan(png_file='%s_%d_gres_manhattan.png' % (plot_prefix , i),
                           percentile=50,
                           plot_bonferroni=True,
                           neg_log_transform=True)
        res.plot_qq('%s_%d_gres_qq.png' % (plot_prefix , i))
        
         


# import sys
#
# def load_dMelGt_genotypes(fn):
#    """
#    Load Dmel GTs
#    """
#    import dataParsers as dp
#    sd = dp.parse_snp_data(fn)
#    return sd
#
#
# def load_dMelGt_phenotypes(fn):
#    """
#    Loads Dmel Phenotyeps now
#    """
#    import phenotypeData as pd
#    phend = pd.parse_phenotype_file(fn)
#
#    return phend
#
# def mixed_mdl_gwas_GxE(fn_gt, fn_pht,
#                       phenotype_id=4, pvalue_file='mm_results.pvals',
#                       manhattan_plot_file='mm_manhattan.png',
#                       qq_plot_file_prefix='mm_qq'):
#    import linear_models as lm
#    import kinship
#    import gwaResults as gr
#
#    for i in xrange(1, 5):
#        sd = load_dMelGt_genotypes(fn_gt)
#        print sd.accessions
#        # Load phenotypes
#        phend = load_dMelGt_phenotypes(fn_pht)
#        print "phenotype: ",
#        print i
#        acc = sd.coordinate_w_phenotype_data(phend, i)
#        envcof = []
#        for rec in sd.accessions:
#            if rec.endswith("m"):
#                envcof.append(1)
#            else:
#                envcof.append(0)
#        K = kinship.calc_ibd_kinship(sd.get_snps())
#        # ADD COFACTORS HERE!!
#        # ARE THERE FUNCTIONS FOR ADDING THIS, WHAT DOES EMMAX expect
#        print len(envcof)
#        print len(sd.accessions)
#        print len(K)
#        mm_results = lm.emmax_w_two_env(sd.get_snps(), phend.get_values(i), K, envcof)
#        # print mm_results.keys()   
#        gtres = mm_results["gt_res"]
#        gtgres = mm_results["gt_g_res"]
#        gres = mm_results["g_res"]
#        # Construct a results object
#        # and plot manhattan plot
#        res = gr.Result(scores=gtres['ps'], snps_data=sd)
#        res.plot_manhattan(png_file=str(i) + "_gtres_" + manhattan_plot_file,
#                           percentile=90,
#                           plot_bonferroni=True,
#                           neg_log_transform=True)
#        res.plot_qq(str(i) + "_gtres_" + qq_plot_file_prefix)
#        res = gr.Result(scores=gtgres['ps'], snps_data=sd)
#        res.plot_manhattan(png_file=str(i) + "_gtgres_" + manhattan_plot_file,
#                           percentile=90,
#                           plot_bonferroni=True,
#                           neg_log_transform=True)
#        res.plot_qq(str(i) + "_gtgres_" + qq_plot_file_prefix)
#        res = gr.Result(scores=gres['ps'], snps_data=sd)
#        res.plot_manhattan(png_file=str(i) + "_gres_" + manhattan_plot_file,
#                           percentile=90,
#                           plot_bonferroni=True,
#                           neg_log_transform=True)
#        res.plot_qq(str(i) + "_gres_" + qq_plot_file_prefix)
#        # Plot a QQ-plot
#        #
#
#
#   
#
#
# def main(fn_gt, fn_pht):
#    mixed_mdl_gwas_GxE(fn_gt, fn_pht)
#    # mixed_model_gwas(fn_gt, fn_pht)
# if __name__ == '__main__':
#    pass


def load_plink_genotypes():
    """
    Loads genotypes in plink format.
    """
    raise NotImplementedError


def load_eigenstrat_genotypes():
    """
    Loads genotypes in eigenstrat format.
    """
    raise NotImplementedError



if __name__ == '__main__':
    _test_GxE_mixed_model_gwas()
    #    linear_regression_gwas()
#    multiple_loci_mixed_model_gwas()
    pass

