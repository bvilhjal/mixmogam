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
    num_causals = 20
    # Simulating (unlinked) genotypes and phenotype pairs w. random positive correlation
    d = simulations.get_simulated_data(num_indivs=num_indivs, num_snps=num_snps,
                                       num_trait_pairs=num_trait_pairs, num_causals=num_causals)
    
    for i in range(num_trait_pairs):
        # The two different phenotypes.
        phen1 = d['trait_pairs'][i][0]
        phen2 = d['trait_pairs'][i][1]
        # Stacking up the two phenotypes into one vector.
        Y = sp.hstack([phen1, phen2])
        
        # The higher genetic correlation, the better the model fit (since we assume genetic correlation is 1).
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
        highlight_loci = sp.array(sd.get_chr_pos_list())[d['causal_indices_list'][i]]
        highlight_loci = highlight_loci.tolist()
        highlight_loci.sort()

        # Plotting stuff
        res = gr.Result(scores=gtres['ps'], snps_data=sd)
        res.plot_manhattan(png_file='%s_%d_gtres_manhattan.png' % (plot_prefix , i),
                           percentile=50,highlight_loci=highlight_loci,
                           plot_bonferroni=True,
                           neg_log_transform=True)
        res.plot_qq('%s_%d_gtres_qq.png' % (plot_prefix , i))
        res = gr.Result(scores=gtgres['ps'], snps_data=sd)
        res.plot_manhattan(png_file='%s_%d_gtgres_manhattan.png' % (plot_prefix , i),
                           percentile=50, highlight_loci=highlight_loci,
                           plot_bonferroni=True,
                           neg_log_transform=True)
        res.plot_qq('%s_%d_gtgres_qq.png' % (plot_prefix , i))
        res = gr.Result(scores=gres['ps'], snps_data=sd)
        res.plot_manhattan(png_file='%s_%d_gres_manhattan.png' % (plot_prefix , i),
                           percentile=50, highlight_loci=highlight_loci,
                           plot_bonferroni=True,
                           neg_log_transform=True)
        res.plot_qq('%s_%d_gres_qq.png' % (plot_prefix , i))
        
        
        
        
        
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

