"""
Examples for how to perform GWAS using mixed models, and stepwise mixed models.

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


def linear_regression_gwas(phenotype_id=5, manhattan_plot_file='lr_manhattan.png', 
                           qq_plot_file='lr_qq.png'):
    """
    Perform GWAS on flowering time (phenotype_id=5 in the phenotype file) 
    in plants grown under 10C conditions. 
    """
    import linear_models as lm
    import gwaResults as gr
    #Load genotypes
    sd = load_a_thaliana_genotypes()
    
    #Load phenotypes
    phend = load_a_thaliana_phenotypes()
    
    #Coordinate phenotype of interest and genotypes.  This filters the genotypes and 
    #phenotypes, leaving only accessions (individuals) which overlap between both, 
    #and SNPs that are polymorphic in the resulting subset.
    sd.coordinate_w_phenotype_data(phend, phenotype_id)
    
    #Perform linear regression GWAS
    lr_results = lm.linear_model(sd.get_snps(), phend.get_values(phenotype_id))
    
    #Construct a results object
    res = gr.Result(scores=lr_results['ps'], snps_data=sd)

    #Plot Manhattan plot
    res.plot_manhattan(png_file=manhattan_plot_file, percentile=90, plot_bonferroni=True,
                        neg_log_transform=True)
    #Plot a QQ-plot
    res.plot_qq(qq_plot_file)


def mixed_model_gwas():
    """
    GWAS 
    """
    pass


def multiple_loci_mixed_model_gwas():
    """
    GWAS 
    """
    pass


def load_human_genotypes():
    """
    """
    pass


def plot_manhattan_plot():
    """
    """
    pass


if __name__=='__main__':
    linear_regression_gwas()
