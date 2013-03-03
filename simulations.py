"""
Code used for simulations of genotypes and phenotypes.

This is potentially useful for debugging.

Author: Bjarni J. Vilhjalmsson
Email: bjarni.vilhjalmsson@gmail.com
"""

import scipy as sp
from scipy import stats
import sys
import random
import snpsdata


def simulate_genotypes(num_indivs=1000, num_snps=1000):
    """
    Simulate genotypes
    """
    snps = sp.random.random((num_snps, num_indivs))
    snps = sp.round_(snps)
    snps = sp.array(snps, dtype='int8')
    snps = snps[sp.sum(snps, 1) > 0]
    positions = []
    chromosomes = []
    for chrom in range(1, 6):
        len_chrom = num_snps / 5
        chromosomes.extend([chrom] * (len_chrom))
        positions.extend(range(1, len_chrom + 1))
    indiv_ids = range(num_indivs)
    sd = snpsdata.construct_snps_data_set(snps, positions, chromosomes, indiv_ids)
    print '\nDone Simulating SNPs.'
    return sd


def simulate_traits(sd, num_trait_pairs=1, h2=0.8, num_causal=100, model='exponential'):
    """
    Return a dictionary object of correlated trait pairs.
    
    The causal markers are randomly chosen, with exponentially distributed effect sizes.
    
    h2: heritability
    """

    snps = sd.getSnps()
    cpl = sd.getChrPosList()
    num_snps = len(snps)
    num_accessions = len(sd.accessions)

    trait_pairs = []
    sample_indices_list = []
    snp_effects_list = []
    num_non_overlapping_list = []
    rho_est_list = []
    trait1_perc_var_list = []
    trait2_perc_var_list = []
    for i in range(num_trait_pairs):
        sys.stdout.write('\b\b\b\b\b\b\b%0.2f%%' % (100.0 * (min(1, (i + 1) / num_trait_pairs))))
        sys.stdout.flush()
            # print (i + 1) / int(num_traits / 100)
        # Simulate trait pair..

        num_non_overlapping = int(round(stats.beta.rvs(4, 2) * num_causal))
        num_non_overlapping_list.append(num_non_overlapping)
        num_causatives = num_causal + num_non_overlapping
        sample_indices = random.sample(range(num_snps), num_causatives)
        chosen_snps = sp.array([snps[i] for i in sample_indices])
        c = sp.random.random_integers(0, 1, (num_causatives, 1))
        chosen_snps = sp.absolute(c - chosen_snps)
        exp_effects = stats.expon.rvs(scale=1, size=(num_causatives, 1))
        # exp_effects = stats.norm.rvs(scale=1, size=(num_causatives, 1))
        snp_effects = chosen_snps * exp_effects
        snp_effects1 = snp_effects[:num_causal]
        snp_effects2 = snp_effects[-num_causal:]
        trait1 = sp.sum(snp_effects1, 0)
        trait2 = sp.sum(snp_effects2, 0)

        gv = sp.var(trait1, ddof=1)
        error = stats.norm.rvs(0, 1, size=num_accessions)
        ev = sp.var(error, ddof=1)
        n_trait1 = trait1 + error * sp.sqrt(((1.0 - h2) / h2) * (gv / ev))
        trait1_perc_var_list.append(sp.var(snp_effects1, 1) / sp.var(n_trait1))
        n_trait1 = (n_trait1 - sp.mean(n_trait1)) / sp.std(n_trait1)
        gv = sp.var(trait2, ddof=1)
        error = stats.norm.rvs(0, 1, size=num_accessions)
        ev = sp.var(error, ddof=1)
        n_trait2 = trait2 + error * sp.sqrt(((1.0 - h2) / h2) * (gv / ev))
        trait2_perc_var_list.append(sp.var(snp_effects2, 1) / sp.var(n_trait2))
        n_trait2 = (n_trait2 - sp.mean(n_trait2)) / sp.std(n_trait2)
        trait_pairs.append((n_trait1, n_trait2))
        sample_indices_list.append(sample_indices)
        snp_effects_list.append(exp_effects)
        rho_est = sp.corrcoef(trait1, trait2)
        rho_est_list.append(rho_est)

        # Variance contributions.

    print '\nDone simulating traits'
    sim_traits_dict = {'trait_pairs':trait_pairs, 'sample_indices_list':sample_indices_list,
            'num_non_overlapping_list':num_non_overlapping_list, 'snp_effects_list':snp_effects_list,
            'rho_est_list':rho_est_list, 'trait1_perc_var_list':trait1_perc_var_list,
            'trait2_perc_var_list':trait2_perc_var_list}
    return sim_traits_dict


def get_simulated_data(num_indivs=1000,num_snps=10000,num_trait_pairs=10):
    sd = simulate_genotypes(num_indivs=num_indivs, num_snps=num_snps)
    trait_dict = simulate_traits(sd, num_trait_pairs=num_trait_pairs)
    trait_dict['sd'] = sd
    return trait_dict




if __name__ == '__main__':
    pass
