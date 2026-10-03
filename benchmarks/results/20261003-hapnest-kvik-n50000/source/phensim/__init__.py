"""phensim: genotype and phenotype simulators for genetic studies."""

__version__ = "1.0.0.dev1"

__all__ = [
    "simulate_independent",
    "simulate_population_structure",
    "simulate_haplotype_blocks",
    "simulate_ar1_blocks",
    "realistic_block_sizes",
    "simulate_coalescent",
    "simulate_by_mutation_rate",
    "simulate_hapnest",
    "iter_hapnest",
    "simulate_trait",
    "simulate_binary_trait",
    "simulate_confounded_trait",
    "simulate_gxe_trait",
    "simulate_correlated_traits",
    "ascertain_case_control",
    "n_eff_case_control",
    "h2_liability",
    "simulate_effects",
    "simulate_effects_pair",
    "genetic_correlation",
    "simulate_sumstats",
    "simulate_sumstats_pair",
    "prepare_blocks",
    "gwas_scan",
    "shake_ld",
    "simulate_pedigree",
    "pedigree_birth_times",
    "kinship_from_pedigree",
    "mendelian_draw",
    "grm",
    "ibs_kinship",
    "iter_loco_kinships",
    "loco_kinships",
    "windowed_kinships",
    "write_plink",
    "HAVE_NUMBA",
    "__version__",
]

_MAP = {
    "simulate_hapnest": ("phensim.hapnest", "simulate_hapnest"),
    "iter_hapnest": ("phensim.hapnest", "iter_hapnest"),
    "simulate_independent": ("phensim.genotypes", "simulate_independent"),
    "simulate_population_structure": ("phensim.genotypes", "simulate_population_structure"),
    "simulate_haplotype_blocks": ("phensim.genotypes", "simulate_haplotype_blocks"),
    "simulate_ar1_blocks": ("phensim.genotypes", "simulate_ar1_blocks"),
    "realistic_block_sizes": ("phensim.genotypes", "realistic_block_sizes"),
    "simulate_coalescent": ("phensim.genotypes", "simulate_coalescent"),
    "simulate_by_mutation_rate": ("phensim.genotypes", "simulate_by_mutation_rate"),
    "simulate_trait": ("phensim.phenotypes", "simulate_trait"),
    "simulate_binary_trait": ("phensim.phenotypes", "simulate_binary_trait"),
    "simulate_confounded_trait": ("phensim.phenotypes", "simulate_confounded_trait"),
    "simulate_gxe_trait": ("phensim.phenotypes", "simulate_gxe_trait"),
    "simulate_correlated_traits": ("phensim.phenotypes", "simulate_correlated_traits"),
    "ascertain_case_control": ("phensim.phenotypes", "ascertain_case_control"),
    "n_eff_case_control": ("phensim.phenotypes", "n_eff_case_control"),
    "h2_liability": ("phensim.phenotypes", "h2_liability"),
    "simulate_effects": ("phensim.sumstats", "simulate_effects"),
    "simulate_effects_pair": ("phensim.sumstats", "simulate_effects_pair"),
    "genetic_correlation": ("phensim.sumstats", "genetic_correlation"),
    "simulate_sumstats": ("phensim.sumstats", "simulate_sumstats"),
    "simulate_sumstats_pair": ("phensim.sumstats", "simulate_sumstats_pair"),
    "prepare_blocks": ("phensim.sumstats", "prepare_blocks"),
    "gwas_scan": ("phensim.sumstats", "gwas_scan"),
    "shake_ld": ("phensim.sumstats", "shake_ld"),
    "simulate_pedigree": ("phensim.pedigree", "simulate_pedigree"),
    "pedigree_birth_times": ("phensim.pedigree", "pedigree_birth_times"),
    "kinship_from_pedigree": ("phensim.pedigree", "kinship_from_pedigree"),
    "mendelian_draw": ("phensim.pedigree", "mendelian_draw"),
    "grm": ("phensim.kinship", "grm"),
    "ibs_kinship": ("phensim.kinship", "ibs_kinship"),
    "iter_loco_kinships": ("phensim.kinship", "iter_loco_kinships"),
    "loco_kinships": ("phensim.kinship", "loco_kinships"),
    "windowed_kinships": ("phensim.kinship", "windowed_kinships"),
    "write_plink": ("phensim.io", "write_plink"),
    "HAVE_NUMBA": ("phensim._numba", "HAVE_NUMBA"),
}


def __getattr__(name):
    if name in _MAP:
        import importlib

        module_name, attr = _MAP[name]
        return getattr(importlib.import_module(module_name), attr)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(__all__)
