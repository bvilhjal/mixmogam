"""Super-efficient mixed linear models for genome-wide association mapping."""

__version__ = "2.0.0.dev6"

__all__ = [
    "LMM",
    "LMFit",
    "Genotypes",
    "Phenotypes",
    "GwasResult",
    "gwas",
    "kinship",
    "lmm",
    "pca",
    "simulate",
    "__version__",
]

_LAZY = {
    "LMM": ("mixmogam.lmm", "LMM"),
    "LMFit": ("mixmogam.lmm", "LMFit"),
    "Genotypes": ("mixmogam.genotypes", "Genotypes"),
    "Phenotypes": ("mixmogam.phenotypes", "Phenotypes"),
    "GwasResult": ("mixmogam.results", "GwasResult"),
    "gwas": ("mixmogam.association", "gwas"),
    "kinship": ("mixmogam.kinship", None),
    "lmm": ("mixmogam.lmm", None),
    "pca": ("mixmogam.pca", None),
    "simulate": ("mixmogam.simulate", None),
}


def __getattr__(name):
    if name in _LAZY:
        import importlib

        module_name, attr = _LAZY[name]
        module = importlib.import_module(module_name)
        return module if attr is None else getattr(module, attr)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(__all__)
