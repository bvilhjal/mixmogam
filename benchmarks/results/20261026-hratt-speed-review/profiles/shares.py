"""Where HRATT's time goes before and after the 26 October speed branch.

usage: python shares.py   (from this directory; reads before/ and after/)

Each profile is one cProfile run of one fit, one thread, produced by
profile_hratt.py: `p3` is the 5,000 x 400,000 founder-block panel of
`20261025T153917Z-kvik-scale`, `hapnest50k` the 50,000 x 20,000 HAPNEST panel
of `20261003-hapnest-kvik-n50000-m20000`. Times are cumulative seconds
(own seconds for decode_table and VBEngine.fit's own row); cProfile inflates
small calls, so shares are indicative. decode_table is called from several
of the other rows, so the rows overlap.
"""
import json
import pstats

NAMES = {
    "hratt fit scores (hold-out + cross-fit fits, Grams)": ("_hratt_fit_scores", "twostep", "cum"),
    "VBEngine.fit, cumulative": ("fit", "_vb", "cum"),
    "VBEngine.fit, own time (GEMMs, sweeps, refresh)": ("fit", "_vb", "tot"),
    "joint Gram preparation": ("_prepare_crossfit_grams", "_vb", "cum"),
    "decode_table (all callers)": ("decode_table", "_standardize", "tot"),
    "prepare_moments": ("prepare_moments", "_standardize", "cum"),
    "_he_alpha": ("_he_alpha", "twostep", "cum"),
    "variance fit (_hratt_variance)": ("_hratt_variance", "twostep", "cum"),
    "Lanczos quadrature (REML)": ("lanczos_quadrature_batch", "_slq", "cum"),
    "predict": ("predict", "_vb", "cum"),
    "_out_of_fold_loco": ("_out_of_fold_loco", "twostep", "cum"),
    "_retro_stats": ("_retro_stats", "twostep", "cum"),
}


def read(path):
    st = pstats.Stats(path)
    out = {}
    for (fn, _, nm), (_, _, tt, ct, _) in st.stats.items():
        for key, (name, mod, kind) in NAMES.items():
            if nm == name and mod in fn:
                out[key] = out.get(key, 0.0) + (ct if kind == "cum" else tt)
    return out, st.total_tt


for label in ("p3-he", "p3-reml", "hapnest50k-he", "hapnest50k-reml"):
    old, t_old = read(f"before/{label}.pstats")
    new, t_new = read(f"after/{label}.pstats")
    iters = json.load(open(f"after/{label}.json"))["extra"].get("loco_iterations")
    print(f"\n{label}: {t_old:.1f} s -> {t_new:.1f} s (LOCO sweeps {iters})")
    for key in NAMES:
        print(f"  {key:52s} {old.get(key, 0):6.2f} -> {new.get(key, 0):6.2f}")
