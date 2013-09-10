"""
A container for functions which aim to analyze or process gwas results, for some aim.


Author: Bjarni J. Vilhjalmsson
Email: bjarni.vilhjalmsson@gmail.com
"""

import scipy as sp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
import itertools as it
import random
import phenotypeData as pd
import math




def calc_median(scores, exp_median=0.5):
    s = sp.copy(scores)
    s.sort()
    median = s[len(s) / 2]
    del s
    return (exp_median - median)

def _estAreaBetweenCurves_(quantiles, expQuantiles):
    area = 0
    for i in range(0, len(quantiles) - 1):
        area += (expQuantiles[i + 1] - expQuantiles[i]) * (abs(quantiles[i + 1] - expQuantiles[i + 1] + quantiles[i] - expQuantiles[i])) / 2.0
    return area

def calc_ks_stats(scores, exp_scores=None):
    from scipy import stats
    if exp_scores:
        (D, p_val) = stats.ks_2samp(scores, exp_scores)
    else:
        (D, p_val) = stats.kstest(scores, stats.uniform.cdf)
    return {'D':D, 'p_val':p_val}

def _getExpectedPvalueQuantiles_(numQuantiles):
    quantiles = []
    for i in range(numQuantiles):
        quantiles.append(float(i) + 0.5 / (numQuantiles + 1))
    return quantiles




def get_log_quantiles(scores, num_dots=1000, max_val=5):
    """
    Uses scipy
    """
    scores = sp.copy(sp.array(scores))
    scores.sort()
    indices = sp.array(10 ** ((-sp.arange(1, num_dots + 1, dtype='single') / (num_dots + 1)) * max_val) \
                * len(scores), dtype='int')
    return -sp.log10(scores[indices])



def simple_log_qqplot(quantiles_list, png_file=None, pdf_file=None, quantile_labels=None, line_colors=None,
            max_val=5, title=None, text=None, plot_label=None, ax=None, **kwargs):
    storeFig = False
    if ax is None:
        f = plt.figure(figsize=(5.4, 5))
        ax = f.add_axes([0.1, 0.09, 0.88, 0.86])
        storeFig = True
    ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.5, linewidth=2.0)
    num_dots = len(quantiles_list[0])
    exp_quantiles = sp.arange(1, num_dots + 1, dtype='single') / (num_dots + 1) * max_val
    for i, quantiles in enumerate(quantiles_list):
        if line_colors:
            c = line_colors[i]
        else:
            c = 'b'
        if quantile_labels:
            ax.plot(exp_quantiles, quantiles, label=quantile_labels[i], c=c, alpha=0.5, linewidth=2.2)
        else:
            ax.plot(exp_quantiles, quantiles, c=c, alpha=0.5, linewidth=2.2)
    ax.set_ylabel("Observed $-log_{10}(p$-value$)$")
    ax.set_xlabel("Expected $-log_{10}(p$-value$)$")
    if title:
        ax.title(title)
    max_x = max_val
    max_y = max(map(max, quantiles_list))
    ax.axis([-0.025 * max_x, 1.025 * max_x, -0.025 * max_y, 1.025 * max_y])
    if quantile_labels:
        fontProp = matplotlib.font_manager.FontProperties(size=10)
        ax.legend(loc=2, numpoints=2, handlelen=0.05, markerscale=1, prop=fontProp, pad=0.018)
    y_min, y_max = plt.ylim()
    if text:
        f.text(0.05 * max_val, y_max * 0.9, text)
    if plot_label:
        f.text(-0.138 * max_val, y_max * 1.01, plot_label, fontsize=14)
    if storeFig == False:
        return
    if png_file != None:
        f.savefig(png_file)
    if pdf_file != None:
        f.savefig(pdf_file, format='pdf')



def get_quantiles(scores, num_dots=1000):
    """
    Uses scipy
    """
    scores = sp.copy(sp.array(scores))
    scores.sort()
    indices = [int(len(scores) * i / (num_dots + 2)) for i in range(1, num_dots + 1)]
    return scores[indices]



def simple_qqplot(quantiles_list, png_file=None, pdf_file=None, quantile_labels=None, line_colors=None,
            title=None, text=None, ax=None, plot_label=None, **kwargs):
    storeFig = False
    if ax is None:
        f = plt.figure(figsize=(5.4, 5))
        ax = f.add_axes([0.11, 0.09, 0.87, 0.86])
        storeFig = True
    ax.plot([0, 1], [0, 1], 'k--', alpha=0.5, linewidth=2.0)
    num_dots = len(quantiles_list[0])
    exp_quantiles = sp.arange(1, num_dots + 1, dtype='single') / (num_dots + 1)
    for i, quantiles in enumerate(quantiles_list):
        if line_colors:
            c = line_colors[i]
        else:
            c = 'b'
        if quantile_labels:
            ax.plot(exp_quantiles, quantiles, label=quantile_labels[i], c=c, alpha=0.5, linewidth=2.2)
        else:
            ax.plot(exp_quantiles, quantiles, c=c, alpha=0.5, linewidth=2.2)
    ax.set_ylabel("Observed $p$-value")
    ax.set_xlabel("Expected $p$-value")
    if title:
        ax.title(title)
    ax.axis([-0.025, 1.025, -0.025, 1.025])
    if quantile_labels:
        fontProp = matplotlib.font_manager.FontProperties(size=10)
        ax.legend(loc=2, numpoints=2, handlelen=0.05, markerscale=1, prop=fontProp, pad=0.018)
    if text:
        f.text(0.05, 0.9, text)
    if plot_label:
        f.text(-0.151, 1.04, plot_label, fontsize=14)
    if storeFig == False:
        return
    if png_file != None:
        f.savefig(png_file)
    if pdf_file != None:
        f.savefig(pdf_file, format='pdf')



def plot_simple_qqplots(png_file_prefix, results, result_labels=None, line_colors=None,
            num_dots=1000, title=None, max_neg_log_val=5):
    """
    Plots both log QQ-plots and normal QQ plots.
    """
    qs = []
    log_qs = []
    for res in results:
        pvals = res.snp_results['scores'][:]
        qs.append(get_quantiles(pvals, num_dots))
        log_qs.append(get_log_quantiles(pvals, num_dots, max_neg_log_val))
    simple_qqplot(qs, png_file_prefix + '_qq.png', quantile_labels=result_labels,
                line_colors=line_colors, num_dots=num_dots, title=title)
    simple_log_qqplot(log_qs, png_file_prefix + '_log_qq.png', quantile_labels=result_labels,
                line_colors=line_colors, num_dots=num_dots, title=title, max_val=max_neg_log_val)

def plot_simple_qqplots_pvals(png_file_prefix, pvals_list, result_labels=None, line_colors=None,
            num_dots=1000, title=None, max_neg_log_val=5):
    """
    Plots both log QQ-plots and normal QQ plots.
    """
    qs = []
    log_qs = []
    for pvals in pvals_list:
        qs.append(get_quantiles(pvals, num_dots))
        log_qs.append(get_log_quantiles(pvals, num_dots, max_neg_log_val))
    simple_qqplot(qs, png_file_prefix + '_qq.png', quantile_labels=result_labels,
                line_colors=line_colors, num_dots=num_dots, title=title)
    simple_log_qqplot(log_qs, png_file_prefix + '_log_qq.png', quantile_labels=result_labels,
                line_colors=line_colors, num_dots=num_dots, title=title, max_val=max_neg_log_val)




if __name__ == '__main__':
    pass
