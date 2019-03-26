"""
Contains classes to handle results of GWAS

Author: Bjarni J. Vilhjalmsson
Email: bjarni.vilhjalmsson@gmail.com
"""


import pdb
import csv
import math
import itertools as it
try:
    import scipy as sp
#    import scipy.stats as st
except Exception, err_str:
    print 'scipy is missing:', err_str
import os
import cPickle
import bisect
#A dictionary for loaded results.. to avoid reloading.
#Use carefully to avoid memory leaks!



class ResultType(object):
    def __init__(self, resultType=None, fileType=None, datasetName=None, resultDir=None, mafCutoff=0, logTransform=None, name=None):
        self.resultType = resultType
        self.fileType = fileType
        self.datasetName = datasetName
        self.resultDir = resultDir
        self.mafCutoff = mafCutoff
        self.logTransform = logTransform
        if self.logTransform == None:
            self.logTransform = self.fileType == ".pvals"
        self.name = name
        if not name and datasetName:
            self.name = resultType + "_" + datasetName
        if resultDir and datasetName:
            self.filenamePrefix = resultDir + resultType + "_" + datasetName + "_"


    def getFileName(self, phed, phenotypeID, secondRun=False):
        print phenotypeID
        if secondRun:
            filename = self.filenamePrefix + str(phed.getPhenotypeName(phenotypeID)) + ".sr" + self.fileType
        else:
            filename = self.filenamePrefix + str(phed.getPhenotypeName(phenotypeID)) + self.fileType
        return filename

    def __str__(self):
        return self.resultType + "_" + self.datasetName




class Result(object):
    """
    Contains information on the result.  (The old object is renamed Result_old, for now..)  Poised to cause problems?
    """
    pickle_attributes = ['snp_results', 'keys', 'orders', 'ranks', 'phen_id', 'result_type', 'name',
                'accessions', 'chromosome_ends']

    def __init__(self, result_file=None, snp_results=None, scores=None, snps_data=None, accessions=None, name=None,
             result_type=None, phen_id=None, positions=None, chromosomes=None, mafs=None, macs=None,
             snps=None, auto_pickling_on=True, keys=None, **snp_results_info):
        """
        A simple init function.
        snps_data is assumed to match the results, (in term of positions, etc)
        
        accessions (when used) should fit the results!
        """

        self.snp_results = {}  #Main data container
        if snp_results:
            self.snp_results = snp_results
        self.phen_id = phen_id
        self.result_type = result_type
        self.name = name
        if keys:
            self.keys = keys
        else:
            self.keys = []
        self.accessions = accessions
        self.chromosome_ends = []
        self.orders = None
        self.ranks = None
        self.auto_pickling_on = auto_pickling_on


        if result_file:
            pickle_file = result_file + '.pickled'
            if self.auto_pickling_on and os.path.isfile(pickle_file):
                try:
                    with open(pickle_file) as f:
                        d = cPickle.load(f)
                    for attr in self.pickle_attributes:
                        setattr(self, attr, d[attr])
#                    all_keys = self.snp_results.keys()[:]
#                    for k in all_keys:
#                        if not k in ['chromosomes', 'positions', 'scores', 'mafs', 'macs']:
#                            del self.snp_results[k]
                    pickle_failed = False
                except Exception, err_str:
                    print 'Loading pickle file failed:', pickle_file
                    print err_str
                    self._load_result_(result_file)
                    pickle_failed = True

            else:
                self._load_result_(result_file)
        else:

            if chromosomes is not None:
                self.keys.append('chromosomes')
                self.snp_results['chromosomes'] = chromosomes
            if positions is not None:
                self.keys.append('positions')
                self.snp_results['positions'] = positions
            if scores is not None:
                self.keys.append('scores')
                self.snp_results['scores'] = scores
            if mafs is not None:
                self.keys.append('mafs')
                self.snp_results['mafs'] = mafs
            if macs is not None:
                self.keys.append('macs')
                self.snp_results['macs'] = macs
        if snps_data:
            self._load_snps_(snps_data)
            if 'mafs' not in self.snp_results:
                self._load_mafs_(snps_data)
        else:
            if snps:
                self.keys.append('snps')
                self.snp_results['snps'] = snps


        #Adding various info...
        if snp_results_info:
            for info in snp_results_info:
                self.keys.append(info)
                self.snp_results[info] = snp_results_info[info]

        if result_file and self.auto_pickling_on:
            if not os.path.isfile(pickle_file) or pickle_failed:
                d = {}
                for attr in self.pickle_attributes:
                    d[attr] = getattr(self, attr)
                with open(pickle_file, 'wb') as f:
                    cPickle.dump(d, f)





    def _load_result_(self, result_file, try_delims=[",", "\t", " ", ]):

        with open(result_file, "r") as f:
            print 'Loading results..'
            line = f.next()
            for delim in try_delims:
                if len(line.split(delim)) > 1:
                    break
            else:
                raise Exception("Appropriate delimiter wasn't found.")

            header = map(str.strip, line.split(delim))
            for key in header:  #Handling old data formats..
                key = key.lower()
                if key == 'chromosome':
                    key = 'chromosomes'
                if key == 'position':
                    key = 'positions'
                if key == 'score':
                    key = 'scores'
                if key == 'marf':
                    key = 'mafs'
                if key == 'maf':
                    key = 'macs'
                self.keys.append(key)

            for key in self.keys:
                self.snp_results[key] = []

            chrom = 0
            chrom_splits = []
            for i, line in enumerate(f):
                l = map(float, map(str.strip, line.split(delim)))
                for v, key in it.izip(l, self.keys):
                    self.snp_results[key].append(v)
                if int(l[0]) != chrom:
                    chrom_splits.append(i)
                    chrom = int(l[0])

            for si in chrom_splits[1:]:
                self.chromosome_ends.append(int(self.snp_results['positions'][si - 1]))
            self.chromosome_ends.append(int(self.snp_results['positions'][-1]))
            self._rank_scores_()



    def _load_snps_(self, snps_data):
        """
        This only loads SNPs...
        """
        self.snp_results['snps'] = snps_data.getSnps()
        self.snp_results['positions'] = snps_data.getPositions()
        self.snp_results['chromosomes'] = snps_data.get_chr_list()
        self.keys.append('chromosomes')
        self.keys.append('positions')
        self.keys.append('snps')
        self.chromosome_ends = snps_data.get_chromosome_ends()


    def _load_mafs_(self, snps_data):
        print 'Loading mafs'
        snps = self.snp_results['snps']
        if snps_data.data_format == 'float':
            self.snp_results['mafs'] = [1] * len(snps)
            self.snp_results['macs'] = [1] * len(snps)
        else:
            maf_d = snps_data.get_mafs()
            self.snp_results['mafs'] = maf_d['marfs']
            self.snp_results['macs'] = maf_d['mafs']
        self.keys.append('mafs')
        self.keys.append('macs')

    def _rank_scores_(self):
        """
        Generates two data structures:
        self.orders (SNP indices ordered by rank)
        self.ranks (ranks of the SNPs)
        """
        scores = self.snp_results['scores']
        rank_ls = zip(scores, range(len(scores)))
        rank_ls.sort()#reverse=True)
        self.orders = []
        for j in range(len(rank_ls)):
            (s, i) = rank_ls[j]
            self.orders.append(i)
            rank_ls[j] = (i, j)

        rank_ls.sort()
        self.ranks = []
        for (i, j) in rank_ls:
            self.ranks.append(j + 1)


    def _sort_by_chr_pos_(self):
        res_ls = zip(self.snp_results['chromosomes'], self.snp_results['positions'],
            range(len(self.snp_results['positions'])))
        res_ls.sort()
        orders = list(zip(*res_ls))[2]

        snp_results = {}
        for key in self.keys:
            snp_results[key] = []
            for i in orders:
                snp_results[key].append(self.snp_results[key][i])
        self.snp_results = snp_results


    def get_quantile(self, quantile):
        if not self.orders:
            self._rank_scores_()
        rank = int(self.num_scores() * quantile)
        return self.snp_results['scores'][self.orders[rank]]



#    def candidate_gene_enrichments(self, cgl=None, cgl_file=None, pval_thresholds=[0.01], gene_radius=20000,
#                methods=['chi_square'], num_perm=500, file_prefix=None,
#                obs_genes_file=None, early_stop_threshold=25, all_genes=None, cand_gene_indices=None):
#        """
#        Performs CGR analysis on this results object.
#        
#        cgl is a list of genes.
#        cgl_file is a file with a list of genes.
#        
#        method: 
#            chi_square - statistics
#            multinomial - statistics
#            gene_perm - permute genes.
#            snps_perm - permute SNPs
#        """
#        import pylab
#        import analyze_gene_enrichment as genr
#
#        if not 'chi_square' in methods:
#            methods.append('chi_square')
#
#        chrom_ends = self.get_chromosome_ends()
#        print 'chrom_ends', chrom_ends
#
#        #Parse cgl file
#        if cgl_file:
#            cgl, cg_tair_ids = load_cand_genes_file(cgl_file)
#        for cg in cgl:
#            print str(cg)
#
#
#        if not all_genes:
#            #Load genes from DB.
#            print 'Fetching all genes'
#            all_genes = get_gene_list(include_intron_exons=False, verbose=False)
#            print 'Fetched %d genes.' % len(all_genes)
#            num_genes = len(all_genes)
#
#        if not cand_gene_indices:
#            #Pre-process cgl.
#            cand_gene_indices = []
#            for i, g in enumerate(all_genes):
#                for cg in cgl:
#                    if g.dbRef == cg.dbRef:
#                        #print g.dbRef, cg.dbRef
#                        cand_gene_indices.append(i)
#                        break
#            num_cand_genes = len(cand_gene_indices)
#            #print num_cand_genes, cand_gene_indices
#
#
#        method_res_dict = {}
#        for m in methods:
#            method_res_dict[m] = {'statistics':[], 'pvals':[]}
#
#        if obs_genes_file:
#            obs_gene_str = ''
#
#        pval_thresholds.sort()
#        last_thres = 1.0
#        log_enrichments = []
#        for pval_threshold in reversed(pval_thresholds):
#            print 'Using p-value threshold % f' % pval_threshold
#            thres = pval_threshold / last_thres
#            print 'Using corrected threshold % f' % thres
#            last_thres = pval_threshold
#
#            #Filter pvalue file
#            self.filter_percentile(1 - thres)
#
#            #pre-process pvalues
#            regions = self.get_regions(gene_radius=gene_radius)
#
#            #Calculate observed candidate gene enrichment. 
#            obs_enrichments = genr.calc_enrichment(all_genes, cand_gene_indices, regions)
#            r1 = obs_enrichments[0] / float(obs_enrichments[1])
#            r2 = (num_cand_genes / float(num_genes))
#            obs_stat = sp.log(r1 / r2)
#            print 'Observed statistics % f' % obs_stat
#            log_enrichments.append(obs_stat)
#
#            #What cand. genes overlap with regions?
#            obs_cg_indices = obs_enrichments[2]
#            if obs_cg_indices:
#                obs_gene_str += str(pval_threshold) + ','
#                tair_ids = [all_genes[cgi].tairID for cgi in obs_cg_indices]
#                obs_gene_str += ','.join(tair_ids)
#                obs_gene_str += '\n'
#
#
#
#            for method in methods:
#                if method == 'chi_square':
#                    chi_sq_pval, chi_sq_stat = genr.get_chi_square_pval(obs_enrichments[0],
#                                            obs_enrichments[1],
#                                            num_cand_genes, num_genes)
#                    method_res_dict[method]['statistics'].append(chi_sq_stat)
#                    method_res_dict[method]['pvals'].append(chi_sq_pval)
#
#                if method == 'multinomial':
#                    pass
#                if method == 'gene_perm':
#                    p_val, perm_stats = genr.get_gene_perm_pval(obs_stat, regions, all_genes,
#                                    cand_gene_indices, num_perm=num_perm,
#                                    early_stop_threshold=early_stop_threshold)
#                    method_res_dict[method]['statistics'].append(perm_stats)
#                    method_res_dict[method]['pvals'].append(p_val)
#
#                if method == 'snps_perm':
#                    p_val, perm_stats = genr.get_snps_perm_pval(obs_stat, regions, all_genes,
#                                    cand_gene_indices, chrom_ends,
#                                    num_perm=num_perm,
#                                    early_stop_threshold=early_stop_threshold)
#
#                    method_res_dict[method]['statistics'].append(perm_stats)
#                    method_res_dict[method]['pvals'].append(p_val)
#
##                        h_res = pylab.hist(perm_stats)
##                        pylab.vlines(obs_stat, 0, max(h_res[0]), colors='r')
##                        pylab.savefig(env.env['tmp_dir'] + 'test.pdf', format='pdf')
#
#
#        if obs_genes_file:
#            with open(obs_genes_file, 'w') as f:
#                f.write(obs_gene_str)
#
#        #Now the plotting of the results.
#        method_name_dict = {'chi_square':'Chi-square test', 'gene_perm':'Candidate gene permutations',
#                    'snps_perm':'SNP positions permutation (chromosome rotation)' }
#
#
#        pval_thresholds.reverse()
#        pos_list = range(len(pval_thresholds))
#        for m in methods:
#            pylab.figure()
#            pvals = []
#            for p in method_res_dict[m]['pvals']:
#                if p != 0:
#                    pvals.append(p)
#                else:
#                    pvals.append(0.5 / num_perm)
#            neg_log_pvals = map(lambda x:-math.log10(x), pvals)
#            pylab.barh(pos_list, neg_log_pvals, align='center', color='g', alpha=0.6)
#            pylab.ylim((-1, len(pval_thresholds)))
#            pylab.yticks(pos_list, map(str, pval_thresholds))
#            pylab.xlabel('Enrichment -log(p-value)')
#            pylab.ylabel('p-value percentile threshold')
#            ymin, ymax = pylab.ylim()
#            xmin, xmax = pylab.xlim()
#            pylab.axvline(-math.log10(0.05), ymin=ymin, ymax=ymax, color='r')
#            pylab.xlim((0, max(-math.log10(0.05), max(neg_log_pvals)) * 1.05))
#            pylab.title(method_name_dict[m])
#            pylab.savefig(file_prefix + '_' + m + '.png', format='png')
#
#
#        pylab.figure()
#        pylab.barh(pos_list, log_enrichments, align='center', color='g', alpha=0.6)
#        pylab.ylim((-1, len(pval_thresholds)))
#        pylab.yticks(pos_list, map(str, pval_thresholds))
#        pylab.xlabel('Enrichment (log[ratio])')
#        pylab.ylabel('p-value percentile threshold')
#        pylab.savefig(file_prefix + '_enrichment_ratio.png', format='png')
#
#        return {'enr_stats':log_enrichments, 'method_res_dict':method_res_dict}
#
#

    def _plot_small_manhattan_(self, pdf_file=None, png_file=None, min_score=0, max_score=None,
                score_type="pvals", ylab="$-$log$_{10}(p-$value$)$", plot_bonferroni=False,
                cand_genes=None, threshold=0, highlight_markers=None, chromosome=None,
                tair_file=None, plot_genes=False, color_map=None, markersize=3, fig_size=(8, 5)):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        tair_ids = None
        min_x = min(self.snp_results['positions'])
        max_x = max(self.snp_results['positions'])
        if plot_genes:
            print 'Retrieving genes from DB.'
            gene_buffer = int((max_x - min_x) / 16) #bp
            gene_dict = get_gene_dict(chromosome, min_x, max_x, only_genes=True)
            tair_ids = sorted(gene_dict.keys())
            print "Found", len(tair_ids), "genes in region:", min_x, "to", max_x
            tid = tair_ids[0]
            gene_line_list = [tid]
            gene_line_map = {tid:1}
            used_lines = set([1])
            free_lines = set()#[2, 3, 4, 5, 6, 7, 8, 9, 10])
            for next_tid in tair_ids[1:]:
                next_g = gene_dict[next_tid]
                s_pos = next_g['start_pos']
                new_gene_line_list = [next_tid]
                for tid in gene_line_list:
                    g = gene_dict[tid]
                    if next_g['start_pos'] - gene_buffer <= g['end_pos']:
                        new_gene_line_list.append(tid)
                    else:

                        gl = gene_line_map[tid]
                        used_lines.remove(gl)
                        free_lines.add(gl)
                #pdb.set_trace()
                if len(free_lines) > 0:
                    gl = min(free_lines)
                    free_lines.remove(gl)
                elif len(used_lines) > 0:
                    gl = max(used_lines) + 1
                else:
                    gl = 1
                used_lines.add(gl)
                gene_line_map[next_tid] = gl
                #pdb.set_trace()
                gene_line_list = new_gene_line_list


#        if tair_file:
#            with open(tair_file, 'w') as f:
#                for gene in tair_genes:
#                    f.write(str(gene) + "\n")


        displayed_unit = 1000.0 #kbs
        scoreRange = max_score - min_score
        f = plt.figure(figsize=fig_size)
        if plot_genes:
            ax = f.add_axes([0.05, 0.4, 0.94, 0.53])
            tair_ax = f.add_axes([0.05, 0, 0.94, 0.32])
        else:
            ax = f.add_axes([0.05, 0.07, 0.94, 0.84])
        starPoints = [[], [], []]
        color = 'b'
        if chromosome:
            if color_map == None:
                color_map = {1:'b', 2:'g', 3:'r', 4:'c', 5:'m'}
            color = color_map[chromosome]

        chrom = self.snp_results['chromosomes'][0]
        positions = map(lambda x: x / displayed_unit, self.snp_results['positions'])
        scores = self.snp_results['scores']
        for s_i, (score, pos) in enumerate(it.izip(scores, positions)):
            if score > max_score:
                starPoints[0].append(pos)
                starPoints[1].append(max_score)
                starPoints[2].append(score)
                score = max_score
            scores[s_i] = score


        max_pos = max(positions)
        min_pos = min(positions)
        x_range = max_pos - min_pos
        ax.axis([min_pos - 0.02 * x_range, max_pos + 0.02 * x_range,
            min_score - 0.05 * scoreRange, max_score + 0.05 * scoreRange])
        ax.plot(positions, scores, ".", markersize=markersize + 1, alpha=0.7, color=color)

        if tair_ids:
            print "Drawing TAIR genes"
            for tid in tair_ids:
                g = gene_dict[tid]
                y_value = -gene_line_map[tid]
                if len(tair_ids) < 50:
                    tair_ax.text((g['start_pos'] - 0.01 * x_range) / displayed_unit, y_value + 0.2, tid, size=5)
#                if len(g.exons) > 0:
#
#                    for i in g.introns:
#                        tair_ax.plot([i.startPos / displayed_unit, i.endPos / displayed_unit],
#                            [y_value, y_value], color=(0.6, 0.6, 0.6), linewidth=1)
#                    for e in g.exons:
#                        tair_ax.plot([e.startPos / displayed_unit, e.endPos / displayed_unit],
#                            [y_value, y_value], color=(0.3, 0.3, 0.3), linewidth=3)
#                else:
                tair_ax.plot([g['start_pos'] / displayed_unit, g['end_pos'] / displayed_unit],
                        [y_value, y_value], color=(0.3, 0.3, 0.3), linewidth=3)
            tair_ax.spines['top'].set_visible(False)
            tair_ax.spines['bottom'].set_visible(False)
            tair_ax.spines['right'].set_visible(False)
            tair_ax.spines['left'].set_visible(False)
            tair_ax.xaxis.set_visible(False)
            tair_ax.yaxis.set_visible(False)
            min_y = -max([gene_line_map[tid] for tid in gene_line_map])
            max_y = 0
            range_y = max_y - min_y
            tair_ax.axis([ax.get_xlim()[0], ax.get_xlim()[1],
                min_y - 0.05 * range_y, max_y + 0.05 * range_y])



        if cand_genes:
            for cg in cand_genes:
                g = gene_dict[cg.tairID]
                ax.axvspan(g['start_pos'] / displayed_unit, g['end_pos'] / displayed_unit, facecolor='#aa11bb',
                        alpha=0.5, linewidth=0)
                tair_ax.axvspan(g['start_pos'] / displayed_unit, g['end_pos'] / displayed_unit, facecolor='#aa11bb',
                        alpha=0.5, linewidth=0)


#        if highlight_markers:
#            ys = []
#            xs = []
#            for c, p, score in highlight_markers:
#                xs.append(p / displayed_unit)
#                if score > max_score:
#                    ax.text(x, max_score * 1.1, str(round(score, 2)), rotation=45, size="small")
#                    ys.append(max_score)
#                else:
#                    ys.append(score)
#            ax.plot(xs, ys, ".", color="#ff9944", markersize=markersize + 3, alpha=0.7)

        if len(starPoints[0]) > 0:
            ax.plot(starPoints[0], starPoints[1], ".", color="#ee9922", markersize=markersize + 1)
            i = 0
            while i < len(starPoints[0]):
                max_point = i
                cur_pos = starPoints[0][i]
                while i < len(starPoints[0]) and abs(starPoints[0][i] - cur_pos) < 1000000:
                    if starPoints[2][i] > starPoints[2][max_point]:
                        max_point = i
                    i += 1
                ax.text(starPoints[0][max_point] - 200000, (starPoints[1][max_point] - 1) * 1.15, str(round(starPoints[2][max_point], 2)), rotation=45, size="small")




        if plot_bonferroni:
            b_threshold = -math.log10(1.0 / (len(scores) * 20.0))
            if threshold :
                ax.plot(list(ax.get_xlim()), [b_threshold, b_threshold], ":")
                threshold = -math.log10(threshold)
                ax.plot(list(ax.get_xlim()), [threshold, threshold], color='#6495ed', linestyle='-.')
            #Bonferroni threshold
            else:
                ax.plot(list(ax.get_xlim()), [b_threshold, b_threshold], color='#000000', linestyle="-.")

        if not ylab:
            if score_type == "pvals":
                ax.set_ylabel('$ - log(p - $value$)$')

            else:
                ax.set_ylabel('score')
        else:
            ax.set_ylabel(ylab)
        ax.set_xlabel("kilobases")
        ax.set_title('Chromsome %d' % chrom)

        if pdf_file:
            plt.savefig(pdf_file, format="pdf")
        if png_file:
            plt.savefig(png_file, format="png", dpi=300, bbox_inches='tight')
        if not (pdf_file or png_file):
            plt.show()

        plt.clf()
        plt.close()



    def _plot_small_manhattan_2_(self, plot_ax, gene_model_ax=None, min_score=0, max_score=None,
                score_type="pvals", ylab="$-$log$_{10}(p-$value$)$", plot_bonferroni=False,
                cand_genes=None, threshold=0, highlight_markers=None, chromosome=None,
                color_map=None, markersize=6, sign_color=None, title=None,
                xlab='kb'):

        b_threshold = -math.log10(1.0 / (len(self.snp_results['positions']) * 20.0))
        tair_ids = None
        min_x = min(self.snp_results['positions'])
        max_x = max(self.snp_results['positions'])
        if gene_model_ax != None or cand_genes != None:
            gene_dict = get_gene_dict(chromosome, min_x, max_x, only_genes=True)
        if gene_model_ax != None:
            print 'Retrieving genes from DB.'
            gene_buffer = int((max_x - min_x) / 16) #bp
            tair_ids = sorted(gene_dict.keys())
            print "Found", len(tair_ids), "genes in region:", min_x, "to", max_x
            tid = tair_ids[0]
            gene_line_list = [tid]
            gene_line_map = {tid:1}
            used_lines = set([1])
            free_lines = set()#[2, 3, 4, 5, 6, 7, 8, 9, 10])
            for next_tid in tair_ids[1:]:
                next_g = gene_dict[next_tid]
                s_pos = next_g['start_pos']
                new_gene_line_list = [next_tid]
                for tid in gene_line_list:
                    g = gene_dict[tid]
                    if next_g['start_pos'] - gene_buffer <= g['end_pos']:
                        new_gene_line_list.append(tid)
                    else:

                        gl = gene_line_map[tid]
                        used_lines.remove(gl)
                        free_lines.add(gl)
                #pdb.set_trace()
                if len(free_lines) > 0:
                    gl = min(free_lines)
                    free_lines.remove(gl)
                elif len(used_lines) > 0:
                    gl = max(used_lines) + 1
                else:
                    gl = 1
                used_lines.add(gl)
                gene_line_map[next_tid] = gl
                #pdb.set_trace()
                gene_line_list = new_gene_line_list

        if max_score == None:
            max_score = max(self.snp_results['scores'])
            if highlight_markers:
                msh = max([s for c, p, s in highlight_markers])
                max_score = max(max_score, msh)

        displayed_unit = 1000.0 #kbs
        scoreRange = max_score - min_score
        starPoints = [[], [], []]
        color = 'b'
        if chromosome:
            if color_map == None:
                color_map = {1:'b', 2:'g', 3:'r', 4:'c', 5:'m'}
            color = color_map[chromosome]

        chrom = self.snp_results['chromosomes'][0]
        positions = map(lambda x: x / displayed_unit, self.snp_results['positions'])
        scores = self.snp_results['scores']
        sign_points = {'positions':[], 'scores':[]}
        for s_i, (score, pos) in enumerate(it.izip(scores, positions)):
            if score > max_score:
                starPoints[0].append(pos)
                starPoints[1].append(max_score)
                starPoints[2].append(score)
                scores[s_i] = max_score
            elif sign_color != None:
                if score >= b_threshold:
                    sign_points['scores'].append(score)
                    sign_points['positions'].append(pos)

        max_pos = max(positions)
        min_pos = min(positions)
        x_range = max_pos - min_pos
        plot_ax.axis([min_pos - 0.02 * x_range, max_pos + 0.02 * x_range,
            min_score - 0.05 * scoreRange, max_score + 0.05 * scoreRange])
        plot_ax.plot(positions, scores, ".", markersize=markersize + 1, alpha=0.7, color=color, mew=0, mec=color)

        if len(sign_points['positions']) > 0:
            plot_ax.plot(sign_points['positions'], sign_points['scores'], ".", color=sign_color,
                markersize=markersize + 1, mew=0, mec=sign_color)
        if tair_ids:
            print "Drawing TAIR genes"
            for tid in tair_ids:
                g = gene_dict[tid]
                y_value = -gene_line_map[tid]
                if len(tair_ids) < 50:
                    gene_model_ax.text((g['start_pos'] - 0.01 * x_range) / displayed_unit, y_value + 0.2, tid, size=5)
                if tid in cand_genes:
                    gene_model_ax.plot([g['start_pos'] / displayed_unit, g['end_pos'] / displayed_unit],
                            [y_value, y_value], color='#aa11cc', linewidth=4)
                    gene_model_ax.text((g['start_pos'] - 5000) / displayed_unit, y_value - 1.5, 'AtHKT1;1', size=8, color='#9911cc')

                else:
                    gene_model_ax.plot([g['start_pos'] / displayed_unit, g['end_pos'] / displayed_unit],
                            [y_value, y_value], color=(0.3, 0.3, 0.3), linewidth=3)
            gene_model_ax.spines['top'].set_visible(False)
            gene_model_ax.spines['bottom'].set_visible(False)
            gene_model_ax.spines['right'].set_visible(False)
            gene_model_ax.spines['left'].set_visible(False)
            gene_model_ax.xaxis.set_visible(False)
            gene_model_ax.yaxis.set_visible(False)
            min_y = -max([gene_line_map[tid] for tid in gene_line_map])
            max_y = 0
            range_y = max_y - min_y
            gene_model_ax.axis([plot_ax.get_xlim()[0], plot_ax.get_xlim()[1],
                min_y - 0.05 * range_y, max_y + 0.05 * range_y])



        if cand_genes:
            for cg in cand_genes:
                g = gene_dict[cg]
                plot_ax.axvspan(g['start_pos'] / displayed_unit, g['end_pos'] / displayed_unit, facecolor='#9911cc',
                        alpha=0.2, linewidth=0)
#                if gene_model_ax != None:
#                    gene_model_ax.axvspan(g['start_pos'] / displayed_unit, g['end_pos'] / displayed_unit, facecolor='#aa11cc',
#                        alpha=0.5, linewidth=0)


#        if highlight_markers:
#            ys = []
#            xs = []
#            for c, p, score in highlight_markers:
#                xs.append(p / displayed_unit)
#                if score > max_score:
#                    plot_ax.text(x, max_score * 1.1, str(round(score, 2)), rotation=45, size="small")
#                    ys.append(max_score)
#                else:
#                    ys.append(score)
#            plot_ax.plot(xs, ys, ".", color="#ff9944", markersize=markersize + 3, alpha=0.7, mew=0,
#                mec="#ff9944")

        if len(starPoints[0]) > 0:
            plot_ax.plot(starPoints[0], starPoints[1], ".", color="#ee9922", markersize=markersize + 1)
            i = 0
            while i < len(starPoints[0]):
                max_point = i
                cur_pos = starPoints[0][i]
                while i < len(starPoints[0]) and abs(starPoints[0][i] - cur_pos) < 1000000:
                    if starPoints[2][i] > starPoints[2][max_point]:
                        max_point = i
                    i += 1
                plot_ax.text(starPoints[0][max_point] - 200000, (starPoints[1][max_point] - 1) * 1.15, str(round(starPoints[2][max_point], 2)), rotation=45, size="small")




        if plot_bonferroni:
            #b_threshold = -math.log10(1.0 / (len(scores) * 20.0))
            if threshold :
                plot_ax.plot(list(plot_ax.get_xlim()), [b_threshold, b_threshold], ":")
                threshold = -math.log10(threshold)
                plot_ax.plot(list(plot_ax.get_xlim()), [threshold, threshold], color='#6495ed', linestyle='-.')
            #Bonferroni threshold
            else:
                plot_ax.plot(list(plot_ax.get_xlim()), [b_threshold, b_threshold], color='#000000', linestyle="-.")

        if not ylab:
            if score_type == "pvals":
                plot_ax.set_ylabel('$ - log(p - $value$)$')

            else:
                plot_ax.set_ylabel('score')
        else:
            plot_ax.set_ylabel(ylab)
        if xlab:
            plot_ax.set_xlabel(xlab)
        if title == None:
            plot_ax.set_title('Chromsome %d' % chrom)




    def get_rare_haplotype_list(self, sd):
        """
        Assumes SNPs are defined..
        """
        #Get non-included accessions..
        #Find the same SNPs as in this object.
        #Sort by haplotypes.



    def plot_manhattan2(self, ax, plot_type='simple', min_score=None, max_score=None, percentile=50,
            type="pvals", ylab="$-$log$_{10}(p-$value$)$", plot_bonferroni=False, b_threshold=None,
            cand_genes=None, threshold=0, highlight_markers=None, tair_file=None, plot_genes=True,
            plot_xaxis=False, highlight_loci=None, neg_log_transform=False, chrom_colormap=None,
            sign_color=None, markersize=6, wo_xtick_labels=False):
        """
        Generates a manhattan using the given axis (ax)...
        """
        num_scores = len(self.snp_results['scores'])

        "Plotting a Manhattan-style plot with %i markers." % num_scores

        if not b_threshold:
            b_threshold = -math.log10(1.0 / (num_scores * 20.0))

        if not chrom_colormap:
            chrom_colormap = {1:'b', 2:'g', 3:'r', 4:'c', 5:'m'}

        chromosome_ends = self.get_chromosome_ends()
        result = self.simple_clone()
        if neg_log_transform:
            result.neg_log_trans()
        chrom_set = set(result.snp_results['chromosomes'])
        chromosomes = list(chrom_set)
        chromosomes.sort()
        if len(chrom_set) == 1:
            percentile = 0.0
        if percentile != 0.0:
            result.filter_percentile(percentile / 100.0)

        if highlight_markers:
            new_h_markers = []
            if len(highlight_markers[0]) == 2:
                indices = result.get_indices(highlight_markers)
                pvals = [result.snp_results['scores'][i] for i in indices]
                for i, (c, p) in enumerate(highlight_markers):
                    s = -math.log10(pvals[i]) if neg_log_transform else pvals[i]
                    new_h_markers.append((c, p, s))
            elif len(highlight_markers[0]) == 3:
                for c, p, pval in highlight_markers:
                    #s = -math.log10(pval) if neg_log_transform else pval
                    s = pval
                    new_h_markers.append((c, p, s))
            highlight_markers = new_h_markers

        if not max_score:
            max_score = max(result.snp_results['scores'])
            if highlight_markers:
                h_scores = [s for c, p, s in highlight_markers]
                max_score = max(max_score, max(h_scores))
        if not min_score:
            if type == "pvals":
                min_score = 0
            else:
                min_score = min(result.snp_results['scores'])


        if cand_genes:
            #processing candidate genes by chromosome
            chr_cand_genes = {}
            for chrom in chromosomes:
                chr_cand_genes[chrom] = []
            for cg in cand_genes:
                chr_cand_genes[int(cg.chromosome)].append(cg)

        if highlight_loci:
            hl_dict = {}
            for chrom in chromosomes:
                hl_dict[chrom] = []
            for c, p in highlight_loci:
                hl_dict[c].append(p)


        scoreRange = max_score - min_score
        chromosome_splits = result.get_chromosome_splits()
        offset = 0
        ticksList1 = []
        ticksList2 = []
        textPos = []
        starPoints = [[], [], []]
        if sign_color != None:
            sign_points = {'positions':[], 'scores':[]}

        chr_offsets = []
        for i, chromosome_end in enumerate(chromosome_ends):
            chr_offsets.append(offset)
            index1 = chromosome_splits[i]
            index2 = chromosome_splits[i + 1]
            scoreList = result.snp_results['scores'][index1:index2]
            posList = result.snp_results['positions'][index1:index2]
            chrom = chromosomes[i]
            newPosList = [offset + pos for pos in posList]

            for s_i, (score, pos) in enumerate(it.izip(scoreList, newPosList)):
                if score > max_score:
                    starPoints[0].append(pos)
                    starPoints[1].append(max_score)
                    starPoints[2].append(score)
                    scoreList[s_i] = max_score
                if sign_color != None and score >= b_threshold:
                    if highlight_markers != None:
                        for c, p, pval in highlight_markers:
                            if abs(c - chrom) < 1 and abs(p - pos + offset) < 1:
                                break
                        else:
                            sign_points['positions'].append(pos)
                            sign_points['scores'].append(score)
                    else:
                        sign_points['positions'].append(pos)
                        sign_points['scores'].append(score)


            #Marking candidate genes
            if cand_genes:
                for cg in chr_cand_genes[chrom]:
                    ax.axvspan(offset + cg.startPos, offset + cg.endPos,
                            facecolor='#555555', edgecolor='#555555', alpha=0.6)
            #Highlighting loci
            if highlight_loci:
                for p in hl_dict[chrom]:
                    ax.axvspan(offset + p - 50000, offset + p + 50000, ec='#660055', fc='#660055',
                        alpha=0.4,)

            #Plotting scores
            ax.plot(newPosList, scoreList, ".", markersize=markersize, alpha=0.6,
                color=chrom_colormap[i + 1], ls='', mew=0, mec=chrom_colormap[i + 1])

            oldOffset = offset
#            textPos.append(offset + chromosome_end / 2 - 2000000)
            offset += chromosome_end
            if plot_xaxis:
                #if len(chromosome_ends) == 23: #This is probably Human!
                ticksList1.append(oldOffset + chromosome_end / 2)
                if wo_xtick_labels:
                    ticksList2.append('')
                else:
                    if chrom == 23:
                        ticksList2.append('X')
                    elif chrom < 16:
                        ticksList2.append(int(chrom))
                    elif chrom % 2 == 0:
                        ticksList2.append(int(chrom))
                    else:
                        ticksList2.append('')

#                else:
#                    for j in range(oldOffset, offset, 4000000):
#                        ticksList1.append(j)
#                    for j in range(0, chromosome_end, 4000000):
#                        if j % 8000000 == 0 and j < chromosome_end - 4000000 :
#                            ticksList2.append(j / 1000000)
#                        else:
#                            ticksList2.append("")


        if highlight_markers:
            ys = []
            xs = []
            for c, p, score in highlight_markers:
                x = chr_offsets[c - 1] + p
                xs.append(x)
                if score > max_score:
                    ax.text(x, max_score * 1.1, str(round(score, 2)), rotation=45, size="small")
                    ys.append(max_score)
                else:
                    ys.append(score)
            ax.plot(xs, ys, ".", color="#ff9944", markersize=markersize + 2, alpha=0.9, ls='', mew=0,
                mec="#ff9944")


        if len(starPoints[0]) > 0:
            ax.plot(starPoints[0], starPoints[1], ".", color="#ee9922", markersize=markersize + 2, mew=0,
                mec="#ee9922")
            i = 0
            while i < len(starPoints[0]):
                max_point = i
                cur_pos = starPoints[0][i]
                while i < len(starPoints[0]) and abs(starPoints[0][i] - cur_pos) < 3000000:
                    if starPoints[2][i] > starPoints[2][max_point]:
                        max_point = i
                    i += 1
                ax.text(starPoints[0][max_point] - 1000000, (starPoints[1][max_point] - 1) * 1.15,
                    str(round(starPoints[2][max_point], 2)), rotation=45, size="small")

        if sign_color:
            if len(sign_points['positions']) > 0:
                ax.plot(sign_points['positions'], sign_points['scores'], ".", color=sign_color,
                    markersize=markersize + 1, mew=0, mec=sign_color)


        if plot_bonferroni:
            if threshold :
                ax.plot([0, sum(result.chromosome_ends)], [b_threshold, b_threshold], ":")
                threshold = -math.log10(threshold)
                ax.plot([0, sum(result.chromosome_ends)], [threshold, threshold], color='#6495ed',
                    linestyle='-.',)
            #Bonferroni threshold
            else:
                ax.plot([0, sum(result.chromosome_ends)], [b_threshold, b_threshold], color='#000000',
                    linestyle="--")

        if plot_xaxis:
            ax.set_xticks(ticksList1)
            ax.set_xticklabels(ticksList2, fontsize='x-small')
            if not wo_xtick_labels:
                #if len(chromosome_ends) == 23: #This is probably Human!
                ax.set_xlabel("Chromosome")
                #else:
                #    ax.set_xlabel("Mb")
        else:
            ax.set_xlabel("bases")

        x_range = sum(result.chromosome_ends)
        ax.axis([-x_range * 0.01, x_range * 1.01, min_score - 0.05 * scoreRange, max_score + 0.05 * scoreRange])




    def plot_manhattan(self, pdf_file=None, png_file=None, min_score=None, max_score=None, percentile=80,
            type="pvals", ylab="$-$log$_{10}(p-$value$)$", plot_bonferroni=False, b_threshold=None,
            cand_genes=None, threshold=0, highlight_markers=None, tair_file=None, plot_genes=False,
            plot_xaxis=True, highlight_loci=None, neg_log_transform=False, markersize=3,
            chrom_col_map=None):

        """
        Plots a 'Manhattan' style GWAs plot.
        """

        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        num_scores = len(self.snp_results['scores'])

        "Plotting a Manhattan-style plot with %i markers." % num_scores

        chromosome_ends = self.get_chromosome_ends()
        result = self.simple_clone()
        if neg_log_transform:
            result.neg_log_trans()
        chrom_set = set(result.snp_results['chromosomes'])
        chromosomes = list(chrom_set)
        chromosomes.sort()
        if len(chrom_set) == 1:
            percentile = 0.0
        if percentile != 0.0:
            result.filter_percentile(percentile / 100.0)

        if highlight_markers:
            new_h_markers = []
            if len(highlight_markers[0]) == 2:
                indices = result.get_indices(highlight_markers)
                scores = [result.snp_results['scores'][i] for i in indices]
                for i, (c, p) in enumerate(highlight_markers):
                    new_h_markers.append((c, p, scores[i]))
            else:
                new_h_markers = [(cps[0], cps[1], cps[2]) for cps in highlight_markers]
            highlight_markers = new_h_markers

        if not max_score:
            max_score = max(result.snp_results['scores'])
            if highlight_markers:
                h_scores = [s for c, p, s in highlight_markers]
                max_score = max(max_score, max(h_scores))
        if not min_score:
            if type == "pvals":
                min_score = 0
            else:
                min_score = min(result.snp_results['scores'])


        if cand_genes:
            #processing candidate genes by chromosome
            chr_cand_genes = {}
            for chrom in chromosomes:
                chr_cand_genes[chrom] = []
            for cg in cand_genes:
                chr_cand_genes[int(cg.chromosome)].append(cg)

        if highlight_loci:
            hl_dict = {}
            for chrom in chromosomes:
                hl_dict[chrom] = []
            for c, p in highlight_loci:
                hl_dict[c].append(p)

        if len(chrom_set) == 1:
            chrom = chrom_set.pop()
            if cand_genes:
                cand_genes = chr_cand_genes[chrom]
            return result._plot_small_manhattan_(pdf_file=pdf_file, png_file=png_file, min_score=min_score,
                        max_score=max_score, ylab=ylab, plot_bonferroni=plot_bonferroni,
                        cand_genes=cand_genes, threshold=threshold,
                        highlight_markers=highlight_markers, chromosome=chrom,
                        tair_file=tair_file, plot_genes=plot_genes, color_map=chrom_col_map,
                        markersize=markersize)


        scoreRange = max_score - min_score
        offset = 0
        chromosome_splits = result.get_chromosome_splits()
        print chromosome_splits

        ticksList1 = []
        ticksList2 = []
        textPos = []
        plt.figure(figsize=(11, 2.8))
        plt.axes([0.045, 0.15, 0.95, 0.71])
        starPoints = [[], [], []]
        chr_offsets = []
        for i, chromosome_end in enumerate(chromosome_ends):
            chr_offsets.append(offset)
            print i, chromosome_splits
            index1 = chromosome_splits[i]
            index2 = chromosome_splits[i + 1]
            scoreList = result.snp_results['scores'][index1:index2]
            posList = result.snp_results['positions'][index1:index2]
            chrom = chromosomes[i]
            newPosList = [offset + pos for pos in posList]

            for s_i, (score, pos) in enumerate(it.izip(scoreList, newPosList)):
                if score > max_score:
                    starPoints[0].append(pos)
                    starPoints[1].append(max_score)
                    starPoints[2].append(score)
                    score = max_score
                scoreList[s_i] = score

            if not chrom_col_map:
                plt.plot(newPosList, scoreList, ".", markersize=markersize, alpha=0.7, mew=0)
            else:
                color = chrom_col_map[chrom]
                plt.plot(newPosList, scoreList, ".", markersize=markersize, alpha=0.7, color=color, mew=0)

            if cand_genes:
                for cg in chr_cand_genes[chrom]:
                    plt.axvspan(offset + cg.startPos, offset + cg.endPos,
                            facecolor='#FF9900', alpha=0.6)

            if highlight_loci:
                for p in hl_dict[chrom]:
                    plt.axvline(offset + p, color='#1166FF', alpha=0.6)

            oldOffset = offset
#            textPos.append(offset + chromosome_end / 2 - 2000000)
            offset += chromosome_end
            if plot_xaxis:
                if len(chromosome_ends) == 23: #This is probably Human!
                    ticksList1.append(oldOffset + chromosome_end / 2)
                    if chrom == 23:
                        ticksList2.append('X')
                    else:
                        ticksList2.append(chrom)
                else:
                    for j in range(oldOffset, offset, 8000000):
                        ticksList1.append(j)
                    for j in range(0, chromosome_end, 8000000):
                        if j % 16000000 == 0 and j < chromosome_end - 8000000 :
                            ticksList2.append(j / 1000000)
                        else:
                            ticksList2.append("")


        plt.plot(starPoints[0], starPoints[1], ".", color="#ee9922", markersize=markersize + 2, mew=0)
        if len(starPoints[0]) > 0:
            i = 0
            while i < len(starPoints[0]):
                max_point = i
                cur_pos = starPoints[0][i]
                while i < len(starPoints[0]) and abs(starPoints[0][i] - cur_pos) < 3000000:
                    if starPoints[2][i] > starPoints[2][max_point]:
                        max_point = i
                    i += 1
                plt.text(starPoints[0][max_point] - 1000000, (starPoints[1][max_point] - 1) * 1.15, str(round(starPoints[2][max_point], 2)), rotation=45, size="small")


        if highlight_markers:
            ys = []
            xs = []
            for c, p, score in highlight_markers:
                x = chr_offsets[c - 1] + p
                xs.append(x)
                if score > max_score:
                    plt.text(x, max_score * 1.1, str(round(score, 2)), rotation=45, size="small")
                    ys.append(max_score)
                else:
                    ys.append(score)
            plt.plot(xs, ys, ".", color="#ff9944", markersize=markersize + 4, alpha=0.8, mew=0)


        if plot_bonferroni:
            if not b_threshold:
                b_threshold = -math.log10(1.0 / (num_scores * 20.0))
            if threshold :
                plt.plot([0, sum(result.chromosome_ends)], [b_threshold, b_threshold], ":")
                threshold = -math.log10(threshold)
                plt.plot([0, sum(result.chromosome_ends)], [threshold, threshold], color='#6495ed', linestyle='-.')
            #Bonferroni threshold
            else:
                plt.plot([0, sum(result.chromosome_ends)], [b_threshold, b_threshold], color='#000000', linestyle="-.")

        x_range = sum(result.chromosome_ends)
        plt.axis([-x_range * 0.01, x_range * 1.01, min_score - 0.05 * scoreRange, max_score + 0.05 * scoreRange])
        if plot_xaxis:
            plt.xticks(ticksList1, ticksList2, fontsize='x-small')
        #pdb.set_trace()
        if not ylab:
            if type == "pvals":
                plt.ylabel('$ - log(p - $value$)$')

            else:
                plt.ylabel('score')
        else:
            plt.ylabel(ylab)
        if plot_xaxis:
            if len(chromosome_ends) == 23: #This is probably Human!
                plt.xlabel("Chromosome")
            else:
                plt.xlabel("Mb")
        else:
            plt.xlabel("bases")


        if pdf_file:
            plt.savefig(pdf_file, format="pdf")
        if png_file:
            plt.savefig(png_file, format="png", dpi=300, bbox_inches='tight')
        if not (pdf_file or png_file):
            plt.show()

        plt.clf()
        plt.close()


    def plot_qq(self, file_prefix, exp_scores=None):
        from mixmogam import analyze_gwas_results as agr
        ks_stat = agr.calc_ks_stats(self.snp_results['scores'], exp_scores=exp_scores)
        exp_median = 0.5
        if exp_scores:
            exp_median = exp_scores[len(exp_scores) / 2]
        p_med = agr.calc_median(self.snp_results['scores'], exp_median=exp_median)
        stat_text = 'D=%0.8f,  M=%0.8f' % (ks_stat['D'], p_med)
        quantiles = agr.get_quantiles(self.snp_results['scores'])
        agr.simple_qqplot([quantiles], file_prefix + '.png', text=stat_text)
        log_quantiles = agr.get_log_quantiles(self.snp_results['scores'])
        agr.simple_log_qqplot([log_quantiles], file_prefix + '_log.png', text=stat_text)



    def get_chromosome_splits(self):
        """
        Returns list of indices (and prev chromosome), for the when the chromosomes
        change in the scores, positions indices.
        """
        last_chrom = -1
        chromosome_splits = []
        for i, chrom in enumerate(self.snp_results['chromosomes']):
            if last_chrom != chrom:
                while last_chrom < chrom:
                    last_chrom += 1
                    chromosome_splits.append(i)
        chromosome_splits.append(i)
        return chromosome_splits


    def neg_log_trans(self):
        """
        Apply - log(x) to the pvalues (scores)
        """
        import math
        f = lambda x: 323.0 if x == 0.0 else -math.log10(x)
        self.snp_results['scores'] = map(f, self.snp_results['scores'])


    def filter_percentile(self, percentile, reversed=False):
        """
        Filter above (or below) percentile.
        """
        sl = self.snp_results['scores'][:]
        sl.sort()
        score_cutoff = sl[int(len(sl) * percentile)]
        self.filter_attr("scores", score_cutoff, reversed=reversed)


    def filter_common_top_scores(self, result, top_fraction=0.1):
        import bisect
        self._rank_scores_()
        result._rank_scores_()
        keep_indices_1 = set()
        keep_indices_2 = set()
        chr_pos_list_1 = self.get_chr_pos_list()
        chr_pos_list_2 = result.get_chr_pos_list()
        for i in self.orders[:int(top_fraction * len(self.orders))]:
            j = bisect.bisect(chr_pos_list_2, chr_pos_list_1[i]) - 1
            if chr_pos_list_2[j] == chr_pos_list_1[i]:
                keep_indices_1.add(i)

        for i in result.orders[:int(top_fraction * len(result.orders))]:
            j = bisect.bisect(chr_pos_list_1, chr_pos_list_2[i]) - 1
            if chr_pos_list_1[j] == chr_pos_list_2[i]:
                keep_indices_2.add(i)

        keep_indices_list_1 = list(keep_indices_1)
        for i in keep_indices_2:
            keep_indices_1.add(bisect.bisect(chr_pos_list_1, chr_pos_list_2[i]) - 1)

        for i in keep_indices_list_1:
            keep_indices_2.add(bisect.bisect(chr_pos_list_2, chr_pos_list_1[i]) - 1)

        print len(keep_indices_1), len(keep_indices_2)
        self.filter_indices(list(keep_indices_1))
        result.filter_indices(list(keep_indices_2))
        return keep_indices_1, keep_indices_2


    def filter_indices(self, indices_to_keep):
        indices_to_keep.sort()
        snp_results = {}
        for info in self.snp_results:
            snp_results[info] = []
        count = len(self.scores)
        for i in indices_to_keep:
            for info in self.snp_results:
                if self.snp_results[info]:
                    snp_results[info].append(self.snp_results[info][i])

        self.snp_results = snp_results
        print "%i scores were removed." % (count - len(self.scores))


    def filter_attr(self, attr_name, attr_threshold, reversed=False, verbose=True, return_clone=False):
        """
        Filter out scores / pvalues etc. which have attr < attr_threshold.

        attr are e.g.
        'mafs', 'macs', 'scores', etc.
        """
        if verbose:
            print "Filtering for attribute ' % s' with threshold: %g" % (attr_name, attr_threshold)
        attr = self.snp_results[attr_name]
        snp_results = {}
        for info in self.snp_results:
            snp_results[info] = []
        count = len(self.snp_results['scores'])
        for i in range(count):
            if reversed:
                if attr[i] <= attr_threshold:
                    for info in self.snp_results:
                        if len(self.snp_results[info]) > 0:
                            snp_results[info].append(self.snp_results[info][i])
            else:
                if attr[i] >= attr_threshold:
                    for info in self.snp_results:
                        if len(self.snp_results[info]) > 0:
                            snp_results[info].append(self.snp_results[info][i])

        num_scores = len(snp_results['scores'])
        if verbose:
            print "%i scores were removed out of %i." % (count - num_scores, count)

        if return_clone:
            return Result(snp_results=snp_results, accessions=self.accessions, name=self.name,
                         result_type=self.result_type, phen_id=self.phen_id, keys=self.keys)
        else:
            self.snp_results = snp_results
            return num_scores



    def get_min_distances(self, chrom_pos_list):
        """
        Return distance statistics on this result object in relation to the chromosome position list.
        """
        cpl = self.get_chr_pos_list()
        cp_array = sp.array(chrom_pos_list)
        cpl_dists = map(sp.absolute, [sp.array(cpt) - cp_array for cpt in cpl])
        min_dists = []
        for cp_dists in cpl_dists:
            curr_min = -1 #different chromosome.
            for cp_dist in cp_dists:
                if cp_dist[0] == 0:
                    if curr_min == -1:
                        curr_min = cp_dist[1]
                    else:
                        curr_min = min(curr_min, cp_dist[1])
            min_dists.append(curr_min)
        return min_dists


    def get_distances(self, chrom_pos_list):
        """
        Return a list of lists of distances.
        """
        cpl = self.get_chr_pos_list()
        cp_array = sp.array(chrom_pos_list)
        cpl_dists = map(sp.absolute, [sp.array(cpt) - cp_array for cpt in cpl])
        return cpl_dists


    def get_gene_analysis(self, gene, bin_radius=[0, 1000, 5000, 10000, 25000, 50000, 100000, 'same_chrom', 'all'], h5py_group=None):
        """
        Returns various statistics on a gene...
         - Distance from min pval
         -  
        """

        min_i = self.arg_min_attr()
        min_chr_pos = sp.array([self.snp_results['chromosomes'][min_i], self.snp_results['positions'][min_i]])
        gene_chr_pos = sp.array([[gene.chromosome, gene.startPos], [gene.chromosome, gene.endPos]])
        d = {'dist_to_min_pval': sp.absolute(min_chr_pos - gene_chr_pos).min(0)}

        bin_dict = {}
        gene.chromosome
        chrom_pos_list = self.get_chr_pos_list()

        for r in bin_radius:
            if r == 'all':
                scores = sp.array(self.snp_results['scores'])
            elif r == 'same_chrom':
                start_i = bisect.bisect(chrom_pos_list, (gene.chromosome, 0))
                stop_i = bisect.bisect(chrom_pos_list, (gene.chromosome + 1, 0))
                scores = sp.array(self.snp_results['scores'][start_i:stop_i])
            else:
                start_pos = gene.startPos - r
                end_pos = gene.endPos + r
                start_i = bisect.bisect(chrom_pos_list, (gene.chromosome, start_pos))
                stop_i = bisect.bisect(chrom_pos_list, (gene.chromosome, end_pos))
                scores = sp.array(self.snp_results['scores'][start_i:stop_i])
            if len(scores) > 0:
                bin_dict[r] = {'min_pval':scores.min(), 'num_snps':len(scores)}
            else:
                bin_dict[r] = {'min_pval':1, 'num_snps':len(scores)}
        d['bin_dict'] = bin_dict
        if h5py_group != None:
            h5py_group.create_dataset('dist_to_min_pval', data=d['dist_to_min_pval'])
            bin_stats = h5py_group.create_group('bin_dict')
            for r in bin_radius:
                brg = bin_stats.create_group(str(r))
                brg.create_dataset('min_pval', data=bin_dict[r]['min_pval'])
                brg.create_dataset('num_snps', data=bin_dict[r]['num_snps'])
        return d



    def get_farthest_w_stronger_association(self, chrom_pos_list, assume_ranked=False):
        """
        Return distance to farthest locus with greater association than some of the given loci.
        (chr_dist,pos_dist)
        """
        if not assume_ranked:
            self._rank_scores_()
        c_indices = self.get_indices(chrom_pos_list)
        c_ranks = [self.ranks[i] for i in c_indices]
        #max_rank_i = sp.argmax(c_ranks)
        stronger_indices = self.orders[:max(c_ranks) - 1]
        if len(stronger_indices) == 0:
            return [0, 0] #chrom,pos 
        cpl = self.get_chr_pos_list()
        stronger_cpl = sp.array([cpl[i] for i in stronger_indices])
        cp = chrom_pos_list[0]
        dist_l = sp.absolute(stronger_cpl - sp.array(cp)).tolist()
        for cp in chrom_pos_list[1:]:
            dist_list = sp.absolute(stronger_cpl - sp.array(cp)).tolist()
            dist_l = [min(d1, d2) for d1, d2 in it.izip(dist_l, dist_list)]
        max_dist = max(dist_l)
        return max_dist










#    def filter_non_segregating_snps(self, ecotype1, ecotype2, accessions=None):
#        """
#        Filter out all SNPs which are not segregating in the two accessions.
#
#        Assumes the accessions map the results objects SNPs. (and that they are defined)
#        """
#        newScores = []
#        newPositions = []
#        newChromosomes = []
#        newMafs = []
#        newMarfs = []
#        new_snps = []
#
#        if accessions:
#            ecotypes = accessions
#        else:
#            ecotypes = self.accessions
#
#        e_i1 = ecotypes.index(ecotype1)
#        e_i2 = ecotypes.index(ecotype2)
#
#        for i in range(len(self.snps)):
#            snp = self.snps[i]
#            if snp[e_i1] != snp[e_i2]:
#                newScores.append(self.scores[i])
#                newPositions.append(self.positions[i])
#                newChromosomes.append(self.chromosomes[i])
#                newMafs.append(self.mafs[i])
#                newMarfs.append(self.marfs[i])
#                new_snps.append(snp)
#
#        self.scores = newScores
#        self.positions = newPositions
#        self.chromosomes = newChromosomes
#        self.mafs = newMafs
#        self.marfs = newMarfs
#        self.snps = new_snps

    def get_scores(self, cloned=True):
        if cloned:
            return self.snp_results['scores'][:]
        else:
            return self.snp_results['scores']

    def get_top_snps_result(self, n):
        """
        returns top n SNPs
        """
        import copy
        result = copy.deepcopy(self) #Cloning
        result.filter_top_snps(n)
        return result


    def get_top_snps(self, n, reversed=False):
        """
        returns top n SNPs
        """
        self._rank_scores_() #Making sure the ranks are updated
        chr_pos_list = []
        if reversed:
            for i in self.orders[-n:]:
                chr_pos_list.append((self.snp_results['chromosomes'][i], self.snp_results['positions'][i]))
        else:
            for i in self.orders[:n]:
                chr_pos_list.append((self.snp_results['chromosomes'][i], self.snp_results['positions'][i]))

        return chr_pos_list


    def get_indices(self, chr_pos_list):
        """
        Return the indices of the given locus.
        """
        import bisect
        indices = []
        for chromosome, position in chr_pos_list:
            cpl = self.get_chr_pos_list()
            i = bisect.bisect(cpl, (chromosome, position))
            if cpl[i - 1] == (chromosome, position):
                indices.append(i - 1)
            else:
                indices.append(-1)
        return indices


#    def get_top_genes(self, n, window_size=5000, conn=None):
#        """
#        Returns a set of (chromosome, start_pos, end_pos), for genes found.
#        """
#        self._rank_scores_() #Making sure the ranks are updated
#        genes = set()
#        snp_ix = []
#        i = 0
#        while len(genes) < n:
#            snp_i = self.orders[i]
#            p = self.snp_results['positions'][snp_i]
#            c = self.snp_results['chromosomes'][snp_i]
#            c_genes = get_gene_list(p - window_size, p + window_size, c, include_intron_exons=False, conn=conn)
#            for g in c_genes:
#                genes.add((g.chromosome, g.startPos, g.endPos, g.tairID))
#            #print len(genes)
#            i += 1
#        print 'found % d genes' % len(genes)
#        return genes


    def get_chr_pos_score_list(self):
        return zip(self.snp_results['chromosomes'], self.snp_results['positions'], self.snp_results['scores'])


    def get_chrom_score_pos_dict(self):
        d = {}
        for chrom in [1, 2, 3, 4, 5]:
            d[chrom] = {'scores':[], 'positions':[], 'macs':[]}
        iterator = it.izip(self.snp_results['chromosomes'], self.snp_results['positions'], \
                self.snp_results['scores'], self.snp_results['macs'])#, self.snp_results['perc_var'])
        for chrom, pos, score, mac in iterator:
            d[chrom]['scores'].append(score)
            d[chrom]['macs'].append(mac)
            d[chrom]['positions'].append(pos)
        return d


    def get_chr_pos_list(self):
        return zip(self.snp_results['chromosomes'], self.snp_results['positions'])


    def get_top_regions(self, n, distance_threshold=25000):
        """
        Returns a list of regions, defined by (chromosome, start_pos, end_pos).
        """
        self._rank_scores_()
        chromosome_ends = self.get_chromosome_ends()

        i = 0 #how many SNPs are needed
        region_count = 0
        regions = [[], [], [], [], []]  #a list of (chromosome,start,stop)
        while sum(map(len, regions)) < n:
            snp_i = self.orders[i]
            snp_pos = self.snp_results['positions'][snp_i]
            chromosome = self.snp_results['chromosomes'][snp_i]
            new_snp_reg = (chromosome, max(snp_pos - distance_threshold, 0), min(snp_pos + distance_threshold, chromosome_ends[chromosome - 1]))
            covered = False
            for reg_i, reg in enumerate(regions[chromosome - 1]):
                if (new_snp_reg[1] < reg[2] and new_snp_reg[2] > reg[2]) or (new_snp_reg[1] < reg[1] and new_snp_reg[2] > reg[1]): #then merge                    
                    regions[chromosome - 1].pop(reg_i) #Removing the region which we're merging with.
                    new_snp_reg = (reg[0], min(reg[1], new_snp_reg[1]), max(reg[2], new_snp_reg[2]))
                elif (new_snp_reg[1] >= reg[1] and new_snp_reg[2] <= reg[2]):
                    covered = True  #It's already covered
                    break
            if not covered:
                regions[chromosome - 1].append(new_snp_reg)
            i += 1
        print "It took %i SNPS to find %i regions." % (i, n)
        regions_list = regions[0] + regions[1] + regions[2] + regions[3] + regions[4]
        return regions_list


    def get_regions(self, gene_radius=20000):
        """
        Returns a list of [[chr,start],[chr,end]] elements, 
        """
        regions = []
        num_scores = len(self.snp_results['scores'])
        iter = enumerate(it.izip(self.snp_results['chromosomes'], self.snp_results['positions']))
        th = sp.array([0, gene_radius])
        (pos_i, t) = iter.next()
        chr_pos = sp.array(t)
        while pos_i < num_scores - 1:
            start_chr_pos = chr_pos
            end_chr_pos = start_chr_pos
            while pos_i < num_scores - 1 and sp.all((chr_pos - end_chr_pos) <= th):
                end_chr_pos = chr_pos
                (pos_i, t) = iter.next()
                chr_pos = sp.array(t)
            if tuple(end_chr_pos) < tuple(start_chr_pos):
                pdb.set_trace()
            if pos_i == num_scores - 1: #Last one..
                if sp.all((chr_pos - end_chr_pos) <= th):
                    regions.append([start_chr_pos - th, chr_pos + th])
                else:
                    regions.append([start_chr_pos - th, end_chr_pos + th])
                    regions.append([chr_pos - th, chr_pos + th])

            else:
                regions.append([start_chr_pos - th, end_chr_pos + th])
        start_chr_pos = chr_pos
        end_chr_pos = start_chr_pos
        regions.append([start_chr_pos - th, end_chr_pos + th])
        print 'Found % d regions' % len(regions)
        return regions


    def count_nearby(self, chrom_pos_list, radius):
        """
        Counts how many of the chrom_pos tuples are within a radius of a SNP in the result.
        """
        cpl = self.get_chr_pos_list()
        count = 0
        for (chrom1, pos1) in chrom_pos_list:
            for (chrom2, pos2) in cpl:
                if chrom1 == chrom2 and abs(pos1 - pos2) <= radius:
                    count += 1
                    break
        return count


    def get_power_analysis(self, caus_chrom_pos_list, window_sizes=[0], debug=False):
        """
        Calculate Power and FDR..
        For all and for each causal individually
        """
        tprs = []
        fdrs = []
        caus_founds_list = []
        for window_size in window_sizes:
            cpl = self.get_chr_pos_list()
            if window_size == 0 and debug:
                pdb.set_trace()
            if len(cpl):
                filtered_cpl = set(cpl)
                caus_founds = []
                for (chrom1, pos1) in caus_chrom_pos_list:
                    caus_found = 0.0
                    for chrom2, pos2 in cpl:
                        if chrom1 == chrom2 and abs(pos1 - pos2) <= window_size:
                            caus_found = 1.0
                            try:
                                filtered_cpl.remove((chrom2, pos2))
                            except Exception:
                                print '(%d,%d) has alreay been filtered' % (chrom2, pos2)
                    caus_founds.append(caus_found)
                tprs.append(float(sum(caus_founds)) / len(caus_chrom_pos_list))
                fdrs.append(float(len(filtered_cpl)) / len(cpl))
                caus_founds_list.append(caus_founds)
            else:
                tprs.append(0.0)
                fdrs.append(0.0)
                caus_founds_list.append([0.0 for cp in caus_chrom_pos_list])
            if window_size == 0 and debug:
                pdb.set_trace()

        return {'tprs':tprs, 'fdrs':fdrs, 'caus_founds':caus_founds_list}



    def get_region_result(self, chromosome, start_pos, end_pos, buffer=0):
        """
        returns a result object with only the SNPs, etc. within the given boundary.
        """
        snp_results = {}
        for k in self.keys:
            snp_results[k] = []
        i = 0
        start_pos -= buffer
        end_pos += buffer

        chromosomes = self.snp_results['chromosomes']
        positions = self.snp_results['positions']
        while i < len(chromosomes) and chromosomes[i] != chromosome:
            i += 1
        while i < len(chromosomes) and positions[i] < start_pos:
            i += 1
        if i == len(chromosomes):
            raise Exception("region results never found!")

        while i < len(chromosomes) and positions[i] <= end_pos and chromosomes[i] == chromosome:
            for k in self.keys:
                snp_results[k].append(self.snp_results[k][i])
            i += 1

        return Result(result_type=self.result_type, snp_results=snp_results, phen_id=self.phen_id,
            accessions=self.accessions, keys=self.keys[:])



    def get_top_region_results(self, n, distance_threshold=25000, buffer=25000):
        reg_results = []
        i = 0
        for reg in self.get_top_regions(n, distance_threshold):
            reg_results.append(self.get_region_result(*reg, buffer=buffer))
        print "Regions were retrieved."
        return reg_results




    def get_chromosome_ends(self):
        if not self.chromosome_ends:
            self._sort_by_chr_pos_()
            i = 0
            chromosome_ends = []
            chromosomes = self.snp_results['chromosomes']
            while i < len(chromosomes):
                curr_chr = chromosomes[i]
                while i < len(chromosomes) and chromosomes[i] == curr_chr:
                    i += 1
                chromosome_ends.append(int(self.snp_results['positions'][i - 1]))
            self.chromosome_ends = chromosome_ends
        return self.chromosome_ends



    def clone(self):
        import copy
        result = copy.deepcopy(self) #Cloning
        return result


    def simple_clone(self):
        snp_results = {}
        for k in self.keys:
            print k
            snp_results[k] = self.snp_results[k][:]
        if self.accessions:
            accessions = self.accessions[:]
        else:
            accessions = None
        result = Result(snp_results=snp_results, accessions=accessions, keys=self.keys[:])
        result.chromosome_ends = self.chromosome_ends[:]
        return result


#
#
#    def na_mafs(self, min_maf=10):
#        """
#        NA scores/pvalues which have maf<minMaf.        
#        """
#        for i in range(0, len(self.scores)):
#            if self.mafs[i] < min_maf:
#                self.scores[i] = "NA"


    def get_num_scores(self):
        return len(self.snp_results['scores'])


    def get_max_snp(self):
        max_val = max(self.snp_results['scores'])
        mi = self.snp_results['scores'].index(max_val)
        return (self.snp_results['snps'][mi], self.snp_results['scores'][mi],
            self.snp_results['chromosomes'][mi], self.snp_results['positions'][mi])


    def min_score(self):
        return min(self.snp_results['scores'])

    def max_score(self):
        return max(self.snp_results['scores'])


    def arg_min_attr(self, attr_name='scores'):
        """
        Returns the index of the minimum attribute.
        """
        return sp.argmin(self.snp_results[attr_name])

    def arg_max_attr(self, attr_name='scores'):
        """
        Returns the index of the maximum attribute.
        """
        return sp.argmax(self.snp_results[attr_name])


    def num_scores(self):
        return len(self.snp_results['scores'])

    def write_to_file(self, filename, additional_columns=None, auto_pickling_on=False, only_pickled=False,
            max_fraction=1.0, columns=None):
        if columns == None:
            columns = ['chromosomes', 'positions', 'scores', 'mafs', 'macs']
        if additional_columns:
            for info in additional_columns:
                if info in self.snp_results:
                    columns.append(info)
        score_threshold = None
        if max_fraction < 1:
            score_threshold = self.get_quantile(max_fraction)
        num_scores = self.num_scores()
#        try:
        if not only_pickled:
            print 'Writing to file: %s'%filename
            with open(filename, "w") as f:
                f.write(','.join(columns) + "\n")
                if score_threshold:
                    for i in range(num_scores):
                        if self.snp_results['scores'][i] < score_threshold:
                            l = [self.snp_results[c][i] for c in columns]
                            l = map(str, l)
                            f.write(",".join(l) + "\n")
                else:
                    for i in range(num_scores):
                        l = [self.snp_results[c][i] for c in columns]
                        l = map(str, l)
                        f.write(",".join(l) + "\n")
        if auto_pickling_on or only_pickled:
            pickle_file = filename + '.pickled'
            d = {}
            for attr in self.pickle_attributes:
                if attr == 'snp_results':
                    snp_results = getattr(self, attr)
                    srd = {}
                    for k in snp_results:
                        if k != 'snps': #Avoid storing SNPs
                            srd[k] = snp_results[k]
                    d[attr] = srd
                    print srd.keys()
                else:
                    d[attr] = getattr(self, attr)
            with open(pickle_file, 'wb') as f:
                cPickle.dump(d, f)






def get_gene_dict(chrom=None, start_pos=None, end_pos=None, only_genes=False):
    from mixmogam import dataParsers as dp
    return dp.parse_tair_gff_file(chrom, start_pos, end_pos, only_genes)





def load_cand_genes_file(file_name, format=1):
    """
    Loads a candidate gene list from a csv file...
    """
    f = open(file_name, "r")
    #print(f.readline())
    reader = csv.reader(f)
    tair_ids = []
    gene_names = []
    if format == 1:
        for row in reader:
            tair_ids.append(row[0].upper())
            if len(row) > 1:
                gene_names.append(row[1])
    f.close()
    return tair_ids





if __name__ == "__main__":
    pass
