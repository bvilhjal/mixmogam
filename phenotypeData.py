"""
A class and functions useful for handling phenotype data.

Author: Bjarni J. Vilhjalmsson
Email: bjarni.vilhjalmsson@gmail.com
"""
import itertools as it
import sys
try:
    import scipy as sp
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception, err_str:
    print 'scipy/matplotlib is missing:', err_str
import warnings


class phenotype_data:
    """
    A class that encapsulates phenotype values and provides basic functionality for these.
    
    This is an update of an older class.
    """

    def __init__(self, phen_dict=None, phenotype_names=[], phen_ids=None):
        if phen_dict:
            self.phen_dict = phen_dict
            self.phen_ids = phen_dict.keys()
            for pid in phen_dict:
                self.phen_dict[pid]['transformation'] = None
                self.phen_dict[pid]['raw_values'] = []
        else:
            if phen_ids:
                self.phen_ids = phen_ids
            else:
                self.phen_ids = range(len(phenotype_names))
            self.phen_dict = {}
            for i, pid in enumerate(self.phen_ids):
                self.phen_dict[pid] = {'name':phenotype_names[i], 'ecotypes':[], 'values':[], 'transformation':None, 'raw_values':[]}


    def num_traits(self):
        return len(self.phen_dict)

    def num_vals(self, pid):
        return len(self.phen_dict[pid]['values'])

    def get_pseudo_heritabilities(self, K, pids=None):
        phers = []
        pvals = []
        if pids == None:
            pids = self.phen_ids
        for pid in pids:
            d = self.get_pseudo_heritability(pid, K)
            phers.append(d['pseudo_heritability'])
            pvals.append(d['pval'])
        return {'phers':phers, 'pvals':pvals}



    def get_blups(self, K, pids=None):
        """
        Returns the BLUP for all traits, along with the pseudo-heritability.
        """
        phers = []
        pvals = []
        blups = []
        blup_residuals = []
        if pids == None:
            pids = self.phen_ids
        for pid in pids:
            d = self.get_blup(pid, K)
            phers.append(d['pseudo_heritability'])
            pvals.append(d['pval'])
            blups.append(d['u_blup'])
            blup_residuals.append(d['blup_residuals'])
        return {'phers':phers, 'pvals':pvals, 'blups':blups, 'blup_residuals':blup_residuals}



    def get_pseudo_heritability(self, pid, K):
        """
        Returns the REML estimate of the heritability.
        
        methods: 'avg' (averages), 'repl' (replicates)
        """
        from scipy import stats
        import linear_models as lm
        phen_vals = self.get_values(pid)
        lmm = lm.LinearMixedModel(phen_vals)
        if len(phen_vals) == len(set(phen_vals)):
            lmm.add_random_effect(K)
        else:
            Z = self.get_incidence_matrix(pid)
            lmm.add_random_effect(Z * K * Z.T)
        r1 = lmm.get_REML()
        ll1 = r1['max_ll']
        rlm = lm.LinearModel(phen_vals)
        ll0 = rlm.get_ll()
        lrt_stat = 2 * (ll1 - ll0)
        pval = stats.chi2.sf(lrt_stat, 1)
        return {'pseudo_heritability':r1['pseudo_heritability'], 'pval':pval}



    def get_blup(self, pid, K):
        """
        Returns the REML estimate for the BLUP and the pseudo-heritability.
        """
        from scipy import stats
        import linear_models as lm
        phen_vals = self.get_values(pid)
        lmm = lm.LinearMixedModel(phen_vals)
        if len(phen_vals) == len(set(phen_vals)):
            lmm.add_random_effect(K)
        else:
            Z = self.get_incidence_matrix(pid)
            lmm.add_random_effect(Z * K * Z.T)
        r1 = lmm.get_REML()
        ll1 = r1['max_ll']
        rlm = lm.LinearModel(phen_vals)
        ll0 = rlm.get_ll()
        lrt_stat = 2 * (ll1 - ll0)
        pval = stats.chi2.sf(lrt_stat, 1)

        #Now the BLUP.
        y_mean = sp.mean(lmm.Y)
        Y = lmm.Y - y_mean
        p_herit = r1['pseudo_heritability']
        delta = (1 - p_herit) / p_herit
#        if K_inverse == None:
#            K_inverse = K.I
#        M = (sp.eye(K.shape[0]) + delta * K_inverse)
#        u_blup = M.I * Y
        M = (K + delta * sp.eye(K.shape[0]))
        u_mean_pred = K * (M.I * Y)
        blup_residuals = Y - u_mean_pred
        return {'pseudo_heritability':r1['pseudo_heritability'], 'pval':pval, 'u_blup':u_mean_pred, 'blup_residuals':blup_residuals}



    def get_broad_sense_heritability(self, pids=None):
        """
        Estimates the broad sense heritability from replicates.
        """
        import linear_models as lm
        if not pids:
            pids = sorted(self.phen_dict.keys())
        bs_herits = []
        bs_avg_herits = []
        bs_pids = []
        bs_herit_pvals = []
        for pid in pids:
            ets = map(int, self.get_ecotypes(pid))
            num_vals = len(ets)
            num_ets = len(set(ets))
            if num_vals == num_ets:
                print "Can't estimate the broad sense heritability when replicates are missing."
                continue
            else:
                avg_repl_num = float(num_vals) / num_ets
                print 'Average number of replicates:', avg_repl_num
            values = self.get_values(pid)
            mod = lm.LinearModel(values)
            res = mod.anova_f_test([sp.array(ets)])
            bs_herit = res['var_perc'][0]
            bs_herit_pval = res['ps'][0]
            bs_herit_pvals.append(bs_herit_pval)
            print 'Heritability:', bs_herit
            print 'Heritability (different from 0) p-value :', bs_herit_pval
            bs_avg_herit = 1.0 - (1.0 - bs_herit) / avg_repl_num
            bs_avg_herits.append(bs_avg_herit)
            print 'Estimated mean value heritability:', bs_avg_herit
            bs_herits.append(bs_herit)
            bs_pids.append(pid)
        return {'bs_herits':bs_herits, 'bs_pids':bs_pids, 'bs_avg_herits':bs_avg_herits, 'bs_herit_pvals':bs_herit_pvals}


    def log_transform(self, pid, method='standard'):
        a = sp.array(self.phen_dict[pid]['values'])
        if method == 'standard':
            vals = sp.log((a - min(a)) + 0.1 * sp.std(a))
        else:
            vals = sp.log(a)
        if not self.phen_dict[pid]['transformation']:
            self.phen_dict[pid]['raw_values'] = self.phen_dict[pid]['values']
            self.phen_dict[pid]['transformation'] = 'log'
        else:
            self.phen_dict[pid]['transformation'] = 'log(' + self.phen_dict[pid]['transformation'] + ')'
        self.phen_dict[pid]['values'] = vals.tolist()
        return True

    def sqrt_transform(self, pid, method='standard'):
        a = sp.array(self.phen_dict[pid]['values'])
        if method == 'standard':
            vals = sp.sqrt((a - min(a)) + 0.1 * sp.std(a))
        else:
            vals = sp.sqrt(a)
        if not self.phen_dict[pid]['transformation']:
            self.phen_dict[pid]['raw_values'] = self.phen_dict[pid]['values']
            self.phen_dict[pid]['transformation'] = 'sqrt'
        else:
            self.phen_dict[pid]['transformation'] = 'sqrt(' + self.phen_dict[pid]['transformation'] + ')'
        self.phen_dict[pid]['values'] = vals.tolist()
        return True


    def ascombe_transform(self, pid, **kwargs):
        a = sp.array(self.phen_dict[pid]['values'])
        vals = 2.0 * sp.sqrt(a + 3.0 / 8.0)
        if not self.phen_dict[pid]['transformation']:
            self.phen_dict[pid]['raw_values'] = self.phen_dict[pid]['values']
            self.phen_dict[pid]['transformation'] = 'ascombe'
        else:
            self.phen_dict[pid]['transformation'] = 'ascombe(' + self.phen_dict[pid]['transformation'] + ')'
        self.phen_dict[pid]['values'] = vals.tolist()
        return True


    def sqr_transform(self, pid, method='standard'):
        a = sp.array(self.phen_dict[pid]['values'])
        if method == 'standard':
            vals = ((a - min(a)) + 0.1 * sp.std(a)) * ((a - min(a)) + 0.1 * sp.std(a))
        else:
            vals = a * a
        if not self.phen_dict[pid]['transformation']:
            self.phen_dict[pid]['raw_values'] = self.phen_dict[pid]['values']
            self.phen_dict[pid]['transformation'] = 'sqr'
        else:
            self.phen_dict[pid]['transformation'] = 'sqr(' + self.phen_dict[pid]['transformation'] + ')'
        self.phen_dict[pid]['values'] = vals.tolist()
        return True

    def exp_transform(self, pid, method='standard'):
        a = sp.array(self.phen_dict[pid]['values'])
        if method == 'standard':
            vals = sp.exp((a - min(a)) + 0.1 * sp.std(a))
        else:
            vals = sp.exp(a)
        if not self.phen_dict[pid]['transformation']:
            self.phen_dict[pid]['raw_values'] = self.phen_dict[pid]['values']
            self.phen_dict[pid]['transformation'] = 'exp'
        else:
            self.phen_dict[pid]['transformation'] = 'exp(' + self.phen_dict[pid]['transformation'] + ')'
        self.phen_dict[pid]['values'] = vals.tolist()
        return True

    def arcsin_sqrt_transform(self, pid, verbose=False):
        a = sp.array(self.phen_dict[pid]['values'])
        if min(a) < 0 or max(a) > 1:
            if verbose:
                print 'Some values are outside of range [0,1], hence skipping transformation!'
            return False
        else:
            vals = sp.arcsin(sp.sqrt(a))
        if not self.phen_dict[pid]['transformation']:
            self.phen_dict[pid]['raw_values'] = self.phen_dict[pid]['values']
            self.phen_dict[pid]['transformation'] = 'arcsin_sqrt'
        else:
            self.phen_dict[pid]['transformation'] = 'arcsin_sqrt(' + self.phen_dict[pid]['transformation'] + ')'
        self.phen_dict[pid]['values'] = vals.tolist()
        return True

    def box_cox_transform(self, pid, lambda_range=(-2.0, 2.0), lambda_increment=0.1, verbose=False, method='standard'):
        """
        Performs the Box-Cox transformation, over different ranges, picking the optimal one w. respect to normality.
        """
        from scipy import stats
        a = sp.array(self.phen_dict[pid]['values'])
        if method == 'standard':
            vals = (a - min(a)) + 0.1 * sp.std(a)
        else:
            vals = a
        sw_pvals = []
        lambdas = sp.arange(lambda_range[0], lambda_range[1] + lambda_increment, lambda_increment)
        for l in lambdas:
            if l == 0:
                vs = sp.log(vals)
            else:
                vs = ((vals ** l) - 1) / l
            r = stats.shapiro(vs)
            if sp.isfinite(r[0]):
                pval = r[1]
            else:
                pval = 0.0
            sw_pvals.append(pval)
        print sw_pvals
        i = sp.argmax(sw_pvals)
        l = lambdas[i]
        if l == 0:
            vs = sp.log(vals)
        else:
            vs = ((vals ** l) - 1) / l
        if not self.phen_dict[pid]['transformation']:
            self.phen_dict[pid]['raw_values'] = self.phen_dict[pid]['values']
            self.phen_dict[pid]['transformation'] = 'box-cox'
        else:
            self.phen_dict[pid]['transformation'] = 'box-cox(' + self.phen_dict[pid]['transformation'] + ')'
        self.phen_dict[pid]['values'] = vs.tolist()
        if verbose:
            print 'optimal lambda was %0.1f' % l
        return True





    def transform_pids(self, pids=None, trans_type='most_normal', method='standard'):
        if not pids:
            pids = self.get_pids()
        return [self.transform(pid, trans_type=trans_type, method=method) for pid in pids]


    def transform(self, pid, trans_type, method='standard', verbose=False):
        if verbose:
            print 'Transforming phenotypes: %s' % trans_type
        if trans_type == 'sqrt':
            self.sqrt_transform(pid, method=method)
        elif trans_type == 'ascombe':
            self.ascombe_transform(pid, method=method)
        elif trans_type == 'log':
            self.log_transform(pid, method=method)
        elif trans_type == 'sqr':
            self.sqr_transform(pid, method=method)
        elif trans_type == 'exp':
            self.exp_transform(pid, method=method)
        elif trans_type == 'arcsin_sqrt':
            self.arcsin_sqrt_transform(pid)
        elif trans_type == 'box_cox':
            self.box_cox_transform(pid, verbose=verbose)
        elif trans_type == 'most_normal':
            trans_type, shapiro_pval = self.most_normal_transformation(pid, verbose=verbose)
        elif trans_type == 'none':
            pass
        else:
            raise Exception('Transformation unknown')
        return trans_type


    def revert_to_raw_values(self, pid):
        if not self.phen_dict[pid]['transformation']:
            warnings.warn('Phenotype values are already raw..')
        else:
            self.phen_dict[pid]['transformation'] = None
            self.phen_dict[pid]['values'] = self.phen_dict[pid]['raw_values']


    def most_normal_transformation(self, pid, trans_types=['none', 'sqrt', 'log', 'sqr', 'exp', 'arcsin_sqrt'],
                perform_trans=True, verbose=False):
        """
        Performs the transformation which results in most normal looking data, according to Shapiro-Wilk's test
        """
        #raw_values = self.phen_dict[pid]['values']
        from scipy import stats
        shapiro_pvals = []
        for trans_type in trans_types:
            if trans_type != 'none':
                if not self.transform(pid, trans_type=trans_type):
                    continue
            phen_vals = self.get_values(pid)
            #print 'sp.inf in phen_vals:', sp.inf in phen_vals
            if sp.inf in phen_vals:
                pval = 0.0
            else:
                r = stats.shapiro(phen_vals)
                if sp.isfinite(r[0]):
                    pval = r[1]
                else:
                    pval = 0.0
            shapiro_pvals.append(pval)
            #self.phen_dict[pid]['values'] = raw_values
            if trans_type != 'none':
                self.revert_to_raw_values(pid)
        argmin_i = sp.argmax(shapiro_pvals)
        trans_type = trans_types[argmin_i]
        shapiro_pval = shapiro_pvals[argmin_i]
        if perform_trans:
            self.transform(pid, trans_type=trans_type)
        if verbose:
            print "The most normal-looking transformation was %s, with a Shapiro-Wilk's p-value of %0.6f" % \
                (trans_type, shapiro_pval)
        return trans_type, shapiro_pval


    def normalize_values(self, pids):
        for pid in pids:
            a = sp.array(self.phen_dict[pid]['values'])
            v = sp.var(self.get_avg_values(pid), ddof=1)
            vals = a / sp.sqrt(v)
            self.phen_dict[pid]['values'] = vals.tolist()


    def na_outliers(self, pids, iqr_threshold):
        raise NotImplementedError


    def filter_phenotypes(self, pids_to_keep):
        """
        Removes phenotypes.
        """
        self.phen_ids = pids_to_keep
        d = {}
        for pid in pids_to_keep:
            if pid in self.phen_dict:
                d[pid] = self.phen_dict[pid]
            else:
                print "skipping pid %d, since it's missing" % pid
        self.phen_dict = d


    def filter_phenotypes_w_few_values(self, min_num_vals=50):
        """
        Removes phenotypes.
        """
        d = {}
        for pid in self.phen_dict:
            if len(self.get_ecotypes(pid)) >= min_num_vals:
                d[pid] = self.phen_dict[pid]
        self.phen_dict = d
        self.phen_ids = sorted(d.keys())


    def filter_near_const_phens(self, min_num_diff=15):
        """
        """
        n1 = self.num_traits()
        pids_to_keep = []
        for pid in self.phen_ids:
            if not self.is_near_constant(pid, min_num_diff):
                pids_to_keep.append(pid)
        self.filter_phenotypes(pids_to_keep)
        n2 = self.num_traits()
        print '%d traits out of %d traits were filtered, leaving %d.' % (n1 - n2, n1, n2)



#    def filter_unique_ecotypes(self, ets, pids):
#        """
#        Removes ecotypes which are not in the given list of ecotypes.
#        """
#        if not pids:
#            pids = self.phen_ids
#        for pid in pids:
#            ecotypes = self.get_ecotypes(pid)
#            indices_to_keep = [i for i, et in enumerate(ecotypes) if et in ets ]
#            self.filter_ecotypes(indices_to_keep, pids=[pid])



    def filter_ecotypes(self, indices_to_keep, random_fraction=1, pids=None):
        """
        Removes the ecotypes from all data.
        """
        import random
        if not pids:
            pids = self.phen_ids
        for pid in pids:
            el = []
            vl = []
            d = self.phen_dict[pid]
            if d['transformation']:
                rvl = []
            if random_fraction < 1:
                indices = range(len(d['ecotypes']))
                indices_to_keep = sorted(random.sample(indices, int(len(d['ecotypes']) * random_fraction)))
            for i in indices_to_keep:
                el.append(d['ecotypes'][i])
                vl.append(d['values'][i])
                if d['transformation']:
                    rvl.append(d['raw_values'][i])
            self.phen_dict[pid]['ecotypes'] = el
            self.phen_dict[pid]['values'] = vl
            if d['transformation']:
                self.phen_dict[pid]['raw_values'] = rvl

    def filter_ecotypes_2(self, ecotypes_to_keep, pids=None):
        if not pids:
            pids = self.phen_ids
        unique_ets = set()
        for pid in pids:
            el = []
            vl = []
            d = self.phen_dict[pid]
            if d['transformation']:
                rvl = []
            for et in ecotypes_to_keep:
                if et in d['ecotypes']:
                    i = d['ecotypes'].index(et)
                    el.append(d['ecotypes'][i])
                    vl.append(d['values'][i])
                    if d['transformation']:
                        rvl.append(d['raw_values'][i])
                    unique_ets.add(et)
            self.phen_dict[pid]['ecotypes'] = el
            self.phen_dict[pid]['values'] = vl
            if d['transformation']:
                self.phen_dict[pid]['raw_values'] = rvl
        return list(unique_ets)


    def order_ecotypes(self, ets_map, pids=None):
        if not pids:
            pids = self.phen_ids
        for pid in pids:
            d = self.phen_dict[pid]
            ets = []
            vals = []
            if d['transformation']:
                rvl = []
            for i in ets_map:
                ets.append(d['ecotypes'][i])
                vals.append(d['values'][i])
                if d['transformation']:
                    rvl.append(d['raw_values'][i])
            self.phen_dict[pid]['ecotypes'] = ets
            self.phen_dict[pid]['values'] = vals
            if d['transformation']:
                self.phen_dict[pid]['raw_values'] = rvl



    def get_name(self, pid):
        return self.phen_dict[pid]['name']

    def get_names(self, pids=None):
        if not pids:
            pids = self.get_pids()
        return [self.phen_dict[pid]['name'] for pid in pids]

    def get_pids(self):
        pids = self.phen_dict.keys()
        pids.sort()
        return pids

    def get_values(self, pid):
        return self.phen_dict[pid]['values']

    def get_values_list(self, pids=None, clone=False):
        if not pids:
            pids = self.get_pids()
        if clone:
            return [self.get_values(pid)[:] for pid in pids]
        else:
            return [self.get_values(pid) for pid in pids]


    def get_avg_values(self, pid):
        d = self.get_avg_value_dict(pid)
        return d['values']

    def get_ecotypes(self, pid):
        return self.phen_dict[pid]['ecotypes']

    def get_incidence_matrix(self, pid):
        ets = sp.array(self.phen_dict[pid]['ecotypes'])
        unique_ets = []
        i = 0
        while i < len(ets):
            et = ets[i]
            unique_ets.append(et)
            while i < len(ets) and ets[i] == et: #The ecotypes are assumed to be sorted
                i += 1
#        unique_ets = sp.mat(sp.unique(ets))
        Z = sp.int8(sp.mat(ets).T == sp.mat(unique_ets))
        #print Z
        return Z


    def _get_ecotype_value_dict_(self, pid):
        ecotypes = self.get_ecotypes(pid)
        values = self.get_values(pid)
        d = {}
        for et in set(ecotypes):
            d[et] = {'values':[], 'rep_num':0}

        for et, val in it.izip(ecotypes, values):
            d[et]['values'].append(val)
            d[et]['rep_num'] += 1
        return d

    def get_avg_value_dict(self, pid):
        """
        Returns the average values, along with the ecotypes, and rep_number
        """
        d = self._get_ecotype_value_dict_(pid)
        ecotypes = d.keys()
        avg_vals = []
        rep_nums = []
        for et in d:
            avg_vals.append(sp.mean(d[et]['values']))
            rep_nums.append(d[et]['rep_num'])
        return {'name':self.get_name(pid) + '_avg', 'ecotypes':ecotypes, 'values':avg_vals, 'rep_nums': rep_nums}



    def convert_to_averages(self, pids=None):
        """
        Replaces phenotypes with replicates with their averages.
        """
        if not pids:
            pids = self.phen_dict.keys()
        for pid in pids:
            if len(set(self.phen_dict[pid]['ecotypes'])) == len(self.phen_dict[pid]['ecotypes']): continue
            phen_name = self.get_name(pid)
            trans = self.phen_dict[pid]['transformation']
            self.phen_dict[pid] = self.get_avg_value_dict(pid)
            self.phen_dict[pid]['name'] = phen_name
            self.phen_dict[pid]['transformation'] = trans



    def plot_histogram(self, pid, title=None , pdf_file=None, png_file=None, x_label=None, p_her=None,
            p_her_pval=None):

        if title:
            plt.figure(figsize=(6, 5.4))
            plt.axes([0.13, 0.11, 0.85, 0.82])
        else:
            plt.figure(figsize=(6, 4.8))
            plt.axes([0.13, 0.11, 0.85, 0.86])
        if x_label:
            plt.xlabel(x_label)
        phen_vals = self.get_values(pid)

        minVal = min(phen_vals)
        maxVal = max(phen_vals)
        x_range = maxVal - minVal
        histRes = plt.hist(phen_vals, bins=round(8 + 2 * sp.log(len(phen_vals))), alpha=0.7)
        y_max = max(histRes[0])
        plt.axis([minVal - 0.035 * x_range, maxVal + 0.035 * x_range, -0.035 * y_max, 1.19 * y_max])
        num_phen_vals = len(phen_vals)
        shapiro_pval = sp.stats.shapiro(phen_vals)[1]
        if p_her:
            if p_her_pval != None:
                st = "Num. of values: %d, herit.: %0.4f, herit. -log(p): %0.4f, transf.: %s" % \
                    (num_phen_vals, p_her, -sp.log10(p_her_pval), str(self.phen_dict[pid]['transformation']))
                plt.text(maxVal - 0.95 * x_range, 1.1 * y_max, st, size="xx-small")
            else:
                st = "Number of values: %d,  Pseudo-heritability: %0.4f,  Transformation: %s" % \
                    (num_phen_vals, p_her, str(self.phen_dict[pid]['transformation']))
                plt.text(maxVal - 0.95 * x_range, 1.1 * y_max, st, size="xx-small")
        else:
            st = "Number of values: %d, Transformation: %s" % (num_phen_vals, str(self.phen_dict[pid]['transformation']))
            plt.text(maxVal - 0.9 * x_range, 1.1 * y_max, st, size="x-small")
        plt.text(maxVal - 0.85 * x_range, 1.02 * y_max, "Shapiro-Wilk normality $p$-value: %0.6f" % shapiro_pval , size="x-small")
        #print max(histRes[0])
        plt.ylabel("Frequency")
        if title:
            plt.title(title)
        if pdf_file:
            plt.savefig(pdf_file, format="pdf")
        if png_file:
            plt.savefig(png_file, format="png", dpi=300)
        elif not pdf_file or png_file:
            plt.show()
        plt.clf()


    def write_to_file(self, file_name, delim=','):
        """
        Writes the object to a file.. (in the new format)
        """
        f = open(file_name, 'w')
        header = ['phenotype_id', 'phenotype_name', 'ecotype_id', 'value', 'replicate_id', ]
        f.write(delim.join(header) + '\n')
        for pid in self.phen_dict:
            d = self.phen_dict[pid]
            phen_name = d['name']
            ets_vals = zip(d['ecotypes'], d['values'])
            ets_vals.sort()
            last_et = -1
            for (et, val) in ets_vals:
                if et != last_et:
                    repl_id = 1
                else:
                    repl_id += 1
                l = map(str, [pid, phen_name, et, val, repl_id])
                f.write(delim.join(l) + '\n')
                last_et = et
        f.close()


    def get_correlations(self, pids=None):
        """
        Returns correlation matrix between traits
        
        All traits are used if pids is left empty.
        """
        import bisect
        if not pids:
            pids = sorted(self.phen_dict.keys())

        num_traits = len(pids)
        corr_mat = sp.ones((num_traits, num_traits))
        for i, pid1 in enumerate(pids):
            pd = self.get_avg_value_dict(pid1)
            ets1 = pd['ecotypes']
            pvs1 = pd['values']
            for j, pid2 in enumerate(pids[:i]):
                pd = self.get_avg_value_dict(pid2)
                ets2 = pd['ecotypes']
                pvs2 = pd['values']
                common_ets = set(ets1).intersection(set(ets2))
                ets_ix1 = map(ets1.index, common_ets)
                ets_ix2 = map(ets2.index, common_ets)
                vs1 = [pvs1[et_i] for et_i in ets_ix1]
                vs2 = [pvs2[et_i] for et_i in ets_ix2]
                corr_mat[i, j] = sp.corrcoef(vs1, vs2)[0, 1]
                corr_mat[j, i] = corr_mat[i, j]
        return corr_mat, pids

    def get_correlation(self, pid1, phed, pid2):
        """
        Returns the correlation with another phenotype_data object 
        """
        assert pid1 in self.phen_dict, 'phenotype ID %d missing in the self phed??' % pid1
        assert pid2 in phed.phen_dict, 'phenotype ID %d missing in the self phed??' % pid2
        pd = self.get_avg_value_dict(pid1)
        ets1 = pd['ecotypes']
        pvs1 = pd['values']
        pd = phed.get_avg_value_dict(pid2)
        ets2 = pd['ecotypes']
        pvs2 = pd['values']
        common_ets = set(ets1).intersection(set(ets2))
        ets_ix1 = map(ets1.index, common_ets)
        ets_ix2 = map(ets2.index, common_ets)
        vs1 = [pvs1[et_i] for et_i in ets_ix1]
        vs2 = [pvs2[et_i] for et_i in ets_ix2]
        return sp.corrcoef(vs1, vs2)[0, 1]




    def is_binary(self, pid):
        return len(sp.unique(self.phen_dict[pid]['values'])) == 2

    def is_constant(self, pid):
        return len(sp.unique(self.phen_dict[pid]['values'])) == 1

    def is_near_constant(self, pid, min_num_diff=10):
        vals = sp.array(self.phen_dict[pid]['values'])
        if sp.std(vals) > 0:
            vals = 50 * (vals - sp.mean(vals)) / sp.std(vals)
            vals = vals - vals.min() + 0.1
            b_counts = sp.bincount(sp.array(sp.around(vals), dtype='int'))
            b = b_counts.max() > len(vals) - min_num_diff
            return b
        else:
            return True




#    def plot_accession_map(self, pid, ecotypes=None, pdf_file=None, png_file=None, map_type='swedish',
#            color_by=None, cmap=None, title='',
#            with_colorbar=True,):
#        """
#        Plot accessions on a map.
#        
#        'color_by' is by default set to be the phenotype values.
#        """
#        import matplotlib
#        matplotlib.use("Agg")
#        import matplotlib.pyplot as plt
#        matplotlib.rcParams['backend'] = 'GTKAgg'
#        if not ecotypes:
#            ecotypes = self.phen_dict[pid]['ecotypes']
#        #eid = get_250K_accession_to_ecotype_dict(dict_key='ecotype_id')
#        eid = get_ecotype_id_info_dict()
#        lats = []
#        lons = []
#        acc_names = []
#        for e in ecotypes:
#            r = eid[int(e)]
#            acc_names.append(r[0])
#            try:
#                latitude = float(r[2])
#                longitude = float(r[3])
##                r = eid[str(e)]
##                latitude = float(r[5])
##                longitude = float(r[6])
#
#            except Exception, err_str:
#                print "Latitude and Longitude, not found?:", err_str
#                print 'Placing them in the Atlantic.'
#                latitude = 40
#                longitude = -20
#
#            lats.append(latitude)
#            lons.append(longitude)
#
#        from mpl_toolkits.basemap import Basemap
#        import numpy as np
#        from pylab import cm
#        if map_type == "global2":
#            plt.figure(figsize=(14, 12))
#            m = Basemap(width=21.e6, height=21.e6, projection='gnom', lat_0=76, lon_0=15)
#            m.drawparallels(np.arange(20, 90, 20))
#            m.drawmeridians(np.arange(-180, 180, 30))
#        elif map_type == 'global':
#
#            plt.figure(figsize=(16, 4))
#            plt.axes([0.02, 0.02, 0.96, 0.96])
#            m = Basemap(projection='cyl', llcrnrlat=10, urcrnrlat=80,
#                    llcrnrlon= -130, urcrnrlon=150, lat_ts=20, resolution='c')
#            m.drawparallels(np.arange(20, 90, 20))
#            m.drawmeridians(np.arange(-180, 180, 30))
#        elif map_type == 'europe':
#
#            plt.figure(figsize=(8, 6))
#            plt.axes([0.02, 0.02, 0.96, 0.96])
#            m = Basemap(projection='cyl', llcrnrlat=35, urcrnrlat=70,
#                    llcrnrlon= -15, urcrnrlon=40, lat_ts=20, resolution='h')
#            m.drawparallels(np.arange(30, 80, 10))
#            m.drawmeridians(np.arange(-20, 100, 10))
#            #m.bluemarble()
#        elif map_type == 'swedish':
#
#            plt.figure(figsize=(5, 6))
#            ax1 = plt.axes([0.05, 0.05, 0.75, 0.9])
#            m = Basemap(width=2800000, height=4000000, projection='cass', llcrnrlat=54, urcrnrlat=65,
#                    llcrnrlon=10, urcrnrlon=25, lat_ts=20, resolution='h', lon_0=17.5, lat_0=59.5, ax=ax1)
#            m.drawparallels(np.arange(30, 80, 10))
#            m.drawmeridians(np.arange(-20, 100, 10))
#            #m.bluemarble()
#        else:
#            raise Exception("map_type is invalid")
#
#        #m.drawmapboundary(fill_color='aqua')
#        m.drawcoastlines(zorder=0)
#        m.fillcontinents(zorder=1)
#        #m.fillcontinents(color='coral',lake_color='aqua')
#
#        xs = []
#        ys = []
#        for lon, lat in zip(lons, lats):
#            x, y = m(*np.meshgrid([lon], [lat]))
#            xs.append(float(x))
#            ys.append(float(y))
#
#        if not color_by:
#            color_vals = self.get_values(pid)
#        else:
#            color_vals = color_by
#        if len(color_vals) != len(self.get_ecotypes(pid)):
#            raise Exception("accessions and color_by_vals values don't match ! ")
#        if not cmap:
#            num_colors = len(set(color_vals))
#            if num_colors <= 10:
#                cmap = cm.get_cmap('jet', num_colors)
#            else:
#                cmap = cm.get_cmap('jet')
#        lws = [0] * len(xs)
#        plt.scatter(xs, ys, s=10, linewidths=lws, c=color_vals, cmap=cmap, alpha=0.7, zorder=2)
#        #plt.plot(xs, ys, 'o', color='r', alpha=0.5, zorder=2,)
#        if with_colorbar:
#            ax2 = plt.axes([0.84, 0.3, 0.05, 0.4])
#            plt.colorbar(ax=ax1, cax=ax2)
#        if title:
#            plt.title(title)
#        if pdf_file:
#            plt.savefig(pdf_file, format="pdf")
#        if png_file:
#            plt.savefig(png_file, format="png")
#        if not pdf_file and not png_file:
#            plt.show()
#
#        return ecotypes, acc_names, lats, lons


    def plot_marker_box_plot(self, pid, marker, m_accessions, m_position=None, m_chromosome=None, plot_file=None,
                plot_format='png', title=None, m_score=None):
        """
        Plots a box plot for the given binary marker and phenotype. 
        
        Assumes the marker is integer based.        
        Assumes the marker and the phenotype accessions are aligned.
        """
        phen_vals = self.get_values(pid)
        if len(m_accessions) != len(phen_vals):
            raise Exception

        nt_counts = sp.bincount(marker)
        if len(nt_counts) > 2:
            import warnings
            warnings.warn("More than 2 alleles, box-plot might be wrong?")

        allele_phen_val_dict = {}
        for nt in set(marker):
            allele_phen_val_dict[nt] = {'values':[], 'ecotypes':[]}

        for i, nt in enumerate(marker):
            allele_phen_val_dict[nt]['values'].append(phen_vals[i])
            if m_accessions:
                allele_phen_val_dict[nt]['ecotypes'].append(m_accessions[i])

        xs = []
        positions = []
        for nt in allele_phen_val_dict:
            positions.append(nt)
            xs.append(allele_phen_val_dict[nt]['values'])
        plt.figure()
        plt.boxplot(xs, positions=positions)
        min_val = min(phen_vals)
        max_val = max(phen_vals)
        val_range = max_val - min_val
        max_pos = max(positions)
        min_pos = min(positions)
        x_range = max_pos - min_pos
        plt.axis([min_pos - 0.5 * x_range, max_pos + 0.5 * x_range, min_val - val_range * 0.3, max_val + val_range * 0.3])
        plt.text(min_pos - 0.45 * x_range, min_val - 0.15 * val_range, "# of obs.: ", color='k')
        for i, (x, pos) in enumerate(it.izip(xs, positions)):
            plt.text(pos - 0.05, min_val - 0.15 * val_range, str(len(xs[i])), color='k')
        if m_score:
            plt.text(min_pos + 0.13 * x_range, max_val + 0.15 * val_range,
                '$-log_{10}$(p-value)/score: %0.2f' % m_score, color='k')
        if title:
            plt.title(title)
        elif m_chromosome and m_position:
            plt.title('%s : chromosome=%d, position=%d' % (self.get_name(pid), m_chromosome, m_position))
        if plot_file:
            plt.savefig(plot_file, format=plot_format)
        else:
            plt.show()
        plt.clf()


    def plot_marker_accessions_hist(self, pid, marker, m_accessions, plot_file=None, plot_format='png',
                m_position=None, m_chromosome=None, title=None, m_score=None):
        """
        A histogram displaying the phenotype values (ordered) on the y-axis, and the accession on the x-axis.
        """
        import matplotlib.cm as cm
        import matplotlib.colors as colors

        color_map = {}
        colors = ['r', 'm', 'b', 'g']
        proxy_rects = []
        labels = []
        for nt in set(marker):
            c = colors.pop()
            color_map[nt] = c
            proxy_rects.append(plt.Rectangle((0, 0), 1, 1, fc=c, alpha=0.6))
            labels.append("'%s' allele" % str(nt))


        phen_values = self.get_values(pid)
        ets = self.get_ecotypes(pid)
        l = zip(phen_values, ets, marker)
        l.sort(reverse=True)

        fig = plt.figure(figsize=(20, 10))
        ax = fig.add_axes([0.07, 0.15, 0.91, 0.82])
        x_range = len(l) - 0.2
        min_y = min(phen_values)
        max_y = max(phen_values)
        y_range = max_y - min_y

        for i, (phen_val, accession, nt) in enumerate(l):
            color = color_map[nt]
            rect = ax.bar(i, phen_val, 0.8, color=color, alpha=0.6)
        ax.axis([-x_range * 0.02, x_range * 1.02, min_y - 0.05 * y_range, max_y + 0.05 * y_range])

        ax.legend(proxy_rects, labels)
        ax.set_ylabel('Phenotype value')
        ax.set_xticks((sp.arange(len(ets)) + 0.4).tolist())
        ax.set_xticklabels(ets, rotation='vertical', fontsize='xx-small')
        ax.set_xlabel('Ecotype IDs')

        fig.savefig(plot_file, format=plot_format, dpi=300)


    def plot_phen_relatedness(self, k, k_accessions, plot_file_prefix, pids=None):
        import kinship
        import pylab
        import scipy as sp
        from scipy import linalg
        if not pids:
            pids = self.get_pids()
        self.convert_to_averages(pids)
        self.filter_ecotypes_2(k_accessions, pids)
        for pid in pids:
            ets = self.get_ecotypes(pid)
            vals = self.get_values(pid)
            k_m = kinship.prepare_k(k, k_accessions, ets)
            c = sp.sum((sp.eye(len(k_m)) - (1.0 / len(k_m)) * sp.ones(k_m.shape)) * sp.array(k_m))
            k_scaled = (len(k) - 1) * k / c
            p_her = self.get_pseudo_heritability(pid, k_m)
            x_list = []
            y_list = []
            for i in range(len(ets)):
                for j in range(i):
                    x_list.append(k_m[i, j])
                    y_list.append(vals[i] - vals[j])
            ys = sp.array(y_list)
            ys = ys * ys
            xs = sp.array(x_list)
            phen_name = self.get_name(pid)
            phen_name = phen_name.replace('<i>', '')
            phen_name = phen_name.replace('</i>', '')
            phen_name = phen_name.replace('+', '_plus_')
            phen_name = phen_name.replace('/', '_div_')
            file_name = plot_file_prefix + '_%d_%s.png' % (pid, phen_name)
            pylab.figure()
            pylab.plot(xs, ys, 'k.', alpha=0.2)
            pylab.xlabel('Relatedness')
            pylab.ylabel('Squared phenotypic difference')
            #Plot regression line
            Y_mat = sp.mat(ys).T
            X_mat = sp.hstack((sp.mat(sp.ones(len(xs))).T, sp.mat(xs).T))
            (betas, residues, rank, s) = linalg.lstsq(X_mat, Y_mat)
            x_min, x_max = pylab.xlim()
            pylab.plot([x_min, x_max], [betas[0] + x_min * betas[1], betas[0] + x_max * betas[1]])
            corr = sp.corrcoef(xs, ys)[0, 1]
            y_min, y_max = pylab.ylim()
            x_range = x_max - x_min
            y_range = y_max - y_min
            pylab.axis([x_min - 0.025 * x_range, x_max + 0.025 * x_range,
                    y_min - 0.025 * y_range, y_max + 0.15 * y_range])
            pylab.text(x_min + 0.1 * x_range, y_max + 0.03 * y_range, 'Correlation: %0.4f' % (corr))
            pylab.text(x_min + 0.5 * x_range, y_max + 0.03 * y_range, 'Pseudo-heritability: %0.4f' % (p_her))
            pylab.savefig(file_name)
            del k_m
            del k_scaled






def parse_phenotype_file(file_name=None, file_object=None, delim=',', file_format='guess', with_db_ids=True):
    """
    Parses a phenotype file, and returns a new phenotype_data object.
    
    File format types:
    
        new - Everything in a long sequence.
        old - Well, old phenotype file_format.
        guess - Guesses the file_format
    """
    phen_dict = {}
    if file_object:
        f = file_object
    else:
        f = open(file_name, 'rU')
    header = f.next()
    if len(header.split(delim)) < 2:
        test_delims = [',', '\t']
        for n_delim in test_delims:
            if len(header.split(n_delim)) > 2:
                delim = n_delim
                break
        else:
            raise Exception('Problems with delimiters', delim, test_delims)
    print 'Using "%s" as delimiter' % delim

    if file_format == 'guess':
        header = map(str.strip, header.split(delim))
        print header
#        if len(header) != 5:
#            print 'Guessing old format.'
#            #print len(header), header
#            file_format = 'old'
        if 'phenotype_id' in header or 'replicate_id' in header:
            print 'Guessing new format.'
            file_format = 'new'
        else:
            print 'Guessing old format.'
            file_format = 'old'
#            v = int(raw_input('File format is ambiguous:\n (1) - old format\n (2) - new format\n (3) - exit\n:'))
#            if v == 1:
#                file_format = 'old'
#            elif v == 2:
#                file_format = 'new'
#            else:
#                sys.exit()

    if file_format == 'old':
        if with_db_ids:
            pids = [int(l.split('_')[0]) for l in header[1:]]
        else:
            pids = range(1, len(header))
        for i, pid in enumerate(pids):
            phen_dict[pid] = {'ecotypes':[], 'values':[], 'name':header[i + 1]}
        for line in f:
            l = map(str.strip, line.split(delim))
            ecotype = l[0]
            del l[0]
            for i, v in enumerate(l):
                pid = pids[i]
                if v != 'NA':
                    phen_dict[pid]['ecotypes'].append(ecotype)
                    phen_dict[pid]['values'].append(float(v))

    elif file_format == 'new':
        last_pid = -1
        for line in f:
            l = line.split(delim)
            pid = int(l[0])
            if pid != last_pid:
                if last_pid != -1:
                    phen_dict[last_pid] = d
                phen_name = l[1]
                d = {'name':phen_name, 'ecotypes':[], 'values':[]}
            d['ecotypes'].append(l[2])
            d['values'].append(float(l[3]))
            last_pid = pid
        phen_dict[last_pid] = d

    #print phen_dict
    if not file_object:
        f.close()

    return phenotype_data(phen_dict=phen_dict, phen_ids=phen_dict.keys())



def get_250K_accession_to_ecotype_dict(fn='at_data/call_method_75_info.tsv', dict_key='nativename'):
    """
    """
    delim = '\t'
    f = open(fn)
    header = f.next().split(delim)
    #print header
    i = header.index(dict_key)
    info_table = []
    ret_dict = {}
    for line in f:
        l = line.split(delim)
        info_table.append(l)
        ret_dict[l[i].strip().lower()] = l
    return ret_dict




if __name__ == '__main__':
   print "Done!"
