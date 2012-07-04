from env import *
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






    def get_db_pid(self, pid, conn=None):
        """
        Retrieves the DB pid, using the phenotype name.
        """
        phen_name = self.get_name(pid)
        import dbutils
        if not conn:
            conn = dbutils.connect_to_default_lookup()
        cursor = conn.cursor()
        sql_statement = "SELECT id FROM stock_250k.phenotype_method WHERE short_name='%s'" % (phen_name)
        print sql_statement
        numRows = int(cursor.execute(sql_statement))
        row = cursor.fetchone()
        if row:
            db_pid = int(row[0])
        else:
            print "No id found in DB for phenotype:%s" % (phen_name)
            db_pid = None
        cursor.close()
        if not conn:
            conn.close()
        return db_pid


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


    def convert_to_tg_ets(self, pids=None):
        """
        Convert the ecotypes to target-ecotypes.
        """
        tg_ets_map = get_tg_ecotype_map()
        if not pids:
            pids = self.phen_dict.keys()
        f = lambda et: str(int(tg_ets_map[int(et)][0])) if int(et) in tg_ets_map else 'NA'
        for pid in pids:
            self.phen_dict[pid]['ecotypes'] = map(f, self.phen_dict[pid]['ecotypes'])


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

    def insert_into_db(self, pids=None, phenotype_scoring='',
               method_description='', growth_condition='', biology_category_id='',
               citations='', data_description='', transformation_description=None,
               method_id=None, data_type=None, comment='', dry_run=False, update_only=False):
        """
        Inserts phenotypes into avg DB, but assumes that values are averages.
        (Assumes no transformation)
        """

        if not pids:
            pids = self.phen_dict.keys()
        import dbutils
        conn = dbutils.connect_to_default_insert("stock_250k")
        cursor = conn.cursor()
        mids = []
        for pid in pids:
            try:
                phen_name = self.get_name(pid)
                cur_method_id = self.get_db_pid(pid, conn=conn)
                phen_vals = self.get_values(pid)
                ecotypes = self.get_ecotypes(pid)
                no_of_accessions = len(ecotypes)
                print cur_method_id
                if cur_method_id:
                    print phen_name, 'is already in DB.'
                    sql_statement = "UPDATE stock_250k.phenotype_method SET no_of_accessions=%d, biology_category_id=%d WHERE short_name='%s'" % (no_of_accessions, biology_category_id, phen_name)
                    print sql_statement
                    cursor.execute(sql_statement)
                    if update_only:
                            print 'Updating values.'
                    else:
                            print 'Deleting existing phenotype values'
                            sql_statement = "DELETE FROM stock_250k.phenotype_avg WHERE method_id=%d" % (cur_method_id)
                            print sql_statement
                            cursor.execute(sql_statement)

                else:
                    if not data_type:
                        if self.is_binary(pid):
                            data_type = 'binary'
                        else:
                            data_type = 'quantitative'
                    print "Inserting phenotype %s into DB." % phen_name
                    if method_id:
                        sql_statement = "INSERT INTO stock_250k.phenotype_method  (id, short_name, only_first_96, no_of_accessions," + \
                             " biology_category_id,phenotype_scoring, method_description, growth_condition, data_description, comment," + \
                             " data_type, transformation_description) VALUES (" + str(method_id) + ", '" + phen_name + "', false, " + \
                             str(no_of_accessions) + ", " + str(biology_category_id) + ", '" + phenotype_scoring + "', '" + method_description + \
                             "', '" + growth_condition + "', '" + data_description + "', '" + comment + "', '" + data_type + "', '" + \
                             str(transformation_description) + "')"
                    else:
                        sql_statement = "INSERT INTO stock_250k.phenotype_method  (short_name, only_first_96, no_of_accessions," + \
                             " biology_category_id,phenotype_scoring, method_description, growth_condition, data_description," + \
                             " comment, data_type, transformation_description) VALUES ('" + phen_name + "', false, " + \
                             str(no_of_accessions) + ", " + str(biology_category_id) + ", '" + phenotype_scoring + "', '" + method_description + \
                             "', '" + growth_condition + "', '" + data_description + "', '" + comment + "', '" + data_type + "', '" + \
                             str(transformation_description) + "')"

                    print sql_statement
                    numRows = int(cursor.execute(sql_statement))
                    row = cursor.fetchone()
                    if row:
                        print row
                    cur_method_id = self.get_db_pid(pid, conn=conn)

                print "Inserting values for phenotype method_id = %d" % cur_method_id
                for e_i, val in zip(ecotypes, phen_vals):
                    sql_statement = "INSERT INTO stock_250k.phenotype_avg  (ecotype_id, value, ready_for_publication, " + \
                            "method_id, transformed_value) VALUES ( " + str(e_i) + ", " + str(val) + ", 0, " + \
                            str(cur_method_id) + ", " + str(val) + " )"
                    try:
                        print sql_statement
                        numRows = int(cursor.execute(sql_statement))
                        row = cursor.fetchone()
                        if row:
                            print row
                    except Exception, err_str:
                        print 'Failed at inserting data:', err_str
                        print 'Trying to update.'
                        sql_statement = "UPDATE stock_250k.phenotype_avg SET value=%f, transformed_value=%f \
                                 WHERE ecotype_id=%d and method_id=%d" % (val, val, int(e_i), cur_method_id)
                        print sql_statement
                        numRows = int(cursor.execute(sql_statement))
                        row = cursor.fetchone()
                        if row:
                            print row


                if not dry_run:
                    print "Committing"
                    conn.commit()
                    mids.append(cur_method_id)
            except Exception, err_str:
                print 'Failed at inserting phentoype ID', pid, ':', err_str

        cursor.close ()
        conn.close ()
        return mids




    def plot_accession_map(self, pid, ecotypes=None, pdf_file=None, png_file=None, map_type='swedish',
            color_by=None, cmap=None, title='',
            with_colorbar=True,):
        """
        Plot accessions on a map.
        
        'color_by' is by default set to be the phenotype values.
        """
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        matplotlib.rcParams['backend'] = 'GTKAgg'
        if not ecotypes:
            ecotypes = self.phen_dict[pid]['ecotypes']
        #eid = get_250K_accession_to_ecotype_dict(dict_key='ecotype_id')
        eid = get_ecotype_id_info_dict()
        lats = []
        lons = []
        acc_names = []
        for e in ecotypes:
            r = eid[int(e)]
            acc_names.append(r[0])
            try:
                latitude = float(r[2])
                longitude = float(r[3])
#                r = eid[str(e)]
#                latitude = float(r[5])
#                longitude = float(r[6])

            except Exception, err_str:
                print "Latitude and Longitude, not found?:", err_str
                print 'Placing them in the Atlantic.'
                latitude = 40
                longitude = -20

            lats.append(latitude)
            lons.append(longitude)

        from mpl_toolkits.basemap import Basemap
        import numpy as np
        from pylab import cm
        if map_type == "global2":
            plt.figure(figsize=(14, 12))
            m = Basemap(width=21.e6, height=21.e6, projection='gnom', lat_0=76, lon_0=15)
            m.drawparallels(np.arange(20, 90, 20))
            m.drawmeridians(np.arange(-180, 180, 30))
        elif map_type == 'global':

            plt.figure(figsize=(16, 4))
            plt.axes([0.02, 0.02, 0.96, 0.96])
            m = Basemap(projection='cyl', llcrnrlat=10, urcrnrlat=80,
                    llcrnrlon= -130, urcrnrlon=150, lat_ts=20, resolution='c')
            m.drawparallels(np.arange(20, 90, 20))
            m.drawmeridians(np.arange(-180, 180, 30))
        elif map_type == 'europe':

            plt.figure(figsize=(8, 6))
            plt.axes([0.02, 0.02, 0.96, 0.96])
            m = Basemap(projection='cyl', llcrnrlat=35, urcrnrlat=70,
                    llcrnrlon= -15, urcrnrlon=40, lat_ts=20, resolution='h')
            m.drawparallels(np.arange(30, 80, 10))
            m.drawmeridians(np.arange(-20, 100, 10))
            #m.bluemarble()
        elif map_type == 'swedish':

            plt.figure(figsize=(5, 6))
            ax1 = plt.axes([0.05, 0.05, 0.75, 0.9])
            m = Basemap(width=2800000, height=4000000, projection='cass', llcrnrlat=54, urcrnrlat=65,
                    llcrnrlon=10, urcrnrlon=25, lat_ts=20, resolution='h', lon_0=17.5, lat_0=59.5, ax=ax1)
            m.drawparallels(np.arange(30, 80, 10))
            m.drawmeridians(np.arange(-20, 100, 10))
            #m.bluemarble()
        else:
            raise Exception("map_type is invalid")

        #m.drawmapboundary(fill_color='aqua')
        m.drawcoastlines(zorder=0)
        m.fillcontinents(zorder=1)
        #m.fillcontinents(color='coral',lake_color='aqua')

        xs = []
        ys = []
        for lon, lat in zip(lons, lats):
            x, y = m(*np.meshgrid([lon], [lat]))
            xs.append(float(x))
            ys.append(float(y))

        if not color_by:
            color_vals = self.get_values(pid)
        else:
            color_vals = color_by
        if len(color_vals) != len(self.get_ecotypes(pid)):
            raise Exception("accessions and color_by_vals values don't match ! ")
        if not cmap:
            num_colors = len(set(color_vals))
            if num_colors <= 10:
                cmap = cm.get_cmap('jet', num_colors)
            else:
                cmap = cm.get_cmap('jet')
        lws = [0] * len(xs)
        plt.scatter(xs, ys, s=10, linewidths=lws, c=color_vals, cmap=cmap, alpha=0.7, zorder=2)
        #plt.plot(xs, ys, 'o', color='r', alpha=0.5, zorder=2,)
        if with_colorbar:
            ax2 = plt.axes([0.84, 0.3, 0.05, 0.4])
            plt.colorbar(ax=ax1, cax=ax2)
        if title:
            plt.title(title)
        if pdf_file:
            plt.savefig(pdf_file, format="pdf")
        if png_file:
            plt.savefig(png_file, format="png")
        if not pdf_file and not png_file:
            plt.show()

        return ecotypes, acc_names, lats, lons


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
        import linear_models as lm
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
            k_m = lm.prepare_k(k, k_accessions, ets)
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



def get_phenotypes_from_db(pids):
    """
    Returns the phenotype using the average table in the DB.
    """
    import dbutils
    conn = dbutils.connect_to_default_lookup()
    cursor = conn.cursor ()
    phen_dict = {}
    for pid in pids:
        #Get the phenotypic values
        sql_statement = "\
        SELECT DISTINCT pa.ecotype_id, pa.value, pm.short_name\
        FROM stock_250k.phenotype_avg pa, stock_250k.phenotype_method pm \
        WHERE pa.method_id=pm.id AND pm.id=%d \
        ORDER BY pa.method_id, pa.ecotype_id" % (pid)
        #print sql_statement
        num_rows = int(cursor.execute(sql_statement))
        num_values = 0
        while(1):
            row = cursor.fetchone()
            if not row:
                break;
            phen_name = row[2]
            if pid not in phen_dict:
                phen_dict[pid] = {'name':phen_name, 'ecotypes':[], 'values':[]}

            if row[1] != None:
                phen_dict[pid]['values'].append(float(row[1]))
                phen_dict[pid]['ecotypes'].append(str(row[0]))
                num_values += 1
        if pid in phen_dict:
            print 'Loaded %d values for phenotype: %s' % (num_values, phen_name)

    conn.close()
    return phenotype_data(phen_dict)




def get_all_phenotypes_from_db(file_name=None):
    """
    Fetches all raw phenotype data from the DB, including replicates.
    
    Checks whether averages are consistent with replicates.
    """
    import dbutils
    conn = dbutils.connect_to_default_lookup()
    cursor = conn.cursor ()

    #Get the phenotypic values
    sql_statement = "\
    SELECT DISTINCT pa.method_id, pa.ecotype_id, pa.value, pm.short_name\
    FROM stock_250k.phenotype_avg pa, stock_250k.phenotype_method pm \
    WHERE pa.method_id=pm.id ORDER BY pa.method_id, pa.ecotype_id"
    #print sql_statement
    numRows = int(cursor.execute(sql_statement))

    phen_dict = {}
    num_values = 0
    while(1):
        row = cursor.fetchone()
        if not row:
            break;
        pid = int(row[0])
        phen_name = row[3]
        if pid not in phen_dict:
            phen_dict[pid] = {'name':phen_name, 'ecotypes':[], 'values':[]}

        if row[2] != None:
            phen_dict[pid]['values'].append(float(row[2]))
            phen_dict[pid]['ecotypes'].append(str(row[1]))
            num_values += 1

    #print 'Values for %d phenotypes were fetched from DB.' % len(phen_dict)
    #print "%d values in total." % num_values

    avg_phend = phenotype_data(phen_dict)

    #Now fetch the replicates
#    sql_statement = "\
#    SELECT DISTINCT p.method_id, p.ecotype_id, p.value, p.replicate, pm.short_name\
#    FROM stock_250k.phenotype p, stock_250k.phenotype_method pm \
#    WHERE p.method_id=pm.id ORDER BY p.method_id, p.ecotype_id"
    sql_statement = "\
    SELECT DISTINCT p.method_id, ei.tg_ecotypeid, p.value, p.replicate, pm.short_name\
    FROM stock_250k.phenotype p, stock_250k.phenotype_method pm, stock.ecotypeid2tg_ecotypeid ei \
    WHERE ei.ecotypeid=p.ecotype_id and p.method_id=pm.id ORDER BY p.method_id, ei.tg_ecotypeid"
    #print sql_statement
    numRows = int(cursor.execute(sql_statement))
    repl_phen_dict = {}
    num_values = 0
    while(1):
        row = cursor.fetchone()
        if not row:
            break;
        pid = int(row[0])
        phen_name = row[4]
        if pid not in repl_phen_dict:
            repl_phen_dict[pid] = {'name':phen_name, 'ecotypes':[], 'values':[], 'replicate':[]}

        if row[2] != None:
            repl_phen_dict[pid]['values'].append(float(row[2]))
            repl_phen_dict[pid]['ecotypes'].append(str(row[1]))
            repl_phen_dict[pid]['replicate'].append(int(row[3]))
            num_values += 1

    #print 'Values for %d phenotypes were fetched from DB.' % len(repl_phen_dict)
    #print "%d values in total." % num_values

    cursor.close ()
    conn.close ()
    #cm72_dict = get_250K_accession_to_ecotype_dict(dict_key='ecotype_id')
    #print cm72_dict
    repl_phend = phenotype_data(repl_phen_dict)
    no_problem_pids = []
    for pid in phen_dict:
        if pid in repl_phen_dict:
            d1 = avg_phend.get_avg_value_dict(pid)
            d2 = repl_phend.get_avg_value_dict(pid)
            ets1 = d1['ecotypes']
            ets2 = d2['ecotypes']
            vls1 = d1['values']
            vls2 = d2['values']
            l1 = zip(ets1, vls1)
            l2 = zip(ets2, vls2)
            l1.sort()
            l2.sort()
            et_set1 = set(ets1)
            et_set2 = set(ets2)
            phen_name = phen_dict[pid]['name']
            #print phen_name, ', pid:', pid, ', et_set1==et_set2 :', et_set1 == et_set2
            common_ets = et_set1.intersection(et_set2)
            num_diffs = 0
            for et in common_ets:
                t = (sp.var(d1['values']) + sp.var(d2['values'])) / 100.0
                i1 = d1['ecotypes'].index(et)
                i2 = d2['ecotypes'].index(et)
                if (d1['values'][i1] - d2['values'][i2]) > t:
                    #print d1['values'][i1], d2['values'][i2]
                    num_diffs += 1
            if not et_set1 == et_set2:
                #print len(et_set1), len(et_set2)
                s = et_set1.symmetric_difference(et_set2)
                #print 'et_set1.difference(et_set2):', et_set1.difference(et_set2)
                #print 'et_set2.difference(et_set1):', et_set2.difference(et_set1)
                if len(s) < 8 and num_diffs == 0:
                    #print 'No major problems'
                    no_problem_pids.append(pid)
                    continue
#                for ei in s:
#                    if str(ei) in cm72_dict:
#                        print cm72_dict[str(ei)]
#                    else:
#                        print ei, 'is not in 250K data!'
            elif num_diffs == 0:
                #print 'No problems?'
                no_problem_pids.append(pid)
    #print no_problem_pids

    merged_phen_dict = {}
    for pid in phen_dict:
        if pid in no_problem_pids:
            d = {'name':phen_dict[pid]['name'], 'values': repl_phen_dict[pid]['values'],
                'ecotypes': repl_phen_dict[pid]['ecotypes']}

            et_set1 = set(phen_dict[pid]['ecotypes'])
            et_set2 = set(repl_phen_dict[pid]['ecotypes'])
            for et in et_set1 - et_set2:
                vi = phen_dict[pid]['ecotypes'].index(et)
                val = phen_dict[pid]['values'][vi]
                d['ecotypes'].append(et)
                d['values'].append(val)
            merged_phen_dict[pid] = d
        else:
            merged_phen_dict[pid] = phen_dict[pid]


    phend = phenotype_data(phen_dict=merged_phen_dict)
    if file_name:
        phend.write_to_file(file_name)
    return phend



def _getFirst96Ecotypes_(host="papaya.usc.edu", user="bvilhjal", passwd="*rri_bjarni@usc"):
    import MySQLdb
    print "Connecting to db, host=" + host
    if not user:
        sys.stdout.write("Username: ")
        user = sys.stdin.readline().rstrip()
    if not passwd:
        import getpass
        passwd = getpass.getpass()
    try:
        conn = MySQLdb.connect (host=host, user=user, passwd=passwd, db="at")
    except MySQLdb.Error, e:
        print "Error %d: %s" % (e.args[0], e.args[1])
        sys.exit (1)
    cursor = conn.cursor ()
    print "Fetching data"
    ecotypes = []
    numRows = int(cursor.execute("select distinct e2te.tg_ecotypeid, c2010.accession_id, e2te.nativename, e2te.stockparent from at.complete_2010_strains_in_stock c2010, stock.ecotypeid2tg_ecotypeid e2te where e2te.ecotypeid=c2010.ecotypeid and e2te.stockparent = c2010.stockparent and c2010.accession_id<97 order by c2010.accession_id"))
    i = 0
    while(1):
        row = cursor.fetchone()
        if not row:
            break;
        ecotypes.append(int(row[0]))
    cursor.close ()
    conn.close ()
    return ecotypes

def get_199_ecotypes():
    """
    Result of a union of all the "192" accessions for the phenotypes to be used for the GWA 2009 paper.
    """
    return ['100000', '5837', '6008', '6009', '6016', '6040', '6042', '6043', '6046', '6064', '6074', '6243', '6709', '6897', '6898', '6899', '6900', '6901', '6903', '6904', '6905', '6906', '6907', '6908', '6909', '6910', '6911', '6913', '6914', '6915', '6916', '6917', '6918', '6919', '6920', '6921', '6922', '6923', '6924', '6926', '6927', '6928', '6929', '6930', '6931', '6932', '6933', '6936', '6937', '6939', '6940', '6942', '6943', '6944', '6945', '6946', '6951', '6956', '6957', '6958', '6959', '6960', '6961', '6962', '6963', '6964', '6965', '6966', '6967', '6968', '6969', '6970', '6971', '6972', '6973', '6974', '6975', '6976', '6977', '6978', '6979', '6980', '6981', '6982', '6983', '6984', '6985', '6988', '7000', '7014', '7033', '7062', '7064', '7081', '7094', '7123', '7147', '7163', '7231', '7255', '7275', '7282', '7296', '7306', '7323', '7346', '7418', '7424', '7438', '7460', '7461', '7477', '7514', '7515', '7516', '7517', '7518', '7519', '7520', '7521', '7522', '7523', '7524', '7525', '7526', '8213', '8214', '8215', '8222', '8230', '8231', '8233', '8235', '8236', '8237', '8239', '8240', '8241', '8242', '8243', '8245', '8247', '8249', '8254', '8256', '8258', '8259', '8264', '8265', '8266', '8270', '8271', '8274', '8275', '8283', '8284', '8285', '8290', '8296', '8297', '8300', '8306', '8310', '8311', '8312', '8313', '8314', '8323', '8325', '8326', '8329', '8334', '8335', '8337', '8343', '8351', '8353', '8354', '8357', '8365', '8366', '8369', '8374', '8376', '8378', '8387', '8388', '8389', '8395', '8411', '8412', '8420', '8422', '8423', '8424', '8426', '8430', '9057', '9058']
    #['100000', '5837', '6008', '6009', '6016', '6024', '6039', '6040', '6042', '6043', '6046', '6064', '6074', '6088', '6243', '6709', '6830', '6897', '6898', '6899', '6900', '6901', '6903', '6904', '6905', '6906', '6907', '6908', '6909', '6910', '6911', '6913', '6914', '6915', '6916', '6917', '6918', '6919', '6920', '6921', '6922', '6923', '6924', '6926', '6927', '6928', '6929', '6930', '6931', '6932', '6933', '6936', '6937', '6938', '6939', '6940', '6942', '6943', '6944', '6945', '6946', '6951', '6956', '6957', '6958', '6959', '6960', '6961', '6962', '6963', '6964', '6965', '6966', '6967', '6968', '6969', '6970', '6971', '6972', '6973', '6974', '6975', '6976', '6977', '6978', '6979', '6980', '6981', '6982', '6983', '6984', '6985', '6988', '7000', '7014', '7033', '7062', '7064', '7081', '7094', '7123', '7147', '7163', '7231', '7255', '7275', '7282', '7296', '7306', '7323', '7327', '7329', '7346', '7418', '7424', '7438', '7460', '7461', '7477', '7514', '7515', '7516', '7517', '7518', '7519', '7520', '7521', '7522', '7523', '7524', '7525', '7526', '8213', '8214', '8215', '8222', '8230', '8231', '8233', '8234', '8235', '8236', '8237', '8239', '8240', '8241', '8242', '8243', '8244', '8245', '8246', '8247', '8249', '8254', '8256', '8257', '8258', '8259', '8264', '8265', '8266', '8270', '8271', '8274', '8275', '8283', '8284', '8285', '8290', '8296', '8297', '8300', '8304', '8306', '8307', '8310', '8311', '8312', '8313', '8314', '8323', '8325', '8326', '8329', '8334', '8335', '8337', '8343', '8351', '8353', '8354', '8356', '8357', '8365', '8366', '8369', '8374', '8376', '8378', '8387', '8388', '8389', '8395', '8411', '8412', '8419', '8420', '8422', '8423', '8424', '8426', '8428', '8430', '9057', '9058']



def get_ecotype_id_info_dict(defaultValue=None):
    """
    Returns a dict containing tuples of (nativename, stockparent, latitude, longitude, country)
    """
    import dbutils
    conn = dbutils.connect_to_default_lookup(db='stock')
    cursor = conn.cursor ()
    #Retrieve the filenames
    print "Fetching data"
    ecotypeDict = {}

    sql_statment = """
SELECT DISTINCT e.id, e.nativename, e.stockparent, e.latitude, e.longitude, c.abbr
FROM stock.ecotype e, stock.site s, stock.address a, stock.country c
WHERE e.siteid = s.id AND s.addressid = a.id AND a.countryid = c.id
"""
    numRows = int(cursor.execute(sql_statment))
    i = 0
    while(1):
        row = cursor.fetchone()
        if not row:
            break;
        latitude = None
        longitude = None
        country = None
        if row[3]:
            latitude = float(row[3])
        if row[4]:
            longitude = float(row[4])
        if row[5]:
            country = row[5]
        ecotypeDict[int(row[0])] = (row[1], row[2], latitude, longitude, country)
        sp = "NA"
        if row[2]:
            sp = row[2]
        #print str(int(row[0]))+","+str(row[1])+","+sp
    cursor.close()
    conn.close()
    return ecotypeDict


def get_tg_ecotype_map():
    """
    Returns a dict containing tuples of (nativename, stockparent, latitude, longitude, country)
    """
    import dbutils
    conn = dbutils.connect_to_default_lookup(db='stock')
    cursor = conn.cursor ()
    #Retrieve the filenames
    print "Fetching data"
    tg_ets_map = {}

    sql_statment = """
SELECT DISTINCT e2e.ecotypeid, e2e.tg_ecotypeid, e2e.name, e2e.nativename, e2e.stockparent
FROM stock.ecotypeid2tg_ecotypeid e2e
"""
    numRows = int(cursor.execute(sql_statment))
    while(1):
        row = cursor.fetchone()
        if not row:
            break;
        tg_ets_map[int(row[0])] = (row[1], row[2], row[3], row[4])
    cursor.close()
    conn.close()
    return tg_ets_map




def get_250K_accession_to_ecotype_dict(call_method=72, dict_key='nativename'):
    """
    """
    delim = '\t'
    fn = env['data_dir'] + 'call_method_%d_info.tsv' % call_method
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



def get_accession_to_ecotype_id_dict(accessions, stockParents=None, defaultValue=None, lower_cases=True):
    import warnings
    warnings.warn("This function is possibly outdated, please update SQL statement before use!")
    import dbutils
    conn = dbutils.connect_to_default_lookup(db='stock')
    cursor = conn.cursor ()
    #Retrieve the filenames
    print "Fetching data in DB."
    accDict = {}
    not_found_count = 0
    for i in range(0, len(accessions)):
        acc = accessions[i]
        acc = unicode(acc, "latin-1")
        #print acc, unicode(acc, "latin-1")
        if stockParents != None and stockParents[i] != "NA":
            sp = stockParents[i]
            sql_statement = "SELECT DISTINCT ei.tg_ecotypeid, ei.nativename \
                     FROM stock.ecotypeid2tg_ecotypeid ei \
                     WHERE ei.nativename like '%s' and ei.stockparent like '%s' " % (acc.lower(), sp)
        else:
            sql_statement = "SELECT DISTINCT ei.tg_ecotypeid, ei.nativename \
                     FROM stock.ecotypeid2tg_ecotypeid ei \
                     WHERE ei.nativename like '%s' " % (acc.lower())
        #print sql_statement
        numRows = int(cursor.execute(sql_statement))
        i = 0
        accession_found = False
        while(1):
            row = cursor.fetchone()
            #print row
            if not row:
                if i == 0:
                    not_found_count += 1
                    print "Accession", acc, "wasn't found."
                break
            else:
                if i == 0:
                    ecotype = int(row[0])
                else:
                    ecotype = min(int(row[0]), ecotype)
                i += 1
        if i:
            if not accDict.has_key(acc.lower()):
                accDict[acc.lower()] = ecotype
                if not lower_cases:
                    accDict[acc.upper()] = ecotype
                    accDict[acc] = ecotype
            elif  accDict[acc.lower()] > ecotype:
                accDict[acc.lower()] = ecotype
                if not lower_cases:
                    accDict[acc.upper()] = ecotype
                    accDict[acc] = ecotype
                #print acc.lower(),":", row[1], int(row[0])
    cursor.close ()
    conn.close ()
    if not_found_count:
        print not_found_count, "accessions weren't found."
    else:
        print "All accessions were found in DB."
    return accDict






def insert_phenotypes_into_db(phen_file, method_description='', comment='', transformation_description='None',
            data_description='', growth_condition='', biology_category_id=4, citations='',
            phenotype_scoring='', method_id=None, data_type=''):
    """
    Insert phenotypes into the  DB. 
    """
    phend = parse_phenotype_file(phen_file)
    phend.convert_to_averages()
    phend.insert_into_db(phenotype_scoring=phenotype_scoring, method_description=method_description,
                 growth_condition=growth_condition, biology_category_id=biology_category_id,
                 citations=citations, data_description=data_description,
                 transformation_description=transformation_description, method_id=method_id,
                 data_type=data_type, comment=comment)




def _castric_plot_():
    fn = '/Users/bjarni.vilhjalmsson/Projects/vincent_castric_coordinates.txt'
    ecotypes = []
    with open(fn) as f:
        for l in f:
            ecotypes.append(l.split('.')[1].strip())
    dummy_phen = [1] * len(ecotypes)
    phen_dict = {1:{'name':'test', 'ecotypes':ecotypes, 'values':dummy_phen}}
    pd = phenotype_data(phen_dict)
    ecotypes, acc_names, lats, lons = pd.plot_accession_map(1, pdf_file=env['tmp_dir'] + 'castric_plain.pdf', png_file=env['tmp_dir'] + 'castric_plain.png',
            with_colorbar=False, map_type='europe')
    with open(env['tmp_dir'] + 'castric_plain.csv', 'w') as f:
        f.write('ecotype_id,accession_name,latitude,longitude\n')
        for t in it.izip(ecotypes, acc_names, lats, lons):
            f.write('%s,%s,%0.4f,%0.4f\n' % t)


def _test_bs_herits_():
    pd = parse_phenotype_file(env['phen_dir'] + 'seed_size.csv')
    pd.get_broad_sense_heritability()


#def _get_107_traits_():
#    ets_199 = get_199_ecotypes()
#    pids = categories_2_phenotypes[1] + categories_2_phenotypes[2] + categories_2_phenotypes[3] + categories_2_phenotypes[4]
#    #pids = [pid for pid in categories_2_phenotypes]
#    phed = parse_phenotype_file(env['phen_dir'] + 'phen_avg_112210.csv')
#    phed.filter_phenotypes(pids)
#    phed.filter_ecotypes_2(ets_199)
#    phed.write_to_file(env['phen_dir'] + '199_phenotypes.csv')
#    import dataParsers as dp
#    sd = dp.load_250K_snps(32)
#    sd.filter_accessions(ets_199)
#    print sd.accessions
#    sd.writeToFile(env['data_dir'] + '199_genotypes.csv')
#    print len(sd.accessions)


def plot_test():
    import dataParsers as dp
    k, k_accessions = dp.load_kinship(return_accessions=True, scaled=False)
    phed = get_all_phenotypes_from_db() #parse_phenotype_file(env['phen_dir'] + 'phen_raw_112210.csv')
    phed.plot_phen_relatedness(k, k_accessions, env['tmp_dir'] + 'rel_plots')



def _test_blups_(pid):
    import dataParsers as dp
    import linear_models as lm
    sd = dp.load_snps_call_method()
    phed = get_phenotypes_from_db([pid])
    phed.convert_to_averages()
    sd.coordinate_w_phenotype_data(phed, pid)
    K = lm.scale_k(sd.get_ibs_kinship_matrix())
    d = phed.get_blup(pid, K)
    c = sp.corrcoef(d['u_blup'].flatten(), d['values'].flatten())[0, 1] ** 2
    print c, d['pseudo_heritability']
    return d


if __name__ == '__main__':
#    insert_phenotypes_into_db(env['phen_dir'] + 'b_dilkes_metabolites.csv', method_description='Metabolites, mass spec?',
#            comment='Data from Brian Dilkes')
    #_castric_plot_()
    #_test_bs_herits_()
    _test_blups_(int(sys.argv[1]))
    print "Done!"
