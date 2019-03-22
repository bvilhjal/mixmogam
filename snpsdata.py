"""
This python library aims to do two things.
1. Offer general wrapper classes around SNPs datasets.
2. Offer basic functions which can aid analysis of the SNPs.

Author: Bjarni J. Vilhjalmsson
Email: bjarni.vilhjalmsson@gmail.com
"""

import sys, warnings
from itertools import *
import bisect
import h5py
import os
import random
import kinship
try:
    import scipy as sp
    sp.seterr(divide='raise')
except Exception, err_str:
    print 'scipy is missing:', err_str


IUPAC_alphabet = ['A', 'C', 'G', 'T', '-', 'Y', 'R', 'W', 'S', 'K', 'M', 'D', 'H', 'V', 'B', 'X', 'N']

"""
Marker type suggestions:
 'Unknown'
 'SNP'
 'Indel'
 'Other'
  ... 
  etc.
  
"""

def get_haplotypes(snps, num_accessions, count_haplotypes=False):
    curr_haplotypes = map(tuple, sp.transpose(sp.array(snps).tolist()))
    hap_set = list(set(curr_haplotypes))
    hap_hash = {}
    for j, h in enumerate(hap_set):
        hap_hash[h] = j
    new_haplotypes = sp.zeros(num_accessions, dtype='int8')
    for j, h in enumerate(curr_haplotypes):
        new_haplotypes[j] = hap_hash[h]
    if count_haplotypes:
        hap_counts = [curr_haplotypes.count(hap) for hap in hap_set]
        return hap_set, hap_counts, new_haplotypes
    else:
        return new_haplotypes


class _SnpsData_(object):
    """
    An abstract superclass.
    """
    def __init__(self, snps, positions, baseScale=None, accessions=None, arrayIds=None, chromosome=None,
            alignment_positions=None, id=None, marker_types=None, alphabet=None, associated_positions=None):
        self.snps = snps  # list[position_index][accession_index]
        self.positions = positions  # list[position_index]
        if accessions:
            self.accessions = accessions  # list[accession_index]
            # self._convert_to_tg_ecotypes_()
        if arrayIds:
            self.arrayIds = arrayIds  # list[accession_index]
        self.chromosome = chromosome
        self.alignment_positions = alignment_positions
        self.id = id
        self.marker_types = marker_types  # Where do these markers come frome, what type are they?  Useful for later analysis.
        self.alphabet = None
        self.missingVal = None
        self.associated_positions = associated_positions




    def merge_data(self, sd, union_accessions=True, allow_multiple_markers=True, error_threshold=0.2):
        """
        Merges data, allowing multiple markers at a position. (E.g. deletions and SNPs.)
        However it merges markers which overlap to a significant degree.
        """
        if union_accessions:
            new_accessions = list(set(self.accessions).union(set(sd.accessions)))
        else:
            new_accessions = self.accessions
        acc_map = []
        for acc in new_accessions:
            try:
                ai1 = self.accessions.index(acc)
            except:
                ai1 = -1
            try:
                ai2 = sd.accessions.index(acc)
            except:
                ai2 = -1
            acc_map.append((ai1, ai2))
        # print acc_map
        # pdb.set_trace()
        index_dict = {}
        j = 0
        last_pos = sd.positions[j]
        for i, pos in enumerate(self.positions):
            curr_pos = last_pos
            while j < len(sd.positions) and curr_pos < pos:
                j += 1
                if j < len(sd.positions):
                    curr_pos = sd.positions[j]
            last_pos = curr_pos
            index_list = []
            while j < len(sd.positions) and curr_pos == pos:
                index_list.append(j)
                j += 1
                if j < len(sd.positions):
                    curr_pos = sd.positions[j]
            if index_list:
                index_dict[i] = index_list


        indices_to_skip = set()
        new_snps = []
        merge_count = 0
        for i, snp1 in enumerate(self.snps):
            new_snp = []
            if i in index_dict:
                index_list = index_dict[i]
                for j in index_list:
                    error_count = 0
                    t_count = 0
                    snp2 = sd.snps[j]
                    for (ai1, ai2) in acc_map:
                        if ai1 != -1 and ai2 != -1:
                            if snp1[ai1] != snp2[ai2] and snp1[ai1] != self.missingVal\
                                        and snp2[ai2] != self.missingVal:
                                error_count += 1
                            t_count += 1
                    if t_count > 0 and error_count / float(t_count) < error_threshold:
                        # print "Merge error is %f"%(error_count/float(t_count))
                        for (ai1, ai2) in acc_map:
                            if ai1 != -1 and ai2 != -1:
                                if snp1[ai1] != self.missingVal:
                                    new_snp.append(snp1[ai1])
                                else:
                                    new_snp.append(snp2[ai2])
                            elif ai1 == -1:
                                new_snp.append(snp2[ai2])
                            else:
                                new_snp.append(snp1[ai1])
                        merge_count += 1
                        indices_to_skip.add(j)
                    elif t_count > 0:
                        print "Not merging since error is %f." % (error_count / float(t_count))
                        for (ai1, ai2) in acc_map:
                            if ai1 != -1:
                                new_snp.append(snp1[ai1])
                            else:
                                new_snp.append(self.missingVal)

            else:
                for (ai1, ai2) in acc_map:
                    if ai1 != -1:
                        new_snp.append(snp1[ai1])
                    else:
                        new_snp.append(self.missingVal)
                # print snp1,new_snp
            if new_snp == []:
                raise Exception
            new_snps.append(new_snp)
        new_positions = self.positions

        for j in range(len(sd.snps)):
            if not j in indices_to_skip:
                snp2 = sd.snps[j]
                new_snp = []
                for (ai1, ai2) in acc_map:
                    if ai2 != -1:
                        new_snp.append(snp2[ai2])
                    else:
                        new_snp.append(self.missingVal)
                if new_snp == []:
                    raise Exception
                new_snps.append(new_snp)
                new_positions.append(sd.positions[j])



        pos_snp_list = zip(new_positions, new_snps)
        pos_snp_list.sort()
        r = map(list, zip(*pos_snp_list))
        self.positions = r[0]
        self.snps = r[1]
        self.accessions = new_accessions
        if len(self.snps) != len(self.positions):
            raise Exception
        print "Merged %d SNPs!" % (merge_count)
        print "Resulting in %d SNPs in total" % len(self.snps)


    def extend(self, snpsd, union_accessions=True):
        """
        Extends the marker data with another snpsd, but assumes they are non overlapping 
        and spatially ordered.  
        """
        if union_accessions:
            new_accessions = list(set(self.accessions).union(set(snpsd.accessions)))
        else:
            raise NotImplementedError

        acc_map = []
        for acc in new_accessions:
            try:
                ai1 = self.accessions.index(acc)
            except:
                ai1 = -1
            try:
                ai2 = snpsd.accessions.index(acc)
            except:
                ai2 = -1
            acc_map.append((ai1, ai2))

        new_snps = []
        new_positions = []
        for i, snp in enumerate(self.snps):
            new_snp = []
            for (ai1, ai2) in acc_map:
                if ai1 == -1:
                    new_snp.append(self.missingVal)
                else:
                    new_snp.append(snp[ai1])
            new_positions.append(self.positions[i])

        for i, snp in enumerate(snpsd.snps):
            new_snp = []
            for (ai1, ai2) in acc_map:
                if ai2 == -1:
                    new_snp.append(self.missingVal)
                else:
                    new_snp.append(snp[ai2])
            new_positions.append(snpsd.positions[i])
        self.snps = new_snps
        self.positions = new_positions
        self.accessions = new_accessions


    def sample_snps(self, random_fraction, seed=None):
        """
        Discards SNPs, leaving a random fraction of them untouched.
        """
        if seed:
            sp.random.seed(seed)
        rs = sp.random.random(len(self.snps))
        new_positions = []
        new_snps = []
        for i in range(len(self.snps)):
            if rs[i] < random_fraction:
                new_positions.append(self.positions[i])
                new_snps.append(self.snps[i])
        self.positions, self.snps = new_positions, new_snps


    def scalePositions(self, baseScale):
        for i in range(0, len(self.positions)):
            self.positions[i] = int(self.positions[i] * baseScale)
        self.baseScale = baseScale


    def addSnp(self, snp):
        self.snps.append(snp)


    def addPos(self, position):
        self.positions.append(position)


    def removeAccessions(self, accessions, arrayIds=None):
        """
        Removes accessions from the data.
        """
        accPair = set()
        for i in range(0, len(accessions)):
            if arrayIds:
                key = (accessions[i], arrayIds[i])
            else:
                key = accessions[i]
            accPair.add(key)
        newAccessions = []
        newArrayIds = []
        for i in range(0, len(self.snps)):
            snp = self.snps[i]
            newSnp = []
            for j in range(0, len(self.accessions)):
                if arrayIds:
                    key = (self.accessions[j], self.arrayIds[j])
                else:
                    key = self.accessions[j]
                if not key in accPair:
                    newSnp.append(snp[j])
                    if i == 0:
                        newAccessions.append(self.accessions[j])
                        if arrayIds:
                            newArrayIds.append(self.arrayIds[j])
            self.snps[i] = newSnp
        self.accessions = newAccessions
        if arrayIds:
            self.arrayIds = newArrayIds




    def removeAccessionIndices(self, indicesToKeep):
        """
        Removes accessions from the data.
        """
        newAccessions = []
        newArrayIds = []
        for i in indicesToKeep:
            newAccessions.append(self.accessions[i])
            if self.arrayIds:
                newArrayIds.append(self.arrayIds[i])
        for i in range(len(self.snps)):
            snp = self.snps[i]
            newSnp = []
            for j in indicesToKeep:
                newSnp.append(snp[j])
            self.snps[i] = newSnp
        self.accessions = newAccessions
        if self.arrayIds:
            # print "removeAccessionIndices: has array IDs: self.arrayIds =",self.arrayIds
            self.arrayIds = newArrayIds
            # print "len(self.arrayIds):",len(self.arrayIds)
        # pdb.set_trace()
        # print "len(self.accessions):",len(self.accessions)


    def filterMonoMorphicSnps(self):
        """
        05/12/08 yh. add no_of_monomorphic_snps_removed
        Removes SNPs from the data which are monomorphic.
        """
        newPositions = []
        newSnps = []
        for i in range(0, len(self.positions)):
            count = 0
            for nt in self.alphabet:
                if nt in self.snps[i]:
                    count += 1
            if count > 1:
                newSnps.append(self.snps[i])
                newPositions.append(self.positions[i])
        numRemoved = len(self.positions) - len(newPositions)
        self.no_of_monomorphic_snps_removed = numRemoved
        self.snps = newSnps
        self.positions = newPositions
        return numRemoved

    def onlyBinarySnps(self):
        """
        Removes all but binary SNPs.  (I.e. monomorphic, tertiary and quaternary alleles SNPs are removed.)
        """
        newPositions = []
        newSnps = []
        for i in range(0, len(self.positions)):
            count = 0
            for nt in self.alphabet:
                if nt in self.snps[i]:
                    count += 1
            if count == 2:
                newSnps.append(self.snps[i])
                newPositions.append(self.positions[i])
        numRemoved = len(self.positions) - len(newPositions)
        self.no_of_nonbinary_snps_removed = numRemoved
        self.snps = newSnps
        self.positions = newPositions
        print "Removed %d non-binary SNPs, leaving %d SNPs in total." % (numRemoved, len(self.snps))

        return numRemoved



    def orderAccessions(self, accessionMapping):
        newAccessions = [None for i in range(0, len(accessionMapping))]
        for (i, j) in accessionMapping:
            newAccessions[j] = self.accessions[i]

        newSnps = []
        for pos in self.positions:
            newSnps.append([self.missingVal for i in range(0, len(accessionMapping))])

        for (i, j) in accessionMapping:
            for posIndex in range(0, len(self.positions)):
                newSnps[posIndex][j] = self.snps[posIndex][i]

        self.accessions = newAccessions
        self.snps = newSnps

        if self.arrayIds:
            newArrayIds = [None for i in range(0, len(self.arrayIds))]
            for (i, j) in accessionMapping:
                newArrayIds[j] = self.arrayIds[i]
            self.arrayIds = newArrayIds




    def filterRegion(self, startPos, endPos):
        """
        Filter all SNPs but those in region.
        """
        i = 0
        while i < len(self.snps) and self.positions[i] < startPos:
            i += 1

        newSnps = []
        newPositions = []
        while i < len(self.snps) and self.positions[i] < endPos:
            newSnps.append(self.snps[i])
            newPositions.append(self.positions[i])
            i += 1
        self.snps = newSnps
        self.positions = newPositions

    def get_region_snpsd(self, start_pos, end_pos):
        """
        Constructs and returns a new object based on the prev. using SNPs only in the given region.
        """
        import copy
        snpsd = copy.deepcopy(self)  # Clone
        i = 0
        while i < len(self.snps) and self.positions[i] < start_pos:
            i += 1

        new_snps = []
        new_positions = []
        while i < len(self.snps) and self.positions[i] < end_pos:
            new_snps.append(self.snps[i])
            new_positions.append(self.positions[i])
            i += 1
        snpsd.positions = new_positions
        snpsd.snps = new_snps
        return snpsd



    def get_snp_phen_pair(self, snp_index, phen_data, phen_index=None, missingVal=None):
        """
        returns a pair of snp and phenotype values, removing NA's
        """
        if not missingVal:
            missingVal = self.missingVal
        snp = self.snps[snp_index]
        phen_vals = phen_data.getPhenVals(phen_index, noNAs=False)
#        print phen_data.accessions
#        print self.accessions
#        print len(phen_vals),len(phen_data.accessions)
        acc_set = list(set(self.accessions).intersection(set(phen_data.accessions)))
        acc_set.sort()  # Not really necessary...
        new_snp = []
        new_phen_vals = []
        for acc in acc_set:
            si = self.accessions.index(acc)
            pi = phen_data.accessions.index(acc)
            if snp[si] != missingVal and phen_vals[pi] != missingVal:
                new_snp.append(snp[si])
                new_phen_vals.append(phen_vals[pi])
        return (new_snp, new_phen_vals)


    def get_mafs(self):
        """
        Returns MAFs and MARFs
        
        Assumes that this is a binary allele.
        """
        mafs = []
        marfs = []
        for snp in self.snps:
            missing_count = list(snp).count(self.missingVal)
            num_nts = len(snp) - missing_count
            nts = set(snp)
            if missing_count:
                nts.remove(self.missingVal)
            c = list(snp).count(nts.pop()) / float(num_nts)
            if c > 0.5:
                c = 1.0 - c
            marfs.append(c)
            mafs.append(int(c * num_nts))
        return {"mafs":mafs, "marfs":marfs}





    def remove_snps(self, remove_indices):
        """
        Filters all SNPs in the indices list.
        """
        remove_indices.sort(reverse=True)
        for i in remove_indices:
            del self.snps[i]
            del self.positions[i]
            if self.alignment_positions:
                del self.alignment_positions[i]


    def filter_snp_indices(self, indices):
        """
        Filters all SNPs but those in the indices list.
        """
        new_snps = []
        new_positions = []
        new_alignm_positions = []
        for i in indices:
            new_snps.append(self.snps[i])
            new_positions.append(self.positions[i])
        self.snps = new_snps
        self.positions = new_positions
        if self.alignment_positions:
            new_alignm_positions = []
            for i in indices:
                new_alignm_positions.append(self.alignment_positions[i])
            self.alignment_positions = new_alignm_positions


    def filter_accessions_by_NAs(self, maxMissing):
        """
        Removes the accessions in the snpsd if it has NA-rates greater than maxMissing. 
        """
        sys.stderr.write("Calculating NA rate ...")
        missingCounts = self.accessionsMissingCounts()

        indices_to_keep = []
        num_snps = float(len(self.snps))
        for i, mc in enumerate(missingCounts):
            if (mc / num_snps) <= maxMissing:
                indices_to_keep.append(i)
        self.removeAccessionIndices(indices_to_keep)


    def filter_na_snps(self, max_na_rate=0.0):
        new_snps = []
        new_positions = []
        num_snps = len(self.snps)
        e_num = float(len(self.accessions))
        for i, snp in enumerate(self.snps):
            if snp.count(self.missingVal) / e_num <= max_na_rate:
                new_snps.append(snp)
                new_positions.append(self.positions[i])
        self.snps = new_snps
        self.positions = new_positions
        print "Removed %d SNPs, leaving %d in total." % (num_snps - len(self.snps), len(self.snps))



    def accessionsMissingCounts(self):
        """
        Returns a list of accessions and their missing value rates.
        """
        missingCounts = [0] * len(self.accessions)

        for i in range(0, len(self.positions)):
            for j in range(0, len(self.accessions)):
                if self.snps[i][j] == self.missingVal:
                    missingCounts[j] += 1

        return missingCounts

    def __str__(self):
        st = "Number of accessions: " + str(len(self.accessions)) + "\n"
        st += "Number of SNPs: " + str(len(self.snps)) + "\n"
        uniqueAccessions = list(set(self.accessions))
        if len(uniqueAccessions) < len(self.accessions):
            st += "Not all accessions are unique. \n" + "Number of unique accessions: " + str(len(uniqueAccessions)) + "\n"
            for acc in uniqueAccessions:
                count = self.accessions.count(acc)
                if count > 1:
                    st += acc + " has " + str(count) + " occurrences.\n\n"
        return st


    def countMissingSnps(self):
        """
        Returns a list of accessions and their missing value rates.
        """
        missingCounts = [0.0 for acc in self.accessions]
        totalCounts = [0.0 for acc in self.accessions]
        for i in range(0, len(self.positions)):
            for j in range(0, len(self.accessions)):
                totalCounts[j] += 1
                if self.snps[i][j] == self.missingVal:
                    missingCounts[j] += 1

        return [missingCounts[i] / totalCounts[i] for i in range(len(self.accessions))], missingCounts


    def print_ecotype_info_to_file(self, filename):
        """
        Prints info on the accessions/ecotype IDs to a file.
        """
        import phenotypeData as pd
        eid = pd.get_ecotype_id_info_dict()
        f = open(filename, 'w')
        f.write((", ".join(["Ecotype_id", "Native_name", "Stock_parent", "latitude", "longitude", "country"])) + "\n")
        for e in self.accessions:
            l = [e] + list(eid[int(e)])
            l = map(str, l)
            f.write((", ".join(l)) + "\n")
        f.close()




class RawDecoder(dict):
    def __missing__(self, key):
        return 'NA'
    def __init__(self, initdict={}):
        for letter in ['A', 'C', 'G', 'T']:
            self[letter] = letter
        for letter in initdict:
            self[letter] = initdict[letter]


class RawSnpsData(_SnpsData_):
    """
    05/11/2008 yh. give default values to all initial arguments so that it could be called without any arguments.
    Similar to the SnpsData class, except it deals with bases (A,C,G,T), instead of allele numbers (0s and 1s).
    
    Alphabet: A,C,G,T, and NA

    """
    # alphabet = ['A','C','G','T']
    # alphabet = ['A','C','G','T','-']




    def __init__(self, snps=None, positions=None, baseScale=None, accessions=None, arrayIds=None,
            chromosome=None, callProbabilities=None, alignment_positions=None, id=None,
            marker_types=None, missing_val='NA'):
        self.snps = snps
        self.positions = positions
        self.accessions = accessions
        # if accessions:
        #    self._convert_to_tg_ecotypes_()

        self.arrayIds = arrayIds
        self.callProbabilities = []  # list[position_index][accession_index]
        if callProbabilities:
            self.callProbabilities = callProbabilities
        self.alignment_positions = alignment_positions
        self.missingVal = missing_val
        self.alphabet = IUPAC_alphabet
        self.marker_types = marker_types  # Where do these markers come frome, what type are they?  Useful for later analysis.
        self.id = id
        self.chromosome = chromosome


    def writeToFile(self, filename, chromosome):
        """
        Writes data to a file.  1 file for each chromosome.

        WARNING OLD, outdated!
        """
        outStr = ""
        fieldStrings = ["Chromosome", "Positions"]
        for acc in self.accessions:
            fieldStrings.append(str(acc))
        outStr += ", ".join(fieldStrings) + "\n"
        for i in range(0, len(self.positions)):
            outStr += str(chromosome) + ", " + str(self.positions[i]) + ", " + ", ".join(self.snps[i]) + "\n"
        f = open(filename, "w")
        f.write(outStr)
        f.close()



#    def impute_data_region(self, reference_snpsd, window_count=100, window_size=None, verbose=False, **kwargs):
#        """
#        Impute any NA SNPs.
#        """
#        import ImputeSNPs as imp
#        import tempfile, copy
#        if window_size:
#            start_pos = max(0, self.positions[0] - window_size)
#            end_pos = self.positions[-1] + window_size
#            snpsd = reference_snpsd.get_region_snpsd(start_pos, end_pos)
#        else:
#            snpsd = copy.deepcopy(reference_snpsd)
#        if (not window_count) and window_size:
#            i = 0
#            while i < len(snpsd.positions) and snpsd.positions[i] < self.positions[0]:
#                i += 1
#            suffix_count = i
#            i = 1
#            while i >= len(snpsd.positions) and snpsd.positions[len(snpsd.positions) - i] > self.positions[-1]:
#                i += 1
#            window_count = max(suffix_count + 1, i)
#
#        snpsd.mergeDataUnion(self, priority=2, unionType=3, verbose=verbose)
#        snpsd.na_ambigious_calls()
#        snpsd.onlyBinarySnps()
#        print snpsd.snps
#        print snpsd.accessions
#        print len(snpsd.accessions), len(set(snpsd.accessions))
#        tmpFile1 = tempfile.mkstemp()
#        os.close(tmpFile1[0])
#        tmpFile2 = tempfile.mkstemp()
#        os.close(tmpFile2[0])
#
#        if verbose:
#            print "Preparing data in", tmpFile1[1]
#        imp.writeAsNputeFile(snpsd, tmpFile1[1])
#        imp.checkNputeFile(tmpFile1[1])
#        nputeCmd = "python " + imp.path_NPUTE + "NPUTE.py -m 0 -w " + str(window_count) + " -i " + str(tmpFile1[1]) + " -o " + str(tmpFile2[1])
#        if verbose:
#            print "Imputing data..."
#            print nputeCmd
#        os.system(nputeCmd)
#        if verbose:
#            print "Imputation done!"
#
#        snpsd = imp.readNputeFile(tmpFile2[1], snpsd.accessions, snpsd.positions)
#        self.mergeDataUnion(snpsd, priority=2)
#        os.remove(tmpFile1[1])
#        os.remove(tmpFile2[1])
#
#
#    def impute_data(self, verbose=True):
#        """
#        Impute any NAs in the data.
#        """
#        import ImputeSNPs as imp
#        import tempfile
#        tmpFile1 = tempfile.mkstemp()
#        os.close(tmpFile1[0])
#        tmpFile2 = tempfile.mkstemp()
#        os.close(tmpFile2[0])
#
#        if verbose:
#            print "Preparing data in", tmpFile1[1]
#        imp.writeAsNputeFile(self, tmpFile1[1])
#        imp.checkNputeFile(tmpFile1[1])
#        nputeCmd = "python " + imp.path_NPUTE + "NPUTE.py -m 0 -w 30 -i " + str(tmpFile1[1]) + " -o " + str(tmpFile2[1])
#        if verbose:
#            print "Imputing data..."
#            print nputeCmd
#        os.system(nputeCmd)
#        if verbose:
#            print "Imputation done!"
# #        pdb.set_trace()
#        snpsd = imp.readNputeFile(tmpFile2[1], self.accessions, self.positions)
#        self.snps = snpsd.snps
#        os.remove(tmpFile1[1])
#        os.remove(tmpFile2[1])


    def compareWith(self, snpsd, withArrayIds=0, verbose=True, heterozygous2NA=False):
        """
        05/10/2008 len(commonSnpsPos) could be zero
        This function performs QC on two datasets.

        Requires accessions to be defined.

        withArrayIds = 0 (no array IDs), =1 the object called from has array IDs, =2 both objects have array IDs 
        """
        basicAlphabet = ['-', 'A', 'C', 'G', 'T']

        if verbose:
            print "Comparing datas"
            print "Number of snps:", len(self.snps), "and", len(snpsd.snps)
        # Find common accession indices
        accessionsIndices = []
        accessionErrorRate = []
        accessionCallRates = [[], []]
        accessionCounts = []
        commonAccessions = []
        arrayIds = []
        for i in range(0, len(self.accessions)):
            acc = self.accessions[i]
            if snpsd.accessions.count(acc):
                j = snpsd.accessions.index(acc)
                accessionsIndices.append([i, j])
                accessionErrorRate.append(0)
                accessionCounts.append(0)
                accessionCallRates[0].append(0)
                accessionCallRates[1].append(0)
                commonAccessions.append(acc)
                if withArrayIds > 0:
                    if withArrayIds == 1:
                        aId = self.arrayIds[i]
                    elif withArrayIds == 2:
                        aId = (self.arrayIds[i], snpsd.arrayIds[j])
                    arrayIds.append(aId)
        commonSnpsPos = []
        snpErrorRate = []
        snpCallRate = [[], []]
        goodSnpsCounts = []
        totalCounts = 0
        totalFails = 0
        i = 0
        j = 0
        while i <= len(self.positions) and j <= len(snpsd.positions):
            if i < len(self.positions):
                pos1 = self.positions[i]
            if j < len(snpsd.positions):
                pos2 = snpsd.positions[j]
            if i < len(self.positions) and pos1 < pos2:
                i = i + 1
            elif j < len(snpsd.positions) and pos2 < pos1:
                j = j + 1
            elif i < len(self.positions) and j < len(snpsd.positions) and pos1 == pos2:
                commonSnpsPos.append(pos1)
                fails = 0
                counts = 0
                missing1 = 0
                missing2 = 0
                for k in range(0, len(accessionsIndices)):
                    accIndices = accessionsIndices[k]
                    snp1 = self.snps[i][accIndices[0]]
                    snp2 = snpsd.snps[j][accIndices[1]]

                    if heterozygous2NA:
                        if snp1 in basicAlphabet and snp2 in basicAlphabet:
                            accessionCounts[k] += 1
                            counts += 1
                            if snp1 != snp2:
                                fails = fails + 1
                                accessionErrorRate[k] += 1
                        else:
                            if snp1 == self.missingVal:
                                accessionCallRates[0][k] += 1
                                missing1 += 1
                            if snp2 == self.missingVal:
                                accessionCallRates[1][k] += 1
                                missing2 += 1
                    else:
                        if snp1 != self.missingVal and snp2 != self.missingVal:
                            accessionCounts[k] += 1
                            counts += 1
                            if snp1 != snp2:
                                fails = fails + 1
                                accessionErrorRate[k] += 1
                        else:
                            if snp1 == self.missingVal:
                                accessionCallRates[0][k] += 1
                                missing1 += 1
                            if snp2 == self.missingVal:
                                accessionCallRates[1][k] += 1
                                missing2 += 1
                goodSnpsCounts.append(counts)
                error = 0
                totalCounts += counts
                totalFails += fails
                if counts > 0:
                    error = float(fails) / float(counts)
                    snpErrorRate.append(error)
                snpCallRate[0].append(missing1 / float(len(accessionsIndices)))
                snpCallRate[1].append(missing2 / float(len(accessionsIndices)))
                i = i + 1
                j = j + 1
            else:
                # One pointer has reached and the end and the other surpassed it.
                break

        for i in range(0, len(accessionErrorRate)):
            if accessionCounts[i] > 0:
                accessionErrorRate[i] = accessionErrorRate[i] / float(accessionCounts[i])
            no_of_common_snps_pos = len(commonSnpsPos)
            if no_of_common_snps_pos > 0:  # 05/08/10 yh
                accessionCallRates[0][i] = accessionCallRates[0][i] / float(no_of_common_snps_pos)
                accessionCallRates[1][i] = accessionCallRates[1][i] / float(no_of_common_snps_pos)
            else:
                accessionCallRates[0][i] = 0
                accessionCallRates[1][i] = 1

        print "In all", len(snpErrorRate), "common snps found"
        print "In all", len(commonAccessions), "common accessions found"
        print "Common accessions IDs :", commonAccessions
        print "Common SNPs positions :", commonSnpsPos
        print "Accessions error rates", accessionErrorRate
        print "Average Accession SNP Error:", sum(accessionErrorRate) / float(len(accessionErrorRate))
        print "SNP error rates", snpErrorRate
        print "Average Snp Error:", sum(snpErrorRate) / float(len(snpErrorRate))


        naCounts1 = [0] * len(accessionsIndices)
        for i in range(0, len(self.positions)):
            for k in range(0, len(accessionsIndices)):
                accIndices = accessionsIndices[k]
                snp = self.snps[i][accIndices[0]]
                if snp == self.missingVal:
                    naCounts1[k] += 1

        naCounts2 = [0] * len(accessionsIndices)
        for i in range(0, len(snpsd.positions)):
            for k in range(0, len(accessionsIndices)):
                accIndices = accessionsIndices[k]
                snp = snpsd.snps[i][accIndices[1]]
                if snp == self.missingVal:
                    naCounts2[k] += 1


        return [commonSnpsPos, snpErrorRate, commonAccessions, accessionErrorRate, accessionCallRates, arrayIds, accessionCounts, snpCallRate, [naCounts1, naCounts2], [totalCounts, totalFails]]



    def convert_to_binary(self):
        """
        An updated and version of the getSnpsData.. uses scipy.  
        """
        raise NotImplementedError



    def getSnpsData(self, missingVal= -1, reference_ecotype='6909', only_binary=True, verbose=True):
        """
        Returns a SnpsData object correspoding to this RawSnpsData object.

        Note that some of the SnpsData attributes are a shallow copy of the RawSnpsData obj.

        Reference ecotype is set to be the 0 allele. 
        """
        decoder = {self.missingVal:missingVal}  #Might cause errors somewhere???!!!
        coding_fun = sp.vectorize(lambda x: decoder[x], otypes=['int8'])

        if reference_ecotype in self.accessions:
            ref_i = self.accessions.index(reference_ecotype)
        else:
            ref_i = 0
            import warnings
            warnings.warn("Given reference ecotype %s wasn't found, using %s as 0-reference." % \
                    (reference_ecotype, self.accessions[ref_i]))
        snps = []
        positions = []
        num_lines = len(self.accessions)
        for snp_i, (snp, pos) in enumerate(izip(self.snps, self.positions)):
            if verbose and snp_i % 100000 == 0:
                print 'Converted %d SNPs.' % snp_i
            unique_nts = sp.unique(snp).tolist()
            if self.missingVal in unique_nts:
                if len(unique_nts) != 3:
                    continue  # Skipping non-binary SNP
                else:
                    unique_nts.remove(self.missingVal)
            else:
                if len(unique_nts) != 2:
                    continue  # Skipping non-binary SNP
            if snp[ref_i] != self.missingVal:
                col0_nt = snp[ref_i]
                decoder[col0_nt] = 0
                unique_nts.remove(col0_nt)
                decoder[unique_nts[0]] = 1
            else:
                decoder[unique_nts[0]] = 0
                decoder[unique_nts[1]] = 1
            snps.append(coding_fun(snp))
            positions.append(pos)

        print 'Removed %d non-binary SNPs out of %d, when converting to binary SNPs.'\
            % (len(self.positions) - len(positions), len(self.positions))

        assert len(snps) == len(positions), 'Somthing odd with the lengths.'
        return SNPsData(snps, positions, accessions=self.accessions, marker_types=self.marker_types, missing_val=missingVal)



    def filterBadSnps(self, snpsd, maxNumError=0):
        """
        05/12/08 add no_of_snps_filtered_by_mismatch
        Filters snps with high rate mismatches, when compared against another snpsd.
        """

        newSnps = []
        newPositions = []
        sys.stderr.write("Comparing datas. Number of snps: %s vs %s. \n" % (len(self.snps), len(snpsd.snps)))
        # Find common accession indices
        accessionsIndices = []
        commonAccessions = []
        for i in range(0, len(self.accessions)):
            acc = self.accessions[i]
            if snpsd.accessions.count(acc):
                j = snpsd.accessions.index(acc)
                accessionsIndices.append([i, j])
                commonAccessions.append(acc)


        commonSnpsPos = []
        i = 0
        j = 0
        while i <= len(self.positions) and j <= len(snpsd.positions):
            if i < len(self.positions):
                pos1 = self.positions[i]
            if j < len(snpsd.positions):
                pos2 = snpsd.positions[j]
            if i < len(self.positions) and pos1 < pos2:
                newSnps.append(self.snps[i])
                newPositions.append(self.positions[i])
                i = i + 1
            elif j < len(snpsd.positions) and pos2 < pos1:
                j = j + 1
            elif i < len(self.positions) and j < len(snpsd.positions) and pos1 == pos2:
                commonSnpsPos.append(pos1)
                fails = 0
                counts = 0
                for k in range(0, len(accessionsIndices)):
                    accIndices = accessionsIndices[k]
                    snp1 = self.snps[i][accIndices[0]]
                    snp2 = snpsd.snps[j][accIndices[1]]
                    if snp1 != self.missingVal and snp2 != self.missingVal:
                        counts += 1
                        if snp1 != snp2:
                            fails = fails + 1
                error = 0
                if counts > 0:
                    error = float(fails) / float(counts)
                if error <= maxNumError:
                    newSnps.append(self.snps[i])
                    newPositions.append(self.positions[i])
                i = i + 1
                j = j + 1
            else:
                # One pointer has reached and the end and the other surpassed it.
                break
        self.no_of_snps_filtered_by_mismatch = len(self.snps) - len(newSnps)
        sys.stderr.write('%s SNPs were filtered\n' % (self.no_of_snps_filtered_by_mismatch))
        self.snps = newSnps
        self.positions = newPositions


    def convertBadCallsToNA(self, minCallProb=0):
        """
        Converts all base calls with call prob. lower than the given one to NAs.
        """
        totalCount = len(self.positions) * len(self.snps[0])
        count = 0
        for i in range(0, len(self.positions)):
            for j in range(0, len(self.snps[i])):
                if self.callProbabilities[i][j] < minCallProb:
                    self.snps[i][j] = self.missingVal
                    count += 1
        fractionConverted = count / float(totalCount)
        return fractionConverted





    def getSnpsNArates(self):
        """
        Returns a list of accessions and their missing value rates.
        """
        naRates = []
        for i in range(0, len(self.positions)):
            missingCounts = 0
            totalCounts = 0
            for j in range(0, len(self.accessions)):
                totalCounts += 1
                if self.snps[i][j] == self.missingVal:
                    missingCounts += 1
            naRates.append(missingCounts / float(totalCounts))

        return naRates


    def filterMissingSnps(self, maxNumMissing=0):
        """
        05/12/2008 add no_of_snps_filtered_by_na
        Removes SNPs from the data which have more than maxNumMissing missing values.
        """
        warnings.warn("This function is deprecated!!")
        newPositions = []
        newSnps = []
        for i in range(0, len(self.positions)):
            missingCount = 0
            for nt in self.snps[i]:
                if nt == self.missingVal:
                    missingCount += 1
            if missingCount <= maxNumMissing:
                newSnps.append(self.snps[i])
                newPositions.append(self.positions[i])
        numRemoved = len(self.positions) - len(newPositions)
        self.no_of_snps_filtered_by_na = numRemoved
        self.snps = newSnps
        self.positions = newPositions
        return numRemoved


    def na_ambigious_calls(self, maxNumMissing=0):
        """
        NAs any ambiguous calls, i.e. anything but A,C,G,T,-
        """
        for snp in self.snps:
            for i, nt in enumerate(snp):
                if not nt in ['A', 'C', 'G', 'T', '-', self.missingVal]:
                    snp[i] = self.missingVal


    def filterMinMAF(self, minMAF=0):
        """
        Removes SNPs from the data which have a minor allele count less than the given one.
        """
        warnings.warn("This function is old, and should be replaced.. please don't use it!!")
        newPositions = []
        newSnps = []
        ntCounts = [0.0] * len(self.alphabet)
        ntl = self.alphabet
        for i in range(0, len(self.positions)):
            snp = self.snps[i]
            totalCount = 0
            for j in range(0, len(self.alphabet)):
                nt = ntl[j]
                c = snp.count(nt)
                ntCounts[j] = c
                totalCount += c
            if totalCount == 0 :
                print "snp:", snp
                print "self.alphabet:", self.alphabet
            maxAF = max(ntCounts)
            maf = 0
            if ntCounts.count(maxAF) < 2:
                for j in range(0, len(self.alphabet)):
                    if ntCounts[j] < maxAF and ntCounts[j] > 0:
                        if ntCounts[j] > maf:
                            maf = ntCounts[j]
            else:
                maf = maxAF

            if minMAF <= maf:
                newSnps.append(self.snps[i])
                newPositions.append(self.positions[i])

        numRemoved = len(self.positions) - len(newPositions)
        self.snps = newSnps
        self.positions = newPositions
        return numRemoved



    def filterMinMARF(self, minMARF=0):
        """
        Removes SNPs from the data which have a Relative MAF less than the given one.
        """
        warnings.warn("This function is old, and should be replaced.. please don't use it!!")
        newPositions = []
        newSnps = []
        ntCounts = [0.0] * len(self.alphabet)
        ntl = self.alphabet
        for i in range(0, len(self.positions)):
            snp = self.snps[i]
            totalCount = 0
            for j in range(0, len(self.alphabet)):
                nt = ntl[j]
                c = snp.count(nt)
                ntCounts[j] = c
                totalCount += c
            if totalCount == 0 :
                print "snp:", snp
                print "self.alphabet:", self.alphabet
            for j in range(0, len(self.alphabet)):
                ntCounts[j] = ntCounts[j] / float(totalCount)
            maxAF = max(ntCounts)
            maf = 0
            if ntCounts.count(maxAF) < 2:
                for j in range(0, len(self.alphabet)):
                    if ntCounts[j] < maxAF and ntCounts[j] > 0.0:
                        if ntCounts[j] > maf:
                            maf = ntCounts[j]
            else:
                maf = maxAF

            if minMARF <= maf:
                newSnps.append(self.snps[i])
                newPositions.append(self.positions[i])

        numRemoved = len(self.positions) - len(newPositions)
        self.snps = newSnps
        self.positions = newPositions
        return numRemoved



    def snpsFilterMAF(self, maf_thresholds):
        self.filterMinMARF(maf_thresholds[0])
        if maf_thresholds[1] < 0.5:
            raise Exception





    def mergeIdenticalAccessions(self, accessionIndexList, priority):
        """
        The priority argument gives particular accessions in the list priority
        if priority is set to 0, then majority (if any) rules.
        """
        pass


class SNPsData(_SnpsData_):
    """
    An alternative to the old SnpsData class, where this uses scipy to speed things up when possible.
    
    """
    alphabet = [-1, 0, 1, 2, 3]  # Here -1 is thought to be missing data value.
    def __init__(self, snps, positions, accessions=None, arrayIds=None, chromosome=None,
            alignment_positions=None, id=None, marker_types=None, missing_val= -1):
        self.snps = snps
        self.positions = positions
        self.accessions = accessions
        self.arrayIds = arrayIds
        self.chromosome = chromosome
        self.alignment_positions = alignment_positions
        self.marker_types = marker_types  # Where do these markers come frome, what type are they?  Useful for later analysis.
        self.id = id
        self.missingVal = missing_val



    def removeAccessionIndices(self, indicesToKeep):
        """
        Removes accessions from the data.
        """
        newAccessions = []
        newArrayIds = []
        for i in indicesToKeep:
            newAccessions.append(self.accessions[i])
            if self.arrayIds:
                newArrayIds.append(self.arrayIds[i])
        num_accessions = len(indicesToKeep)
        for i in range(len(self.snps)):
            self.snps[i] = self.snps[i][indicesToKeep]
#            snp = self.snps[i]
#            newSnp = sp.empty(num_accessions,dtype='int8')
#            for j,k in enumerate(indicesToKeep):
#                newSnp[j] = snp[k]
#            self.snps[i] = newSnp
        self.accessions = newAccessions
        if self.arrayIds:
            # print "removeAccessionIndices: has array IDs: self.arrayIds =",self.arrayIds
            self.arrayIds = newArrayIds
            # print "len(self.arrayIds):",len(self.arrayIds)
        # pdb.set_trace()
        # print "len(self.accessions):",len(self.accessions)


    def onlyBinarySnps(self):
        """
        Removes all but binary SNPs.  (I.e. monomorphic, tertiary and quaternary alleles SNPs are removed.)
        """
        new_positions = []
        new_snps = []
        for i, (snp, pos) in enumerate(izip(self.snps, self.positions)):
            if len(sp.unique(snp)) == 2:
                new_snps.append(snp)
                new_positions.append(pos)
        num_removed = len(self.positions) - len(new_positions)
        self.no_of_nonbinary_snps_removed = num_removed
        self.snps = new_snps
        self.positions = new_positions
        # print "Removed %d non-binary SNPs, leaving %d SNPs in total." % (num_removed, len(self.snps))
        return num_removed


    def remove_monomorphic_snps(self):
        """
        Removes all fixed SNPs.  (I.e. monomorphic.)
        """
        new_positions = []
        new_snps = []
        for i, (snp, pos) in enumerate(izip(self.snps, self.positions)):
            if len(sp.unique(snp)) > 1:
                new_snps.append(snp)
                new_positions.append(pos)
        num_removed = len(self.positions) - len(new_positions)
        self.snps = new_snps
        self.positions = new_positions
        # print "Removed %d non-binary SNPs, leaving %d SNPs in total." % (num_removed, len(self.snps))
        return num_removed


    def haplotize(self, snp_window=None, base_window=None):

        if snp_window:

            haplotypes_list = []
            i = 0
            num_accessions = len(self.accessions)
            while i < snp_window:
                cur_snps = self.snps[:i + snp_window + 1]
                haplotypes_list.append(get_haplotypes(cur_snps, num_accessions))
                i += 1
            while i < len(self.snps) - snp_window:
                cur_snps = self.snps[i - snp_window:i + snp_window + 1]
                haplotypes_list.append(get_haplotypes(cur_snps, num_accessions))
                i += 1
                if i % 100 == 0: print i
            while i < len(self.snps):
                cur_snps = self.snps[i - snp_window:]
                haplotypes_list.append(get_haplotypes(cur_snps, num_accessions))
                i += 1
            self.snps = haplotypes_list

        elif base_window:
            raise NotImplementedError
        else:
            raise Exception('Window size missing')


    def get_mafs(self, w_missing=False, type='binary'):
        """
        Returns MAFs and MARFs
        
        (Uses numpy.bincount)
        
        types supported: classes (e.g. binary data), diploid_ints, ..
        """


        def _get_maf_(snp):  # For missing data coded as -1s
            l = sp.bincount(snp)
            maf = max(l)
            for m in l[1:]:
                if m != 0 and m < maf:
                    maf = m
            return maf

        mafs = []
        marfs = []
        num_nts = len(self.snps[0])
        if w_missing:
            for snp in self.snps:
                if self.missingVal in snp:
                    missing_count = list(snp).count(self.missingVal)
                    num_nts = len(snp) - missing_count
                    nts = set(snp)
                    nts.remove(self.missingVal)
                    c = list(snp).count(nts.pop()) / float(num_nts)
                    if c > 0.5:
                        c = 1.0 - c
                    marfs.append(c)
                    mafs.append(int(c * num_nts))
                else:
                    l = sp.bincount(snp)
                    maf = min(l)
                    mafs.append(maf)
                    marfs.append(maf / float(num_nts))
        else:
            if type in ['binary', 'int']:
                for snp in self.snps:
                    l = sp.bincount(snp)
                    maf = min(l)
                    mafs.append(maf)
                    marfs.append(maf / float(num_nts))
            elif type == 'diploid_int':
                for snp in self.snps:
                    bin_counts = sp.bincount(snp)
                    bin_counts = sp.bincount(snp, minlength=3)
                    l = sp.array([bin_counts[0], bin_counts[2]]) + bin_counts[1] / 2.0
                    maf = l.min()
                    mafs.append(maf)
                    marfs.append(maf / float(num_nts))
            else:
                raise NotImplementedError

        return {"mafs":mafs, "marfs":marfs}




    def filter_mac(self, min_mac=15, w_missing=False, data_format='binary'):
        """
           Filter minor allele count SNPs.
        """
        print 'Filtering SNPs with MAC<%d, assuming %s data format' % (min_mac, data_format)
        new_snps = []
        new_positions = []
        if w_missing:
            raise NotImplementedError
        if data_format in ['binary', 'int']:
            for snp, pos in izip(self.snps, self.positions):
                bc = sp.bincount(snp)
                if len(bc) > 1 and bc.min() >= min_mac:
                    new_snps.append(snp)
                    new_positions.append(pos)
        elif data_format == 'diploid_int':
            for snp, pos in izip(self.snps, self.positions):
                bin_counts = sp.bincount(snp, minlength=3)
                l = sp.array([bin_counts[0], bin_counts[2]]) + bin_counts[1] / 2.0
                if l.min() >= min_mac:
                    new_snps.append(snp)
                    new_positions.append(pos)


        print 'Removed %d SNPs out of %d, leaving %d SNPs.' % (len(self.positions) - len(new_positions),
                                len(self.positions), len(new_positions))
        self.positions = new_positions
        self.snps = new_snps




    def merge_data(self, sd, acc_merge_type='intersection', error_threshold=0.1, discard_error_threshold=0.1):
        """
        Merges data, possibly allowing multiple markers at a position. (E.g. deletions and SNPs.)
        However it merges markers which overlap to a significant degree (error_threshold).
        
        (Uses scipy SNPs)
        """
        perc_overlap = len(set(self.accessions).intersection(set(sd.accessions))) \
                / float(len(set(self.accessions).union(set(sd.accessions))))
        print "Percentage of overlapping accessions %s" % perc_overlap
        if acc_merge_type == 'union':
            new_accessions = list(set(self.accessions).union(set(sd.accessions)))
        elif acc_merge_type == 'intersection':
            new_accessions = list(set(self.accessions).intersection(set(sd.accessions)))
        else:
            new_accessions = self.accessions
        acc_map = []
        for acc in new_accessions:
            try:
                ai1 = self.accessions.index(acc)
            except:
                ai1 = -1
            try:
                ai2 = sd.accessions.index(acc)
            except:
                ai2 = -1
            acc_map.append((ai1, ai2))

        num_accessions = len(new_accessions)
        missing_val_snp = sp.array(sp.repeat(self.missingVal, num_accessions), dtype='int8')

        # To handle multiple markers at the same position
        index_dict = {}
        j = 0
        last_pos = sd.positions[j]
        for i, pos in enumerate(self.positions):
            curr_pos = last_pos
            while j < len(sd.positions) and curr_pos < pos:
                j += 1
                if j < len(sd.positions):
                    curr_pos = sd.positions[j]
            last_pos = curr_pos
            index_list = []
            while j < len(sd.positions) and curr_pos == pos:
                index_list.append(j)
                j += 1
                if j < len(sd.positions):
                    curr_pos = sd.positions[j]
            if index_list:
                index_dict[i] = index_list


        indices_to_skip = set()  # Markers which are merged in the second SNPsData
        new_snps = []
        new_positions = []
        merge_count = 0
        for i, snp1 in enumerate(self.snps):
            new_snp = sp.array(sp.repeat(self.missingVal, num_accessions), dtype='int8')
            if i in index_dict:  # If there are markers at the same position.
                index_list = index_dict[i]
                for j in index_list:
                    error_count = 0
                    t_count = 0  # total count
                    snp2 = sd.snps[j]
                    for (ai1, ai2) in acc_map:
                        if ai1 != -1 and ai2 != -1:
                            if snp1[ai1] != snp2[ai2] and snp1[ai1] != self.missingVal\
                                        and snp2[ai2] != self.missingVal:
                                error_count += 1
                            t_count += 1
                    merge_error = error_count / float(t_count) if t_count > 0 else 1
                    if merge_error <= error_threshold or merge_error > discard_error_threshold:
                        indices_to_skip.add(j)

                        if merge_error <= error_threshold:
                            print "Merge error is %f" % merge_error
                            for ni, (ai1, ai2) in enumerate(acc_map):
                                if ai1 != -1 and ai2 != -1:
                                    if snp1[ai1] != self.missingVal:
                                        new_snp[ni] = snp1[ai1]
                                    else:
                                        new_snp[ni] = snp2[ai2]
                                elif ai1 == -1:
                                    new_snp[ni] = snp2[ai2]
                                else:
                                    new_snp[ni] = snp1[ai1]
                            merge_count += 1
                        else:
                            print 'Removing SNP at position %d since they have an error of %f'\
                                 % (self.positions[i], merge_error)

                    elif t_count > 0:
                        print "Not merging since error is %f" % merge_error
                        for ni, (ai1, ai2) in enumerate(acc_map):
                            if ai1 != -1:
                                new_snp[ni] = snp1[ai1]

            else:  # There were no markers at this position in the other snps data.
                for ni, (ai1, ai2) in enumerate(acc_map):
                    if ai1 != -1:
                        new_snp[ni] = snp1[ai1]


            if sp.any(new_snp != missing_val_snp):  # Some are not are missing
                new_snps.append(new_snp)
                new_positions.append(self.positions[i])

        print 'Inserting %d non-overlapping SNPs into the self snps data.' % len(sd.snps)
        for j in range(len(sd.snps)):
            if not j in indices_to_skip:  # There were no markers at this position in the other snps data.
                snp2 = sd.snps[j]
                new_snp = sp.array(sp.repeat(self.missingVal, num_accessions), dtype='int8')
                for ni, (ai1, ai2) in enumerate(acc_map):
                    if ai2 != -1:
                        new_snp[ni] = snp2[ai2]
                if new_snp == []:
                    raise Exception
                new_snps.append(new_snp)
                new_positions.append(sd.positions[j])


        print 'Sorting SNPs by positions..'
        pos_snp_list = zip(new_positions, range(len(new_positions)))
        pos_snp_list.sort()
        r = map(list, zip(*pos_snp_list))
        self.positions = r[0]
        self.snps = [new_snps[i] for i in r[1]]
        self.accessions = new_accessions
        if len(self.snps) != len(self.positions):
            raise Exception
        print "Merged %d SNPs!" % (merge_count)
        print "Resulting in %d SNPs in total" % len(self.snps)




class SnpsData(_SnpsData_):
    """
    A class for SNPs data.  It uses 0, 1 (and 2, 3 if there are more than 2 alleles) to represent SNPs.  
    -1 is used if the allele data is missing.

    This class should really be named marker data, as it is used for more than just SNP markers.

    Contains various functions that aid the analysis of the data.
    """
    # freqs = [] #list[position_index1][position_index2-position_index1+1] Linkage frequencies. 
    # baseScale = 1 #Scaling for positions
    alphabet = [0, 1, 2, 3]
    def __init__(self, snps, positions, baseScale=None, accessions=None, arrayIds=None, chromosome=None,
            alignment_positions=None, id=None, marker_types=None, missing_val= -1):
        self.snps = snps
        self.positions = positions
        if baseScale:
            self.scalePositions(baseScale)
        self.accessions = accessions
        self.arrayIds = arrayIds
        self.chromosome = chromosome
        self.alignment_positions = alignment_positions
        self.marker_types = marker_types  # Where do these markers come frome, what type are they?  Useful for later analysis.
        self.id = id
        self.missingVal = missing_val


    def clone(self):
        """
        Filter all SNPs but those in region.
        """
        newSNPs = []
        newPositions = []
        newAccessions = []
        for acc in self.accessions:
            newAccessions.append(acc)


        for i in range(0, len(self.positions)):
            new_snp = []
            for j in range(0, len(self.accessions)):
                new_snp.append(self.snps[i][j])
            newSNPs.append(new_snp)
            newPositions.append(self.positions[i])

        return SnpsData(newSNPs, newPositions, accessions=newAccessions, arrayIds=self.arrayIds, chromosome=self.chromosome)


    def remove_redundant_snps(self, r2_threshold=0.8, w_missing=False):
        """
        Removes any redundant SNPs as measured by correlation.  Returns information on what SNPs are correlated.
        """
        r2s = self.calc_r2_matrix(w_missing=w_missing)
        snps_indices = range(len(self.snps))
        new_snps_indices = []
        redundant_snps_indices_list = []
        redundant_pos_list = []
        while snps_indices:
            si1 = snps_indices.pop()
            new_snps_indices.append(si1)
            redundant_snps_indices = []
            redundant_positions = []
            for si2 in snps_indices[:]:
                if r2s[si1, si2] >= r2_threshold:
                    redundant_snps_indices.append(si2)
                    redundant_positions.append(self.positions[si2])
                    snps_indices.remove(si2)

            redundant_pos_list.append(redundant_positions)
            redundant_snps_indices_list.append(redundant_snps_indices)
        print "In total", len(self.snps), ", thereof", len(self.snps) - len(new_snps_indices), "are redundant."
        self.filter_snp_indices(new_snps_indices)
        self.redundant_positions = redundant_pos_list
        return redundant_pos_list



    def calcFreqs(self, windowSize, innerWindowSize=0):
        """
        Returns a list of two loci comparison frequencies with in a window.
        """
        freqs = []
        delta = 0
        if len(self.snps) > 0:
            delta = 1.0 / float(len(self.snps[0]))
        for i in xrange(0, len(self.snps) - 1):
            l = []
            j = i + 1
            while j < len(self.snps) and (self.positions[j] - self.positions[i]) < innerWindowSize :
                j = j + 1
            while j < len(self.snps) and (self.positions[j] - self.positions[i]) <= windowSize:
                jfreqs = [0.0] * 4  # The frequencies, of a pair of alleles.
                snp1 = self.snps[i]
                snp2 = self.snps[j]
                count = 0
                for k in xrange(0, len(snp1)):
                    val = snp1[k] * 2 + snp2[k]
                    jfreqs[val] = jfreqs[val] + delta
                l.append(jfreqs)
                j = j + 1
            freqs.append(l)
        self.freqs = freqs
        return self.freqs

    def calc_all_freqs(self, w_missing=False):
        """
        Returns a matrix of two loci comparison frequencies.
        """
        if w_missing:
            freqs = []
            for i, snp1 in enumerate(self.snps):
                l = []
                for j in range(i):
                    jfreqs = [0.0] * 4  # The frequencies, of a pair of alleles.
                    snp2 = self.snps[j]
                    count = 0.0
                    for nt1, nt2 in zip(snp1, snp2):
                        if nt1 != self.missingVal and nt2 != self.missingVal:
                            count += 1.0
                            val = nt1 * 2 + nt2
                            jfreqs[val] += 1.0
                    jfreqs = [f / count for f in jfreqs]
                    l.append(jfreqs)
                freqs.append(l)
            return freqs
        else:
            freqs = []
            delta = 0
            if len(self.snps) > 0:
                delta = 1.0 / float(len(self.snps[0]))
            for i, snp1 in enumerate(self.snps):
                l = []
                for j in range(i):
                    jfreqs = [0.0] * 4  # The frequencies, of a pair of alleles.
                    snp2 = self.snps[j]
                    for k in range(len(snp1)):
                        val = snp1[k] * 2 + snp2[k]
                        jfreqs[val] = jfreqs[val] + delta
                    l.append(jfreqs)
                freqs.append(l)

            return freqs

    def calcFreqsUnbiased(self, windowSize, innerWindowSize=0):  # Returns a list of two loci comparison frequencies with in a window.
        """ 
        Uses a distribution of frequencies that is not dependent on the lenght of the sequence.  (Alot of data is disregarded.)
        """
        numPairs = 0  # Counts the number of pairs
        freqs = []
        delta = 0
        if len(self.snps) > 0:
            delta = 1.0 / float(len(self.snps[0]))
        for i in xrange(0, len(self.snps) - 1):
            if (self.positions[len(self.snps) - 1] - self.positions[i]) >= windowSize:
                j = i + 1
                l = []
                while j < len(self.snps) and (self.positions[j] - self.positions[i]) < innerWindowSize :
                    j = j + 1
                while  j < len(self.snps) and (self.positions[j] - self.positions[i]) <= windowSize:
                    jfreqs = [0.0] * 4
                    snp1 = self.snps[i]
                    snp2 = self.snps[j]
                    count = 0
                    for k in xrange(0, len(snp1)):
                        val = snp1[k] * 2 + snp2[k]
                        jfreqs[val] = jfreqs[val] + delta
                    l.append(jfreqs)
                    j = j + 1
                    numPairs = numPairs + 1
                freqs.append(l)
            else:
                break
        self.freqs = freqs
        return self.freqs
        # return numPairs


    def calcFreqsSimple(self):
        freqs = []
        delta = 0
        if len(self.snps) > 0:
            delta = 1.0 / float(len(self.snps[0]))
        for i in xrange(0, len(self.snps) - 1):
            l = []
            for j in xrange(i + 1, len(self.snps)):
                jfreqs = [0.0] * 4
                snp1 = self.snps[i]
                snp2 = self.snps[j]
                count = 0
                for k in xrange(0, len(snp1)):
                    val = snp1[k] * 2 + snp2[k]
                    jfreqs[val] = jfreqs[val] + delta
                l.append(jfreqs)
            freqs.append(l)
        self.freqs = freqs
        return self.freqs


    def snpsFilter(self):
        """
        Splits the SNPs up after their allele frequency ..
        """
        snpsDatas = [SnpsData([], []), SnpsData([], []), SnpsData([], []), SnpsData([], []), SnpsData([], []), SnpsData([], []), SnpsData([], []), SnpsData([], []), SnpsData([], []), SnpsData([], [])]

        l = len(self.snps[0])
        for j in xrange(0, len(self.snps)):
            snp = self.snps[j]
            c = snp.count(0) / float(l)
            """
            if c>0.5:
            c = 1-c
            if c ==    0.5:
            c =0.499
            """
            if c == 1:
                c = 0.99999
            i = int(c * 10)
            snpsDatas[i].addSnp(snp)
            snpsDatas[i].addPos(self.positions[j])
        return snpsDatas


    def snpsFilterMAF(self, mafs, type='classes'):
        """
        Filters all snps with MAF not in the interval out of dataset.
        """
        newsnps = []
        newpos = []
        # print self.snps
        l = len(self.snps[0])
        if l == 0:
            print self.snps
        for j in xrange(0, len(self.snps)):
            snp = self.snps[j]
            c = snp.count(0) / float(l)
            if c > 0.5:
                c = 1 - c
            if c > mafs[0] and c <= mafs[1]:
                newsnps.append(snp)
                newpos.append(self.positions[j])
        if len(newsnps) == 0:
            print "Filtered out all snps from", len(self.snps), " what to do what to do?"
        del self.snps
        self.snps = newsnps
        del self.positions
        self.positions = newpos


    def snpsFilterRare(self, threshold=0.1):
        """Filters all snps with MAF of less than threshold out of dataset."""
        newsnps = []
        newpos = []
        # print self.snps
        l = len(self.snps[0])
        if l == 0:
            print self.snps
        for j in xrange(0, len(self.snps)):
            snp = self.snps[j]
            c = snp.count(0) / float(l)
            if c > 0.5:
                c = 1 - c
            if c > threshold:
                newsnps.append(snp)
                newpos.append(self.positions[j])
        # if len(newsnps)==0:  
            # print "Filtered out all snps from",len(self.snps)," what to do what to do?"
        del self.snps
        self.snps = newsnps
        del self.positions
        self.positions = newpos


    def _genRecombFile(self, filename, windowSize, maxNumPairs):
        n = len(self.snps[0])  # number of individuals/accessions
        numPairs = self.calcFreqsUnbiased(windowSize)
        if numPairs != 0:
            filterProb = float(maxNumPairs) / numPairs
        else:
            return numPairs
        f = open(filename, 'w')
        numPairs = 0
        for i in range(0, len(self.freqs)):
            for j in range(0, len(self.freqs[i])):
                if random.random() <= filterProb:
                    numPairs = numPairs + 1
                    st = str(i + 1) + " " + str(j + i + 2) + " " + str(self.positions[j + i + 1] - self.positions[i]) + " u "  # u denotes unknown as opposed to ad ancestral derived
                    st = st + str(int(self.freqs[i][j][0] * n + 0.5)) + " " + str(int(self.freqs[i][j][1] * n + 0.5)) + " "
                    st = st + str(int(self.freqs[i][j][2] * n + 0.5)) + " " + str(int(self.freqs[i][j][3] * n + 0.5)) + " 0 0 0 0 " + str(100 - n) + "\n"
                    f.write(st)
        f.close()
        return numPairs

#    def estimateRecomb(self, windowSize, maxNumPairs=10000, tempfile1="tmp1", tempfile2="tmp2", meanTract=200, numPoints=50):
#        num = self._genRecombFile(tempfile1, windowSize, maxNumPairs)
#        if num < 1:
#            return [0, 0, 0]
#        os.system(homedir + "Projects/programs/Maxhap/maxhap 1 " + homedir + "Projects/programs/Maxhap/h100rho  .0008 10 .1 0.01 500 " + str(numPoints) + " " + str(meanTract) + " < " + tempfile1 + " > " + tempfile2)
#        f = open(tempfile2, 'r')
#        lines = f.readlines()
#        npairs = float(lines[1].split()[2])
#        i = 2
#        while(lines[i].split()[0] == "Warning:"):
#            i = i + 1
#        rho = float(lines[i].split()[1])
#        ratio = float(lines[i].split()[2])
#        f.close()
#        return [rho, ratio, npairs]


    def meanAF(self):
        """ Mean allele frequency. """
        if len(self.snps):
            l = float(len(self.snps[0]))
            c = 0
            for j in xrange(0, len(self.snps)):
                snp = self.snps[j]
                snpsc = snp.count(0)
                if snpsc < (l / 2.0):
                    c = c + snpsc / l
                else:
                    c = c + abs((l / 2.0) - snpsc) / l
            return c / len(self.snps)
        else:
            return 0

    def EHH(self, snp1, snp2):
        """ Calculates the EHH between two SNPs"""
#        data = self.snps[snp1:snp2]
        haplotypes = []
        haplotypecount = []
        for i in range(0, len(self.snps[0])):
            haplotype = []
            for j in range(snp1, snp2 + 1):
                haplotype.append(self.snps[j][i])
            if not haplotype in haplotypes:
                haplotypes.append(haplotype)
                haplotypecount.append(1.0)
            else:
                k = haplotypes.index(haplotype)
                haplotypecount[k] = haplotypecount[k] + 1.0
        s = 0.0
        for i in range(0, len(haplotypes)):
            if haplotypecount[i] > 1:
                s = s + haplotypecount[i] * (haplotypecount[i] - 1)
        s = s / (len(self.snps[0]) * (len(self.snps[0]) - 1))
        return s

    def totalEHH(self, windowSize, innerWindowSize):
        """ 
        Lenght indep mean EHH statistics.. (Note: no data filtering!)
        """
        ehhcount = 0
        ehh = 0.0
        for i in xrange(0, len(self.snps) - 1):
            if (self.positions[len(self.snps) - 1] - self.positions[i]) >= windowSize:
                j = i + 1
                l = []
                while j < len(self.snps) and (self.positions[j] - self.positions[i]) < innerWindowSize :
                    j = j + 1
                while  j < len(self.snps) and (self.positions[j] - self.positions[i]) <= windowSize:
                    ehh = ehh + self.EHH(i, j)
                    ehhcount = ehhcount + 1
                    j = j + 1
            else:
                break
        return [ehh, ehhcount]








class SNPsDataSet:
    """
    A class that encompasses multiple _SnpsData_ chromosomes objects (chromosomes), and can deal with them as a whole.

    This object should eventually replace the snpsdata lists.. 
    """
    snpsDataList = None
    chromosomes = None
    accessions = None

    def __init__(self, snpsds, chromosomes, id=None, call_method=None, data_format=None):
        self.snpsDataList = snpsds
        self.chromosomes = chromosomes
        self.accessions = self.snpsDataList[0].accessions
        self.array_ids = self.snpsDataList[0].arrayIds
        self.id = id
        self.missing_val = snpsds[0].missingVal
        self.call_method = call_method
        self.data_format = data_format  # binary, diploid_ints, floats, int
        if not id and snpsds[0].id:
                self.id = id
        for i in range(1, len(self.chromosomes)):
            if self.accessions != self.snpsDataList[i].accessions:
                raise Exception("Accessions (or order) are different between SNPs datas")
        self.is_binary = list(snpsds[0].snps[0]).count(0) or list(snpsds[0].snps[0]).count(1)





#    def add_to_db(self, short_name, method_description='', data_description='', comment='', **kwargs):
#        """
#        Other possible keyword args are parent_id, accession_set_id ,imputed ,unique_ecotype 
#        """
#        conn = dbutils.connect_to_papaya()
#        cursor = conn.cursor()
#
#
#        #Checking whether data is already inserted...
#        sql_statement = "SELECT id FROM stock_250k.call_method WHERE short_name='%s';" % (short_name)
#        print sql_statement
#        cursor.execute(sql_statement)
#        row = cursor.fetchone()
#        print row
#        if row:
#            print "Data is already inserted in DB.  File will however be updated."
#            call_method_id = int(row[0])
#            filename = '/Network/Data/250k/db/dataset/call_method_%d.tsv' % call_method_id
#        else:
#
#            #Inserting data
#            sql_statement = "INSERT INTO stock_250k.call_method (short_name,method_description,data_description,comment"
#            for k in kwargs:
#                sql_statement += ',%s' % k
#            sql_statement += ") VALUES ('%s','%s','%s','%s'" % (short_name, method_description, data_description, comment)
#            for k in kwargs:
#                sql_statement += ',%s' % str(kwargs[k])
#            sql_statement += ');'
#            print sql_statement
#            cursor.execute(sql_statement)
#
#            #Getting method_id
#            sql_statement = "SELECT id FROM stock_250k.call_method WHERE short_name='%s';" % (short_name)
#            print sql_statement
#            cursor.execute(sql_statement)
#            row = cursor.fetchone()
#            call_method_id = int(row[0])
#            filename = '/Network/Data/250k/db/dataset/call_method_%d.tsv' % call_method_id
#
#            #Updating the filename
#            sql_statement = "UPDATE stock_250k.call_method SET filename='%s' WHERE id=%d" % (filename, call_method_id)
#            print sql_statement
#            cursor.execute(sql_statement)
#
#            print "Committing transaction (making changes permanent)."
#            conn.commit()
#
#        print "Call method id is %d" % call_method_id
#        #Generating Yu's file...
#        self.write_to_file_yu_format(filename)
#
#
#
#        #Generate DB call files...
#        self._generate_db_call_files_(call_method=call_method_id, cursor=cursor, conn=conn)
#
#        #Add SNPs to DB?
#
#
#        cursor.close()
#        conn.close()
#
#        return call_method_id
#
#
#
#
#    def _generate_db_call_files_(self, call_method=None, array_ids=None, file_dir='/Network/Data/250k/db/calls/', cursor=None, conn=None):
#        import os
#        import warnings
#
#        if not array_ids:
#            if not self.array_ids:
#                raise Exception('Array IDs are missing.')
#            else:
#                array_ids = self.array_ids
#        if not call_method:
#            raise Exception("Call method is missing!!")
#        chr_pos_snp_list = self.getChrPosSNPList()
#
#        #Create the call method directory
#        file_dir = file_dir + 'method_' + str(call_method) + '/'
#        if not os.path.lexists(file_dir):
#            os.mkdir(file_dir)
#        else:
#            warnings.warn('Directory already exists: %s' % file_dir)
#
#        #Connect to DB, if needed
#        if not cursor:
#            import dbutils
#            conn = dbutils.connect_to_papaya()
#            cursor = conn.cursor()
#
#        for i, aid in enumerate(array_ids):
#            print "Inserting genotype files into DB."
#            #Checking if it is already in the DB.
#            sql_statement = "SELECT id FROM stock_250k.call_info WHERE array_id=%s AND method_id=%d;"\
#                    % (aid, call_method)
#            print sql_statement
#            cursor.execute(sql_statement)
#            row = cursor.fetchone()
#            if row:
#                'Already in DB. File will be updated.'
#                call_info_id = int(row[0])
#            else:
#                #Insert info
#                sql_statement = "INSERT INTO stock_250k.call_info (array_id, method_id) VALUES (%s,%d);"\
#                        % (aid, call_method)
#                print sql_statement
#                cursor.execute(sql_statement)
#
#                sql_statement = "SELECT id FROM stock_250k.call_info WHERE array_id=%s AND method_id=%d;"\
#                        % (aid, call_method)
#                print sql_statement
#                cursor.execute(sql_statement)
#                row = cursor.fetchone()
#                if row:
#                    call_info_id = int(row[0])
#                print "Committing transaction (making changes permanent)."
#                conn.commit()
#
#            try:
#                #Write to a designated place.
#                file_name = file_dir + str(call_info_id) + "_call.tsv"
#                sql_statement = "UPDATE stock_250k.call_info \
#                         SET filename='%s'\
#                         WHERE id=%d" % (file_name, call_info_id)
#                print sql_statement
#                cursor.execute(sql_statement)
#
#                #Generate file in the right place
#                f = open(file_name, 'w')
#                f.write('SNP_ID\t%s\n' % aid)
#                for (c, p, s) in chr_pos_snp_list:
#                    f.write('%d_%d\t%s\n' % (c, p, s[i]))
#                f.close()
#            except Exception, err_str:
#                print "Couldn't generate call info file, error message:%s" % err_str
#        print "Closing connection."
#        #Close connection
#        if not cursor:
#            cursor.close()
#            conn.close()
#        print "Remember to copy files to papaya, i.e. everything in directory: %s" % file_dir



    def writeToFile(self, filename, delimiter=",", missingVal="NA", accDecoder=None,
            withArrayIds=False, decoder=None, callProbFile=None, binary_format=False):
        """
        Writes data to a file. 
        
        Note that there is no decoder dictionary option here..
        """

        print "Writing data to file:", filename
        numSnps = 0
        for i in range(0, len(self.chromosomes)):
            numSnps += len(self.snpsDataList[i].positions)


        # outStr = "NumSnps: "+str(numSnps)+", NumAcc: "+str(len(accessions))+"\n"
        if withArrayIds:
            outStr = ", ".join(["-", "-"] + self.snpsDataList[0].arrayIds) + "\n"
        else:
            outStr = ""
        fieldStrings = ["Chromosome", "Positions"]
        if accDecoder:
            for acc in self.snpsDataList[i].accessions:
                fieldStrings.append(str(accDecoder[acc]))
        else:
            for acc in self.snpsDataList[i].accessions:
                fieldStrings.append(str(acc))
        outStr += delimiter.join(fieldStrings) + "\n"
        if binary_format:
            f = open(filename, "wb")
        else:
            f = open(filename, "w")
        f.write(outStr)
        if decoder:
            for chromosome, snpsd in izip(self.chromosomes, self.snpsDataList):
                sys.stdout.write(".")
                sys.stdout.flush()
                for pos, snp in izip(snpsd.positions, snpsd.snps):
                    outStr = str(chromosome) + delimiter + str(pos)
                    for nt in snp:
                        outStr += delimiter + str(decoder[snp])
                    outStr += "\n"
                    f.write(outStr)
#                for j in range(0,len(self.snpsDataList[i].positions)):
#                    outStr =""
#                    outStr += str(self.chromosomes[i])+delimiter+str(self.snpsDataList[i].positions[j])
#                    for k in range(0, len(self.snpsDataList[0].accessions)):
#                        outStr += delimiter+str(decoder[self.snpsDataList[i].snps[j][k]])
#                    outStr +="\n"
#                    f.write(outStr)
        else:
            for chromosome, snpsd in izip(self.chromosomes, self.snpsDataList):
                sys.stdout.write(".")
                sys.stdout.flush()
                for pos, snp in izip(snpsd.positions, snpsd.snps):
                    outStr = '%d%s%d%s%s\n' % (chromosome, delimiter, pos, delimiter,
                                delimiter.join(map(str, snp.tolist())))
                    f.write(outStr)
#            for i in range(0,len(self.chromosomes)):
#                sys.stdout.write(".")
#                sys.stdout.flush()
#                for j in range(0,len(self.snpsDataList[i].positions)):
#                    outStr =""
#                    outStr += str(self.chromosomes[i])+delimiter+str(self.snpsDataList[i].positions[j])
#                    snp = self.snpsDataList[i].snps[j]
#                    if len(snp) != len(self.snpsDataList[0].accessions):
#                        print "len(snp):",len(snp),", vs. len(self.snpsDataList[0].accessions):",len(self.snpsDataList[0].accessions)
#                        raise Exception("The length didn't match")
#                    for nt in snp:
#                        outStr += delimiter+str(nt)
#                    outStr +="\n"
#                    f.write(outStr)
        f.close()
        print ""

        if callProbFile:
            if withArrayIds:
                outStr = "-, -, " + ", ".join(self.snpsDataList[0].arrayIds) + "\n"
            else:
                outStr = ""
            f = open(callProbFile, "w")
            outStr += delimiter.join(fieldStrings) + "\n"
            f.write(outStr)
            f.flush()
            for i in range(0, len(self.chromosomes)):
                outStr = ""
                snpsd = self.snpsDataList[i]
                self.snpsDataList[i] = []
                for j in range(0, len(snpsd.positions)):
                    outStr += str(self.chromosomes[i]) + delimiter + str(snpsd.positions[j])
                    for k in range(0, len(snpsd.accessions)):
                        outStr += delimiter + str(snpsd.callProbabilities[j][k])
                    outStr += "\n"
                del snpsd
                f.write(outStr)
                f.flush()
            f.close()


#
#    def write_to_file_yu_format(self, filename):
#        """
#        Writes data to a file in Yu's format. (Requires array IDs)
#        
#        Only works with raw sequence data.. (IUPAC nucleotides)
#        """
#        import yu_snp_key as yk
#        print "transposing chr_pos_snp_list"
#        chr_pos_snp_list = map(list, zip(*self.getChrPosSNPList()))
#        chrs = chr_pos_snp_list[0]
#        positions = chr_pos_snp_list[1]
#        snps = chr_pos_snp_list[2]
#        assert len(snps) == len(chrs) == len(positions), "SNPs, chromosomes, and positions not with same lenght"
#
#        print "transposing SNPs"
#        snps = map(list, zip(*snps))
#
#        array_ids = self.array_ids
#        ecotypes = self.accessions
#        assert len(snps) == len(array_ids) == len(ecotypes), "SNP, array IDs, and ecotype IDs not with same lenght"
#
#        print "len(array_ids) len(ecotypes):", len(array_ids), len(ecotypes)
#
#        print "Writing data to file:", filename
#        f = open(filename, "w")
#        str_list = ['ecotype_id', 'array_id']
#        str_list.extend([str(c) + "_" + str(p) for (c, p) in zip(chrs, positions)])
#        f.write('\t'.join(str_list) + '\n')
#        for ei, ai, nts in zip(ecotypes, array_ids, snps):
#            str_list = [str(ei), str(ai)]
#            str_list.extend([str(yk.nt_2_number[nt]) for nt in nts])
#            f.write('\t'.join(str_list) + '\n')
#        f.close()


    def coordinate_w_phenotype_data(self, phend, pid, coord_phen=True, verbose=False):

        """
        Deletes accessions which are not common, and sorts the accessions, removes monomorphic SNPs, etc.
        """
        print "Coordinating SNP and Phenotype data."
        ets = phend.phen_dict[pid]['ecotypes']
        # Checking which accessions to keep and which to remove.
#        common_ets = list(set(self.accessions + ets))
#        common_ets.sort()

        sd_indices_to_keep = set()  # []
        pd_indices_to_keep = []

        for i, acc in enumerate(self.accessions):
            for j, et in enumerate(ets):
                if et == acc:
                    sd_indices_to_keep.add(i)
                    pd_indices_to_keep.append(j)

        if verbose:
            initial_indices = set(range(len(ets)))
            missing_indices = list(initial_indices.difference(set(pd_indices_to_keep)))
            missing_indices.sort()
            if len(missing_indices) > 0:
                print "You phenotyped accessions that were not identified in the SNP-file."
                print [ets[i] for i in missing_indices ]

#
#        for i, acc in enumerate(self.accessions):
#            if common_ets[bisect.bisect(common_ets, acc) - 1] == acc:
#                    sd_indices_to_keep.add(i)
#        for j, et in enumerate(ets):
#            if common_ets[bisect.bisect(common_ets, et) - 1] == et:
#                    pd_indices_to_keep.append(j)
        sd_indices_to_keep = list(sd_indices_to_keep)
        sd_indices_to_keep.sort()



        # Filter accessions that do not have phenotype values (from the genotype data).
        print "Filtering genotype data"
        # if len(sd_indices_to_keep) != len(self.accessions):
        self.filter_accessions_indices(sd_indices_to_keep)
        if coord_phen:
            num_values = len(phend.phen_dict[pid]['ecotypes'])
            print "Filtering phenotype data."
            phend.filter_ecotypes(pd_indices_to_keep, pids=[pid])  # Removing accessions that don't have genotypes or phenotype values
            ets = phend.phen_dict[pid]['ecotypes']
            print "Out of %d, leaving %d values." % (num_values, len(ets))
        # Ordering accessions according to the order of accessions in the genotype file

#        if ets != self.accessions:
#            l = zip(ets, range(len(ets)))
#            l.sort()
#            l = map(list, zip(*l))
#            ets_map = l[1]
#            phend.order_ecotypes(ets_map, pids=[pid])


        if self.data_format == 'binary':
            print 'Filtering non-binary SNPs'
            total_num = 0
            removed_num = 0
            for snpsd in self.snpsDataList:
                total_num += len(snpsd.snps)
                removed_num += snpsd.onlyBinarySnps()
            print 'Removed %d non-binary SNPs out of %d SNPs' % (removed_num, total_num)
        elif self.data_format in ['int', 'diploid_int']:
            print 'Filtering monomorhpic SNPs'
            total_num = 0
            removed_num = 0
            for snpsd in self.snpsDataList:
                total_num += len(snpsd.snps)
                removed_num += snpsd.remove_monomorphic_snps()
            print 'Removed %d monomorphic SNPs out of %d SNPs' % (removed_num, total_num)
        return {'pd_indices_to_keep':pd_indices_to_keep, 'n_filtered_snps':removed_num}



    def get_ibs_kinship_matrix(self, debug_filter=1, snp_dtype='int8', dtype='single'):
        """
        Calculate the IBS kinship matrix. 
        (un-scaled)
        
        Currently it works only for binary kinship matrices.
        """
        import kinship
        print 'Starting kinship calculation'
        snps = self.getSnps(debug_filter)
        return kinship.calc_ibs_kinship(snps, snps_data_format=self.data_format, snp_dtype=snp_dtype,
                                        dtype=dtype)





    def get_local_n_global_kinships(self, focal_chrom_pos=None, window_size=25000, chrom=None, start_pos=None,
                    stop_pos=None, kinship_method='ibd', global_kinship=None, verbose=False):
        """
        Returns local and global kinship matrices.
        """
        if focal_chrom_pos != None:
            chrom, pos = focal_chrom_pos
            start_pos = pos - window_size
            stop_pos = pos + window_size

        local_snps, global_snps = self.get_region_split_snps(chrom, start_pos, stop_pos)
        if verbose:
            print 'Found %d local SNPs' % len(local_snps)
            print 'and %d global SNPs' % len(global_snps)
        if kinship_method == 'ibd':
            local_k = kinship.calc_ibd_kinship(local_snps) if len(local_snps) else None
            if global_kinship == None:
                global_k = kinship.calc_ibd_kinship(global_snps) if len(global_snps) else None
        elif kinship_method == 'ibs':
            local_k = kinship.calc_ibs_kinship(local_snps) if len(local_snps) else None
            if global_kinship == None:
                global_k = kinship.calc_ibs_kinship(global_snps) if len(global_snps) else None
        else:
            raise NotImplementedError
        if global_kinship != None:
            if len(local_snps):
                global_k = (global_kinship * self.num_snps() - local_k * len(local_snps)) / len(global_snps)
            else:
                local_k = None
                global_k = global_kinship
        return {'local_k':local_k, 'global_k':global_k, 'num_local_snps':len(local_snps),
            'num_global_snps':len(global_snps)}



    def get_chrom_vs_rest_kinships(self, chrom=None, kinship_method='ibd', global_kinship=None, verbose=False):
        """
        Returns local and global kinship matrices.
        """
        local_snps, global_snps = self.get_chrom_split_snps(chrom)
        if verbose:
            print 'Found %d SNPs on chromosome %d' % (len(local_snps), chrom)
            print 'and %d SNPs on other chromosomes' % len(global_snps)
        if kinship_method == 'ibd':
            local_k = self._calc_ibd_kinship_(local_snps, num_dots=0) if len(local_snps) else None
            if global_kinship == None:
                global_k = self._calc_ibd_kinship_(global_snps, num_dots=0) if len(global_snps) else None
        elif kinship_method == 'ibs':
            local_k = self._calc_ibs_kinship_(local_snps, num_dots=0) if len(local_snps) else None
            if global_kinship == None:
                global_k = self._calc_ibs_kinship_(global_snps, num_dots=0) if len(global_snps) else None
        else:
            raise NotImplementedError
        if global_kinship != None:
            if len(local_snps):
                global_k = (global_kinship * self.num_snps() - local_k * len(local_snps)) / len(global_snps)
            else:
                local_k = None
                global_k = global_kinship
        return {'local_k':local_k, 'global_k':global_k, 'num_local_snps':len(local_snps),
            'num_global_snps':len(global_snps)}






    def get_region_split_snps(self, chrom, start_pos, end_pos):
        """
        Returns two SNP sets, one with the SNPs within the given region, 
        and the other with the remaining SNPs.
        """
        global_snps = []
        local_snps = []
        chr_pos_l = self.get_chr_pos_list(cache_list=True)
        start_i = bisect.bisect(chr_pos_l, (chrom, start_pos))
        stop_i = bisect.bisect(chr_pos_l, (chrom, end_pos))
        snps = self.get_snps(cache=True)
        local_snps = snps[start_i:stop_i]
        global_snps = snps[:start_i] + snps[stop_i:]
        return local_snps, global_snps

    def get_chrom_split_snps(self, chrom):
        c_ends = self.get_chromosome_ends()
        ci = self.chromosomes.index(chrom)
        return self.get_region_split_snps(chrom, 0, c_ends[ci] + 1)


    def get_region_split_kinships(self, chrom_pos_list, kinship_method='ibd', global_kinship=None, verbose=False):
        """
        Returns local and global kinship matrices.
        """
        region_snps = self.get_regions_split_snps(chrom_pos_list)
        num_snps_found = len(region_snps)
        if verbose:
            print 'Found %d SNPs in the regions' % num_snps_found
        if kinship_method == 'ibd':
            regional_k = self._calc_ibd_kinship_(region_snps, num_dots=0) if num_snps_found else None
        elif kinship_method == 'ibs':
            regional_k = self._calc_ibs_kinship_(region_snps, num_dots=0) if num_snps_found else None
        else:
            raise NotImplementedError
        if global_kinship != None and regional_k != None:
            global_k = (global_kinship * self.num_snps() - regional_k * num_snps_found) \
                    / (self.num_snps() - num_snps_found)
        else:
            global_k = None
        return {'regional_k':regional_k, 'global_k':global_k, 'num_snps_found':num_snps_found}



    def get_regions_split_snps(self, chrom_pos_list):
        """
        Returns the SNP sets, wherethe SNPs are in multiple regions, 
        and the other is the remaining set. 
        """
        chr_pos_l = self.get_chr_pos_list()
        snps = self.get_snps()
        region_snps = []
        for chrom, start_pos, end_pos in chrom_pos_list:
            start_i = bisect.bisect(chr_pos_l, (chrom, start_pos))
            stop_i = bisect.bisect(chr_pos_l, (chrom, end_pos))
            region_snps.extend(snps[start_i:stop_i])
        return region_snps

#    def get_ibd_kinship_matrix_old(self, debug_filter=1, num_dots=10000, with_correction=True,
#                snp_dtype='int8', dtype='single'):
#        """
#        Calculate the IBD kinship matrix, as described in (Yang et al., Nat. Genetics, 2010) 
#        (un-scaled)
#        """
#        if self.data_format != 'binary':
#            raise NotImplementedError
#        print 'Starting IBD kinship calculation, it prints %d dots.' % num_dots
#        snps = self.getSnps(debug_filter)
#        num_snps = len(snps)
#        num_lines = len(self.accessions)
#        print 'Allocating K matrix'
#        k_mat = sp.zeros((num_lines, num_lines), dtype=dtype)
#
#        print 'Calculating IBD kinship... one SNP at a time'
#        #Do one SNP at a time to save memory...
#        for snp_i, snp in enumerate(snps):
#            p = sp.mean(snp)
#            norm_snp = sp.mat((snp - p) / sp.std(snp))
#            M = norm_snp.T * norm_snp
#            if with_correction:
#                c = p * (1 - p * (4 - 6 * p + 3 * p * p))
#                cor_norm_snp = (snp - p) / sp.sqrt(c)
#                for i in range(0, num_lines):
#                    M[i, i] = cor_norm_snp[i] ** 2
#            k_mat += M#norm_snp.T * norm_snp
#            if num_snps >= num_dots and (snp_i + 1) % (num_snps / num_dots) == 0: #Print dots
#                sys.stdout.write('.')
#                sys.stdout.flush()
#        k_mat = k_mat / len(snps)
#        return k_mat


    def get_ibd_kinship_matrix(self, debug_filter=1, dtype='single'):
        import kinship
        print 'Starting IBD calculation'
        snps = self.getSnps(debug_filter)
        cov_mat = kinship.calc_ibd_kinship(snps, len(self.accessions), dtype=dtype)
        print 'Finished calculating IBD kinship matrix'
        return cov_mat


    def get_snp_cov_matrix(self, debug_filter=1, num_dots=10, dtype='single'):
        print 'Initializing'
        num_lines = len(self.accessions)
        chunk_size = num_lines
        num_snps = self.num_snps()
        num_splits = num_snps / chunk_size
        cov_mat = sp.zeros((num_lines, num_lines))
        snps = self.getSnps(debug_filter)

        print 'Estimating the accession means'
        accession_sums = sp.zeros(num_lines)
        for chunk_i, i in enumerate(range(0, self.num_snps(), chunk_size)):
            snps_array = sp.array(snps[i:i + chunk_size])
            snps_array = snps_array.T
            norm_snps_array = (snps_array - sp.mean(snps_array, 0)) / sp.std(snps_array, 0)
            accession_sums += sp.sum(norm_snps_array, 1)
        accession_means = accession_sums / num_snps
        print accession_means

        print 'Starting covariance calculation'
        for chunk_i, i in enumerate(range(0, self.num_snps(), chunk_size)):
            snps_array = sp.array(snps[i:i + chunk_size])
            snps_array = snps_array.T
            norm_snps_array = (snps_array - sp.mean(snps_array, 0)) / sp.std(snps_array, 0)
            x = sp.mat(norm_snps_array.T - accession_means)  # The only difference from the IBD matrix is that we subtract the accession means.
            cov_mat += x.T * x
            if num_splits >= num_dots and (chunk_i + 1) % int(num_splits / num_dots) == 0:  # Print dots
                sys.stdout.write('.')
                sys.stdout.flush()
        cov_mat = cov_mat / self.num_snps()
        print 'Finished calculating covariance matrix'
        return cov_mat



    def get_macs(self):
        r = self.get_mafs()
        return r['mafs']



    def get_normalized_snps(self, debug_filter=1, dtype='single'):
        print 'Normalizing SNPs'
        snps = self.getSnps(debug_filter)
        snps_array = sp.array(snps)
        snps_array = snps_array.T
        # Normalizing (subtracting by)
        norm_snps_array = (snps_array - sp.mean(snps_array, 0)) / sp.std(snps_array, 0)
        print 'Finished normalizing them'
        return norm_snps_array



    def convert_data_format(self, target_data_format='binary'):
        """
        Converts the underlying raw data format to a binary one, i.e. A,C,G,T,NA,etc. are converted to 0,1,-1
        """
        if self.data_format == target_data_format:
            import warnings
            warnings.warn("Data appears to be already in %s format!" % target_data_format)
        else:
            if self.data_format == 'nucleotides' and target_data_format == 'binary':
                snpsd_list = []
                for snpsd in self.snpsDataList:
                    snpsd_list.append(snpsd.getSnpsData())
                self.snpsDataList = snpsd_list
                self.data_format = 'binary'
                self.missing_val = self.snpsDataList[0].missingVal
            else:
                raise NotImplementedError


    def haplotize(self, snp_window=None, base_window=None):
        """
        Converts the data-format to a haplotype numbering format
        """
        print 'Haplotizing!'
        assert (snp_window or base_window), 'snp_window or base_window arguments are missing.'
        for i, sd in enumerate(self.snpsDataList):
            sd.haplotize(snp_window=snp_window, base_window=base_window)
            print i



    def impute_missing(self, reference_data):
        if not len(self.snpsDataList) == len(reference_data.snpsDataList):
            raise Exception("Reference data for imputation doesn't have equal chromosome number.")
        for i, snpsd in enumerate(self.snpsDataList):
            snpsd.inpute_data_region(reference_data.snpsDataList[i])




    def get_snps(self, random_fraction=None, cache=False):
        if cache:
            try:
                return self.snps
            except Exception:
                pass
        snplist = []
        if random_fraction:
            import random
            for snpsd in self.snpsDataList:
                for snp in snpsd.snps:
                    if random.random() < random_fraction:
                        snplist.append(snp)
        else:
            for snpsd in self.snpsDataList:
                snplist.extend(snpsd.snps)
                if cache:
                        self.snps = snplist
        return snplist


    def getSnps(self, random_fraction=None):
        return self.get_snps(random_fraction=None)


    def num_snps(self):
        num_snps = 0
        for snpsd in self.snpsDataList:
            num_snps += len(snpsd.snps)
        return num_snps


    def plot_tree(self, tree_file, verbose=True, kinship_method='ibs'):

        if verbose:
            print "Calculating kinship matrix"
        if kinship_method == 'ibs':
            K = self.get_ibs_kinship_matrix()
        if kinship_method == 'ibd':
            K = self.get_ibd_kinship_matrix()
        plot_tree(K, tree_file, self.accessions, verbose=verbose)


    def plot_snp_map(self, chromosome, position, pdf_file=None, png_file=None, map_type='global',
            color_by=None, cmap=None, title='', eid=None):
        """
        Plot accessions on a map.
        
        'color_by' is by default set to be the phenotype values.
        """
        import phenotypeData as pd
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        # matplotlib.rcParams['backend'] = 'GTKAgg'
        if eid == None:
            eid = pd.get_ecotype_id_info_dict()
        lats = []
        lons = []
        acc_names = []
        for e in self.accessions:
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
            # m.bluemarble()
        elif map_type == 'sweden':

            plt.figure(figsize=(2.4, 4))
            plt.axes([0.02, 0.02, 0.96, 0.96])
            m = Basemap(projection='merc', llcrnrlat=55, urcrnrlat=67,
                    llcrnrlon=10, urcrnrlon=25, lat_ts=10, resolution='i')
            m.drawparallels(np.arange(45, 75, 5))
            m.drawmeridians(np.arange(5, 30, 5))
            # m.bluemarble()
        else:
            raise Exception("map_type is invalid")

        # m.drawmapboundary(fill_color='aqua')
        m.drawcoastlines()
        m.fillcontinents()
        m.drawcountries()
        # m.fillcontinents(color='green', lake_color='blue')

        xs = []
        ys = []
        for lon, lat in zip(lons, lats):
            x, y = m(*np.meshgrid([lon], [lat]))
            xs.append(float(x))
            ys.append(float(y))

        if not color_by:
            color_vals = self.get_snp_at(chromosome, position)
        else:
            color_vals = color_by
        assert len(color_vals) == len(self.accessions), "accessions and color_by_vals values don't match ! "
        if not cmap:
            num_colors = len(set(color_vals))
            if num_colors <= 10:
                cmap = cm.get_cmap('jet', num_colors)
            else:
                cmap = cm.get_cmap('jet')
        lws = [0] * len(xs)
        plt.scatter(xs, ys, s=10, linewidths=lws, c=color_vals, cmap=cmap, alpha=0.7, zorder=2)
        # plt.plot(xs, ys, 'o', color='r', alpha=0.5, zorder=2,)
        if title:
            plt.title(title)
        if pdf_file:
            plt.savefig(pdf_file, format="pdf", dpi=400)
        if png_file:
            plt.savefig(png_file, format="png", dpi=400)
        if not pdf_file and not png_file:
            plt.show()

        return self.accessions, lats, lons


#    def get_snps(self, random_fraction=None, region=None):
#        snplist = []
#        if random_fraction:
#            import random
#            for snpsd in self.snpsDataList:
#                for snp in snpsd.snps:
#                    if random.random() < random_fraction:
#                        snplist.append(snp)
#        elif region:
#            chrom, start_pos, end_pos = region
#            snpsd = self.snpsDataList[chrom - 1]
#            ps_gen = izip(snpsd.positions, snpsd.snps)
#            p, s = ps_ge.next()
#            while p < start_pos:
#                p, s = ps_ge.next()
#            #Now found
#            while p < end_pos:
#                snplist.append(s)
#                pos_list = []
#                p, s = ps_ge.next()
#
#        else:
#            for snpsd in self.snpsDataList:
#                snplist += snpsd.snps
#        return snplist


    def get_snp_at(self, chromosome, position):
        """
        Returns the SNP at the given position, if it exits.
        """
        c_i = self.chromosomes.index(chromosome)
        sd = self.snpsDataList[c_i]
        i = 0
        while sd.positions[i] < position:
            i += 1
        if sd.positions[i] == position:
            print 'Found the SNP.'
            return sd.snps[i]
        else:
            print "Didn't find the SNP on chromosome %d, at position %d" % (chromosome, position)
            return None


    def get_positions(self):
        poslist = []
        for snpsd in self.snpsDataList:
            for pos in snpsd.positions:
                poslist.append(pos)
        return poslist

    def getPositions(self):
        return self.get_positions()

    def get_top_correlated_snp(self, snp, r2_threshold=0.5):
        sample_snp_chr_pos_marf = []
        sample_r2s = []
        print_r2s = []
        print_snp_chr_pos_marf = []
        for snpsd, chromosome in zip(self.snpsDataList, self.chromosomes):
            marfs = snpsd.get_mafs()['marfs']
            r2s = snpsd.calc_r2_with_snp(snp)
            for r2, sample_snp, sample_snp_pos, marf in zip(r2s, snpsd.snps, snpsd.positions, marfs):
                if r2 > r2_threshold:
                    sample_r2s.append(r2)
                    sample_snp_chr_pos_marf.append((sample_snp, chromosome, sample_snp_pos, marf))
                if r2 > 0.1:
                    print_r2s.append(r2)
                    print_snp_chr_pos_marf.append((sample_snp, chromosome, sample_snp_pos, marf))

        l = zip(print_r2s, print_snp_chr_pos_marf)
        l.sort()
        print l[-10:]
        return sample_snp_chr_pos_marf


    def get_all_snp_w_info(self):
        """
        Returns the SNPs along with some info..
        """
        sample_snp_chr_pos_marf = []
        for snpsd, chromosome in zip(self.snpsDataList, self.chromosomes):
            marfs = snpsd.get_mafs()['marfs']
            sample_snp_chr_pos_marf += zip(snpsd.snps, [chromosome] * len(snpsd.snps), snpsd.positions, marfs)
        return sample_snp_chr_pos_marf


    def get_chr_pos_list(self, cache_list=False):
        if cache_list:
            try:
                chr_pos_list = self.chr_pos_list
                if len(chr_pos_list) > 0:
                    return chr_pos_list
                else:
                    raise Exception
            except Exception:
                pass
        chr_pos_list = []
        for chrom, snpsd in izip(self.chromosomes, self.snpsDataList):
            chr_pos_list += zip([chrom] * len(snpsd.positions), snpsd.positions)
        if cache_list:
            self.chr_pos_list = chr_pos_list
        return chr_pos_list

    def getChrPosList(self):
        return self.get_chr_pos_list()

    def get_chr_list(self):
        chr_list = []
        for c, snpsd in izip(self.chromosomes, self.snpsDataList):
            chr_list.extend([c] * len(snpsd.positions))
        return chr_list


    def get_chromosome_ends(self):
        chr_ends = []
        for snpsd in self.snpsDataList:
            chr_ends.append(snpsd.positions[-1])
        return chr_ends


    def get_genome_length(self):
        return sum(self.get_chromosome_ends())

    def get_mafs(self):
        """
        Returns the mafs and marfs as a dictionary.
        
        types: classes, diploid_ints
        """
        maf_list = []
        marf_list = []
        for snpsd in self.snpsDataList:
            r = snpsd.get_mafs(type=self.data_format)
            maf_list.extend(r["mafs"])
            marf_list.extend(r["marfs"])
        print "Finished calculating MAFs."
        return {"mafs":maf_list, "marfs":marf_list}


    def get_chr_pos_snp_list(self):
        chr_pos_snp_list = []
        for i in range(0, len(self.snpsDataList)):
            snpsd = self.snpsDataList[i]
            chr = i + 1
            for j in range(0, len(snpsd.positions)):
                pos = snpsd.positions[j]
                snp = snpsd.snps[j]
                chr_pos_snp_list.append((chr, pos, snp))
        return chr_pos_snp_list

    def getChrPosSNPList(self):
        return self.get_chr_pos_snp_list()



    def get_region_pos_snp_dict(self, chromosome, start_pos=None, end_pos=None):
        """
        Returns a dict containing a list of positions and snps.
        """
        positions = []
        snps = []
        for snpsd in self.snpsDataList:
            if snpsd.chromosome == chromosome:
                break
        assert snpsd.chromosome == chromosome, "SNPs data with appropriate chromosomes wasn't found"
        i = 0
        while i < len(snpsd.positions) and snpsd.positions[i] < start_pos:
            i += 1
        if end_pos:
            while i < len(snpsd.positions) and snpsd.positions[i] < end_pos:
                positions.append(snpsd.positions[i])
                snps.append(snpsd.snps[i])
                i += 1
        else:
            positions = snpsd.positions[i:]
            snps = snpsd.snps[i:]
        return {'positions':positions, 'snps':snps}



    def get_cand_genes_snp_priors(self, cand_genes, radius=10000, num_exp_causal=1.0, cg_prior_fold_incr=50.0,
                    method_type='sum_all_priors'):
        """
        Returns SNP priors
        """
        chr_pos_list = self.getChrPosList()
        num_snps = float(len(chr_pos_list))
        i = 0
        gene_chr_pos_list = sorted([(int(gene.chromosome), gene.startPos, gene.endPos) for gene in cand_genes])
        num_cgs = len(gene_chr_pos_list)
        cg_snp_indices = []
        for gene_chr, gene_start_pos, gene_end_pos in gene_chr_pos_list:
            start_i = bisect.bisect(chr_pos_list, (gene_chr, gene_start_pos - radius - 1))
            stop_i = bisect.bisect(chr_pos_list, (gene_chr, gene_end_pos + radius + 1))
            cg_snp_indices.extend(range(start_i, stop_i))
        cg_snp_indices = sorted(list(set(cg_snp_indices)))
        if method_type == 'sum_all_priors':
            pi_0 = num_exp_causal / num_snps  # Basis prior
            pi_1 = pi_0 * cg_prior_fold_incr  # Cand. gene prior
        elif method_type == 'sum_base_priors':
            pi_0 = num_exp_causal / (num_snps - num_cgs + cg_prior_fold_incr * num_cgs)  # Basis prior
            pi_1 = pi_0 * cg_prior_fold_incr  # Cand. gene prior
        else:
            raise NotImplementedError
        snp_priors = sp.repeat(pi_0, num_snps)
        for i in cg_snp_indices:
            snp_priors[i] = pi_1
        print max(snp_priors), min(snp_priors)
        return snp_priors.tolist()




    def get_snp_priors(self, cpp_list=None, cand_genes=None, radius=25000, num_exp_causal=10.0, cg_prior_fold_incr=10):
        """
        Takes a list of SNPs/markers with some priors, and extrapolates that to the SNPs in the data.
        """
        cpp_list.sort()
        l = map(list, zip(*cpp_list))
        priors = l[2]
        cp_list = self.getChrPosList()
        snp_priors = []
        snp_i = 0
        p_i = 0


        if cpp_list != None:
            # Do chromosome by chromosome..
            for chrom in [1, 2, 3, 4, 5]:
                p_i = bisect.bisect(cpp_list, (chrom, 0, 0))
                snp_i = 0
                positions = self.snpsDataList[chrom - 1].positions
                num_snps = len(positions)
                pos = positions[snp_i]
                chrom_1, pos_1, prior_1 = cpp_list[p_i]
                chrom_2, pos_2, prior_2 = cpp_list[p_i + 1]
                while snp_i < num_snps and p_i < len(cpp_list) - 2 and chrom_2 == chrom:
                    chrom_1, pos_1, prior_1 = cpp_list[p_i]
                    chrom_2, pos_2, prior_2 = cpp_list[p_i + 1]
                    if chrom_2 == chrom:
                        while snp_i < num_snps - 1 and pos <= pos_1:
                            snp_priors.append(prior_1)
                            snp_i += 1
                            pos = positions[snp_i]

                        while snp_i < num_snps - 1 and pos_1 < pos <= pos_2:
                            d = pos_2 - pos_1
                            s = (pos - pos_1) / float(d)
                            snp_priors.append(prior_1 * (1 - s) + prior_2 * s)
                            snp_i += 1
                            pos = positions[snp_i]

                    p_i += 1
                    chrom_1, pos_1, prior_1 = cpp_list[p_i]
                    chrom_2, pos_2, prior_2 = cpp_list[p_i + 1]

                if chrom_2 != chrom:  # The chromosome is ending
                    while snp_i < num_snps:
                        snp_priors.append(prior_1)
                        snp_i += 1

                elif p_i >= len(cpp_list) - 2:  # It's finishing
                    while snp_i < num_snps:
                        snp_priors.append(prior_2)
                        snp_i += 1
            snp_priors = sp.array(snp_priors)
            snp_priors = 1000 * (snp_priors - snp_priors.min()) / (snp_priors.max() - snp_priors.min()) + 1
        else:
            snp_priors = sp.repeat(1.0, self.num_snps())


        if cand_genes != None:
            print 'Now for candidate genes'
            chr_pos_list = self.getChrPosList()
            num_snps = len(chr_pos_list)
            i = 0
            gene_chr_pos_list = sorted([(int(gene.chromosome), gene.startPos, gene.endPos) for gene in cand_genes])
            print gene_chr_pos_list
            num_cgs = len(gene_chr_pos_list)
            cg_snp_indices = []
            for gene_chr, gene_start_pos, gene_end_pos in gene_chr_pos_list:
                start_i = bisect.bisect(chr_pos_list, (gene_chr, gene_start_pos - radius - 1))
                stop_i = bisect.bisect(chr_pos_list, (gene_chr, gene_end_pos + radius + 1))
                cg_snp_indices.extend(range(start_i, stop_i))
            cg_snp_indices = sorted(list(set(cg_snp_indices)))
            snp_priors[cg_snp_indices] = snp_priors[cg_snp_indices] * cg_prior_fold_incr

        snp_priors = (num_exp_causal * snp_priors / sum(snp_priors)).tolist()
        print max(snp_priors), min(snp_priors)

        return snp_priors




#    def get_pc(self, pc_num=1, random_fraction=0.1):
#        """
#        Returns the pc_num'th principal components of the genotype 
#        """
#        import random
#        import rpy, util
#
#        if not self.is_binary:
#            print "Converting the snps data to binary format."
#            self.convert_2_binary()
#
#        snps = []
#        for sd in self.snpsDataList:
#            snps.extend(random.sample(sd.snps, int(random_fraction * len(sd.snps))))
#        genotypes = map(list, zip(*snps))
#
#
#        for genotype in genotypes:
#            sd = util.calcSD(genotype)
#            for j in range(len(genotype)):
#                genotype[j] = genotype[j] / sd
#
#        genotypes = sp.transpose(sp.array(genotypes))
#        #print genotypes
#        pc = rpy.r.princomp(genotypes)
#        pc_sorted = zip(list(pc["scores"][pc_num - 1]), self.accessions)
#        pc_sorted.sort()
#        print pc_sorted
#        pc = list(pc["scores"][pc_num - 1])
#        self.pc = pc
#        return pc


    def updateRegions(self, regionList):
        """
        Deprecated 11/11/08 - Bjarni
        """
        c_i = 0
        i = 0
        rl_i = 0  # region list index
        while c_i < len(self.chromosomes) and rl_i < len(regionList):
            region = regionList[rl_i]
            snpsd = self.snpsDataList[c_i]
            cp1 = (c_i + 1, snpsd.positions[i])
            cp_start = (region.chromosome, region.startPos)
            while cp1 < cp_start:
                if i < len(snpsd.positions):
                    i += 1
                else:
                    c_i += 1
                    snpsd = self.snpsDataList[c_i]
                    i = 0
                cp1 = (c_i + 1, snpsd.positions[i])
            cp_end = (region.chromosome, region.endPos)
            while cp1 <= cp_end:
                """Update current region!"""
                region.snps.append(snpsd.snps[i])
                region.snps_indices.append((c_i, i))

                i += 1
                if i < len(snpsd.positions):
                    cp1 = (c_i + 1, snpsd.positions[i])
                else:
                    c_i += 1
                    i = 0
                    break

            rl_i += 1

    def filter_na_snps(self, max_na_rate=0.0):
        for snpsd in self.snpsDataList:
            snpsd.filter_na_snps(max_na_rate=max_na_rate)

    def remove_snps_indices(self, indices_to_remove):
        remove_indices = [[] for i in range(len(self.snpsDataList))]
        indices_to_remove.sort()
        offset = 0
        i = 0
        for snpsd_i, snpsd in enumerate(self.snpsDataList):
            max_i = offset + len(snpsd.snps)
            i2r = indices_to_remove[i]
            while i < len(indices_to_remove) and i2r < max_i:
                    remove_indices[snpsd_i].append(i2r - offset)
                    i += 1
                    i2r = indices_to_remove[i]
            offset = max_i
        for i, snpsd in enumerate(self.snpsDataList):
            snpsd.remove_snps(remove_indices[i])


    def filter_na_accessions(self, max_na_rate=0.2, verbose=False):
        accessions_na_counts = [0.0 for acc in self.accessions]
        total_snps = 0
        for snpsd in self.snpsDataList:
            for i, c in enumerate(snpsd.countMissingSnps()[1]):
                accessions_na_counts[i] += c
            total_snps += float(len(snpsd.snps))
        accessions_na_counts = [accessions_na_counts[i] / total_snps for i in range(len(accessions_na_counts))]
        # print accessions_na_counts
        acc_to_keep = []
        for i, na_rate in enumerate(accessions_na_counts):
            if na_rate <= max_na_rate:
                acc_to_keep.append(self.accessions[i])
        # print len(acc_to_keep)
        self.filter_accessions(acc_to_keep)



    def filter_maf_snps(self, maf, maf_ub=1):
        for snpsd in self.snpsDataList:
            snpsd.snpsFilterMAF([maf, maf_ub], type=self.data_format)


    def filter_mac_snps(self, mac_threshold=15):
        for snpsd in self.snpsDataList:
            snpsd.filter_mac(mac_threshold, data_format=self.data_format)


    def filter_monomorphic_snps(self):
        n_filtered_snps = 0
        for snpsd in self.snpsDataList:
            n_filtered_snps += snpsd.filterMonoMorphicSnps()
        return n_filtered_snps


    def filter_snps(self, snps_to_keep):
        """
        Filter snps, leaving the snps in the given order.
        """
        assert len(snps_to_keep) != 0, "Can't remove all snps."
        
        for i, snps in enumerate(snps_to_keep):
            snps_indices_to_keep = []
            for snp in snps:
                snps_indices_to_keep.append(bisect.bisect(self.snpsDataList[i].positions, snp) - 1)
            self.snpsDataList[i].filter_snp_indices(snps_indices_to_keep)
        
   
        
        
    def filter_accessions(self, accessions_to_keep):
        """
        Filter accessions, leaving the remaining accession in the given order.
        """
        assert len(accessions_to_keep) != 0, "Can't remove all ecotypes."
        ecotypes = accessions_to_keep
        acc_indices_to_keep = []
        for et in ecotypes:
            try:
                i = self.accessions.index(et)
                acc_indices_to_keep.append(i)
            except:
                continue
        # pdb.set_trace()
        # acc_indices_to_keep.sort()
        self.filter_accessions_indices(acc_indices_to_keep)


    def filter_accessions_indices(self, acc_indices_to_keep):
        num_accessions = len(self.accessions)
        for i, snpsd in enumerate(self.snpsDataList):
            # print i, len(snpsd.accessions)
            snpsd.removeAccessionIndices(acc_indices_to_keep)
        self.accessions = self.snpsDataList[0].accessions
        self.array_ids = self.snpsDataList[0].arrayIds
        print "Removed %d accessions, leaving %d in total." % (num_accessions - len(acc_indices_to_keep), len(acc_indices_to_keep))


    def filter_for_countries(self, country_codes, complement=False):
        import phenotypeData as pd
        ei_dict = pd.get_ecotype_id_info_dict()
        acc_indices_to_keep = []
        if complement:
            for i, ei in enumerate(self.accessions):
                if ei_dict[int(ei)][4] not in country_codes:
                    acc_indices_to_keep.append(i)
        else:
            for i, ei in enumerate(self.accessions):
                if ei_dict[int(ei)][4] in country_codes:
                    acc_indices_to_keep.append(i)
        self.filter_accessions_indices(acc_indices_to_keep)

        for ei in self.accessions:
            print ei, ei_dict[int(ei)]


    def sample_snps(self, random_fraction, seed=None):
        """
        Samples a random fraction of the SNPs.
        """

        for snpsd in self.snpsDataList:
            snpsd.sample_snps(random_fraction, seed=seed)


    def get_region_snpsd(self, chr, start_pos=None, end_pos=None):
        """
        Modifies original object, beware!
        """
        c_i = self.chromosomes.index(chr)
        print chr, start_pos, end_pos
        snpsd = self.snpsDataList[c_i]
        if start_pos != None:
            new_snps = []
            new_positions = []
            i = 0
            while i < len(snpsd.positions) - 1 and snpsd.positions[i] < start_pos:
                i += 1
            while i < len(snpsd.positions) - 1 and snpsd.positions[i] < end_pos:
                new_snps.append(snpsd.snps[i])
                new_positions.append(snpsd.positions[i])
                i += 1
            snpsd.positions = new_positions
            snpsd.snps = new_snps
        return snpsd

    def merge_snps_data(self, sd, acc_merge_type='intersection', error_threshold=0.1, discard_error_threshold=0.1):
        """
        Merges data using _SnpsData_.merge_data
        """
        if self.chromosomes != sd.chromosomes:
            raise Exception
        else:
            self.new_snps_data_list = []
            for sd1, sd2, chromosome in zip(self.snpsDataList, sd.snpsDataList, self.chromosomes):
                print "Merging data on chromosome %s." % (str(chromosome))
                sd1.merge_data(sd2, acc_merge_type=acc_merge_type, error_threshold=error_threshold,
                        discard_error_threshold=discard_error_threshold)
                self.new_snps_data_list.append(sd1)
            self.accessions = self.new_snps_data_list[0].accessions
            self.snpsDataList = self.new_snps_data_list








def plot_tree(K, tree_file, ets, verbose=True, label_values=None):
        import scipy.cluster.hierarchy as hc
        import pylab
        import phenotypeData
        e_dict = phenotypeData.get_ecotype_id_info_dict()
        # print e_dict
        labels = []
        for et_i, et in enumerate(ets):
                try:
                        s1 = unicode(e_dict[int(et)][0], 'iso-8859-1')
                        if label_values != None:
                                s2 = unicode(' %s' % str(label_values[et_i]))
                        else:
                                s2 = unicode('(%0.1f,%0.1f)' % (e_dict[int(et)][2], e_dict[int(et)][3]))
                        s = s1 + s2
                except Exception, err_s:
                        print err_s
                        s = str(et)
                labels.append(s)
        if verbose:
                print "Plotting tree for SNPs:"
        Z = hc.average(K)
        pylab.figure(figsize=(24, 15))
        pylab.axes([0.03, 0.08, 0.96, 0.91])
        dend_dict = hc.dendrogram(Z, leaf_font_size=7, labels=labels)
        xmin, xmax = pylab.xlim()
        xrange = xmax - xmin
        ymin, ymax = pylab.ylim()
        yrange = ymax - ymin
        pylab.axis([xmin - 0.01 * xrange, xmax + 0.01 * xrange, ymin - 0.02 * yrange, ymax + 0.02 * yrange])
        pylab.savefig(tree_file, format='pdf')
        pylab.clf()
        if verbose:
                print "Done plotting tree, saved in file:", tree_file, "\n"




def construct_snps_data_set(snps, positions, chromosomes, indiv_ids, data_format='binary'):
    """
    Assumes the SNPs are sorted by chromosome and position.
    """
    last_i = 0
    snpsds = []
    chrom_list = list(set(chromosomes))
    chrom_list.sort()
    for chrom in chrom_list[:-1]:
        i = bisect.bisect(chromosomes, chrom)
        snpsds.append(SNPsData(snps=snps[last_i:i], positions=positions[last_i:i],
                    chromosome=chrom, accessions=indiv_ids))
        last_i = i
    assert snps.dtype == sp.dtype('int8'), "Type doesn't match the data format."
    snpsds.append(SNPsData(snps=snps[last_i:], positions=positions[last_i:],
                chromosome=chrom_list[-1], accessions=indiv_ids))
    return SNPsDataSet(snpsds, chrom_list, data_format=data_format)



if __name__ == "__main__":
    pass

