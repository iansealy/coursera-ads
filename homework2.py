#!/usr/bin/env python

"""Algorithms for DNA Sequencing - Programming Homework 2"""

from bm_preproc import BoyerMoore
from kmer_index import Index
import bisect


def main():
    chr1 = readGenome('chr1.GRCh38.excerpt.fasta')
    p = 'GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG'
    p_bm = BoyerMoore(p)
    print(naive(p, chr1)[2])
    print(naive(p, chr1)[1])
    print(boyer_moore(p, p_bm, chr1)[2])
    p = 'GGCGCGGTGGCTCACGCCTGTAAT'
    print(len(approximate_match(p, chr1, 2)[0]))
    print(approximate_match(p, chr1, 2)[1])
    print(approximate_match_subseq(p, chr1, 2, 3)[1])


def boyer_moore(p, p_bm, t):
    """ Do Boyer-Moore matching. p=pattern, t=text,
        p_bm=BoyerMoore object for p """
    i = 0
    occurrences = []
    comparisons = 0
    alignments = 0
    while i < len(t) - len(p) + 1:
        alignments += 1
        shift = 1
        mismatched = False
        for j in range(len(p)-1, -1, -1):
            comparisons += 1
            if p[j] != t[i+j]:
                skip_bc = p_bm.bad_character_rule(j, t[i+j])
                skip_gs = p_bm.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                break
        if not mismatched:
            occurrences.append(i)
            skip_gs = p_bm.match_skip()
            shift = max(shift, skip_gs)
        i += shift
    return occurrences, comparisons, alignments
    
    
def naive(p, t):
    occurrences = []
    comparisons = 0
    alignments = 0
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        alignments += 1
        match = True
        for j in range(len(p)):  # loop over characters
            comparisons += 1
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences, comparisons, alignments


def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome


def approximate_match(p, t, n):
    segment_length = int(round(len(p) / (n+1)))
    all_matches = set()
    p_idx = Index(t, segment_length)
    idx_hits = 0
    for i in range(n+1):
        start = i*segment_length
        end = min((i+1)*segment_length, len(p))
        matches = p_idx.query(p[start:end])
        
        # Extend matching segments to see if whole p matches
        for m in matches:
            idx_hits += 1
            if m < start or m-start+len(p) > len(t):
                continue
            
            mismatches = 0
            
            for j in range(0, start):
                if not p[j] == t[m-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break
            for j in range(end, len(p)):
                if not p[j] == t[m-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break
            
            if mismatches <= n:
                all_matches.add(m - start)
    return list(all_matches), idx_hits


def approximate_match_subseq(p, t, n, ival):
    segment_length = int(round(len(p) / (n+1)))
    all_matches = set()
    p_idx = SubseqIndex(t, segment_length, ival)
    idx_hits = 0
    for i in range(n+1):
        start = i
        matches = p_idx.query(p[start:])
        
        # Extend matching segments to see if whole p matches
        for m in matches:
            idx_hits += 1
            if m < start or m-start+len(p) > len(t):
                continue
            
            mismatches = 0
            
            for j in range(0, len(p)):
                if not p[j] == t[m-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break
            
            if mismatches <= n:
                all_matches.add(m - start)
    return list(all_matches), idx_hits


class SubseqIndex(object):
    """ Holds a subsequence index for a text T """

    def __init__(self, t, k, ival):
        """ Create index from all subsequences consisting of k characters
            spaced ival positions apart.  E.g., SubseqIndex("ATAT", 2, 2)
            extracts ("AA", 0) and ("TT", 1). """
        self.k = k  # num characters per subsequence extracted
        self.ival = ival  # space between them; 1=adjacent, 2=every other, etc
        self.index = []
        self.span = 1 + ival * (k - 1)
        for i in range(len(t) - self.span + 1):  # for each subseq
            self.index.append((t[i:i+self.span:ival], i))  # add (subseq, offset)
        self.index.sort()  # alphabetize by subseq

    def query(self, p):
        """ Return index hits for first subseq of p """
        subseq = p[:self.span:self.ival]  # query with first subseq
        i = bisect.bisect_left(self.index, (subseq, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != subseq:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits

if __name__ == '__main__':
    main()
