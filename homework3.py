#!/usr/bin/env python

"""Algorithms for DNA Sequencing - Programming Homework 3"""

from collections import defaultdict


def main():
    chr1 = readGenome('chr1.GRCh38.excerpt.fasta')
    p = 'GCTGATCGATCGTACG'
    print(approximate_match(p, chr1))
    p = 'GATTTACCAGATTGAG'
    print(approximate_match(p, chr1))
    seqs, quals = readFastq('ERR266411_1.for_asm.fastq')
    edges, suffixes = overlap_graph(seqs, 30)
    print(edges)
    print(suffixes)


def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome


def approximate_match(p, t):
    # Create distance matrix
    D = []
    for i in range(len(p)+1):
        D.append([0]*(len(t)+1))
       
    # Initialize first row and column of matrix
    for i in range(len(p)+1):
        D[i][0] = i
    for i in range(len(t)+1):
        D[0][i] = 0
       
    # Fill in the rest of the matrix
    for i in range(1, len(p)+1):
        for j in range(1, len(t)+1):
            distHor = D[i][j-1] + 1
            distVer = D[i-1][j] + 1
            if p[i-1] == t[j-1]:
                distDiag = D[i-1][j-1]
            else:
                distDiag = D[i-1][j-1] + 1
            D[i][j] = min(distHor, distVer, distDiag)
       
    # Edit distance is the value in the bottom right corner of the matrix
    return min(D[-1])


def readFastq(filename):
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline()  # skip name line
            seq = fh.readline().rstrip()  # read base sequence
            fh.readline()  # skip placeholder line
            qual = fh.readline().rstrip()  # base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities


def overlap_graph(reads, k):
    # Make index
    index = defaultdict(set)
    for read in reads:
        for i in range(len(read) - k + 1):
            index[read[i:i+k]].add(read)

    # Make graph
    graph = defaultdict(set)
    for r in reads:
        for o in index[r[-k:]]:
            if r != o:
                if overlap(r, o, k):
                    graph[r].add(o)

    edges = 0
    for read in graph:
        edges += len(graph[read])
    return(edges, len(graph))


def overlap(a, b, min_length=3):
    """ Return length of longest suffix of 'a' matching
        a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists,
        return 0. """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's suffx in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match

if __name__ == '__main__':
    main()
