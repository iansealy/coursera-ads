#!/usr/bin/env python

"""Algorithms for DNA Sequencing - Programming Homework 1"""


def main():
    virus = readGenome('lambda_virus.fa')
    print(len(naive_with_rc('AGGT', virus)))
    print(len(naive_with_rc('TTAA', virus)))
    print(naive_with_rc('ACTAAGT', virus)[0])
    print(naive_with_rc('AGTCGA', virus)[0])
    print(len(naive_2mm('TTCAAGCC', virus)))
    print(naive_2mm('AGGAGGTT', virus)[0])
    sequences, qualities = readFastq('ERR037900_1.first1000.fastq')
    print(lowest_quality_base(qualities))


def naive(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences


def reverseComplement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t


def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome


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


def naive_with_rc(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        for seq in (p, reverseComplement(p)):
            match = True
            for j in range(len(seq)):  # loop over characters
                if t[i+j] != seq[j]:  # compare characters
                    match = False
                    break
            if match:
                occurrences.append(i)  # all chars matched; record
                break
    return occurrences


def naive_2mm(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        mismatches = 0
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                mismatches += 1
                if mismatches > 2:
                    break
        if mismatches <= 2:
            occurrences.append(i)  # all chars matched; record
    return occurrences


def lowest_quality_base(qs):
    total = [0] * len(qs[0])
    for q in qs:
        for i, phred in enumerate(q):
            total[i] += phred33ToQ(phred)
    return total.index(min(total))


def phred33ToQ(qual):
    return ord(qual) - 33

if __name__ == '__main__':
    main()
