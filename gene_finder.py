# -*- coding: utf-8 -*-
"""
SoftDes 2017 Mini Project 1: Gene Finder

@author: Oliver Steele

"""

import random
from itertools import takewhile

from amino_acids import aa_table  # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###

START_CODON = 'ATG'
STOP_CODONS = ['TAG', 'TAA', 'TGA']

# Initialize the dictionary. This is obfuscated hackery for such a short dictionary.
# `{'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}` would suffice.
NUCLEOTIDE_COMPLEMENTS = {n1: n2 for i in [-1, 1] for n1, n2 in zip(*['AC', 'TG'][::i])}
"""A dict of {nucleotide: complement}"""


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    >>> get_complement('T')
    'A'
    >>> get_complement('G')
    'C'
    """
    return NUCLEOTIDE_COMPLEMENTS[nucleotide]


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    return ''.join(get_complement(n) for n in dna[::-1])


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    >>> rest_of_ORF("ATGAG") # return whole string, even to incomplete codon
    'ATGAG'
    """
    codons = (dna[i:i + 3] for i in range(0, len(dna), 3))
    return ''.join(takewhile(lambda codon: codon not in STOP_CODONS, codons))


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    >>> find_all_ORFs_oneframe("AATGCATTAG") # non-aligned
    []
    >>> find_all_ORFs_oneframe("ATGATGTAGATGAAATAG") # nested
    ['ATGATG', 'ATGAAA']
    """
    frames = []
    i = 0
    while i + 2 < len(dna):
        if dna[i:i + 3] == START_CODON:
            frame = rest_of_ORF(dna[i:])
            frames.append(frame)
            i += len(frame)
        else:
            i += 3
    return frames


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    >>> find_all_ORFs("ATGCA")
    ['ATGCA']
    """
    return [orf for i in range(3) for orf in find_all_ORFs_oneframe(dna[i:])]


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    return find_all_ORFs(dna) + find_all_ORFs(get_reverse_complement(dna))


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("TAG") # no ORF
    >>> longest_ORF("ATGCCCTGAATGTAG") # longest ORF is first
    'ATGCCC'
    >>> longest_ORF("ATGTAGATGCCCTGA") # longest ORF is second
    'ATGCCC'
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    orfs = find_all_ORFs_both_strands(dna)
    return max(orfs, key=len) if orfs else None


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        >>> longest_ORF_noncoding("ATGCGAATGTAGCATCAAA", 100) > 10
        True

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    return max(len(longest_ORF(shuffle_string(dna)) or '') for _ in range(num_trials))


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    return ''.join(aa_table[dna[i:i + 3]] for i in range(0, len(dna) - 2, 3))


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    threshold = longest_ORF_noncoding(dna, 1500)
    return [coding_strand_to_AA(orf)
            for orf in find_all_ORFs_both_strands(dna)
            if len(orf) > threshold]


if __name__ == "__main__":
    import sys
    import doctest

    if sys.argv[1:]:
        for arg in sys.argv[1:]:
            if arg.endswith('.fa'):
                dna = load_seq(arg)
                print('Finding genes in {}-nucleotide strand in file {}...'.format(len(dna), arg))
                print(gene_finder(dna))
            else:
                print('test', arg)
                doctest.run_docstring_examples(globals()[arg], globals())
    else:
        doctest.testmod()
