# -*- coding: utf-8 -*-
"""
SoftDes 2017 Mini Project 1: Gene Finder

@author: Oliver Steele

"""

import random
import amino_acids   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way.
    >>> random.seed(1)
    >>> shuffle_string('abcdef')
    'aedfcb'
    >>> shuffle_string('abcdef')
    'dfaebc'
    """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


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
    for complement_pair in [['A', 'T'], ['C', 'G']]:
        n1 = complement_pair[0]
        n2 = complement_pair[1]
        if nucleotide == n1:
            return n2
        if nucleotide == n2:
            return n1
    1/0  #  raise an error, to make it instantly evident that something broke.
    # The description of the error will be misleading, but the stack trace will take you to this line
    # and also show who called the function.

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
    complements = ''
    for n in dna:
        complements = get_complement(n) + complements
    return complements

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
    stop_codons = ['TAG', 'TAA', 'TGA']
    sequence = ''
    for i in range(0, len(dna), 3):
        codon = dna[i:i + 3]
        if codon in stop_codons:
            break
        sequence = sequence + codon
    return sequence


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
    i = 0
    frames = []
    start_codon = 'ATG'
    while i < len(dna):
        if dna[i:i + 3] == start_codon:
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
    orfs = []
    for i in xrange(3):
        orfs.extend(find_all_ORFs_oneframe(dna[i:]))
    return orfs


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
    orfs.sort(key=len, reverse=True)
    if orfs:
        return orfs[0]
    else:
        return None


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        >>> random.seed(1)
        >>> longest_ORF_noncoding("ATGCGAATGTAGCATCAAA", 100)
        19

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    max_len = 0
    for _ in xrange(num_trials):
        orf = longest_ORF(shuffle_string(dna))
        if orf:
            max_len = max(max_len, len(orf))
    return max_len


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
    aa = ''
    for i in xrange(0, len(dna) - 2, 3):
        codon = dna[i:i + 3]
        aa_index = None
        for j in range(len(amino_acids.codons)):
            if codon in amino_acids.codons[j]:
                aa_index = j
                break
        aa = aa + amino_acids.aa[aa_index]
    return aa


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        >>> random.seed(1)
        >>> gene_finder("ATGTCATTGCGTGTGAGACAGATTGATCGTCGCGAATGGCTATTGGCGCAAACCGCGACAGAATGCCAGCGCCATGGCCGGGAA" \
                        "GCGACGCTGGAATATCCGACGCGACAGGGAATGTGGGTTCGGTTGAGCGATGCAGAAAAACGGTGGTCGGCCTGGATTAAACCT" \
                        "GGGGACTGGCTTGAGCATGTCTCTCCCGCTCTGGCTGGGGCGGCGGTTTCTGCTGGCGCTGAGCACCTGGTCGTTCCCTGGCTT")
        ['MSLRVRQIDRREWLLAQTATECQRHGREATLEYPTRQGMWVRLSDAEKRWSAWIKPGDWLEHVSPALAGAAVSAGAEHLVVPWL']

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    threshold = longest_ORF_noncoding(dna, 1500)
    aas = []
    for orf in find_all_ORFs_both_strands(dna):
        if len(orf) > threshold:
            aas.append(coding_strand_to_AA(orf))
    return aas

if __name__ == "__main__":
    import sys
    import doctest

    if sys.argv[1:]:
        for arg in sys.argv[1:]:
            if arg.endswith('.fa'):
                dna = load_seq(arg)
                print 'Finding genes in {}-nucleotide strand in file {}...'.format(len(dna), arg)
                print gene_finder(dna)
            else:
                print 'test', arg
                doctest.run_docstring_examples(globals()[arg], globals())
    else:
        doctest.testmod()
