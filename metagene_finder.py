# -*- coding: utf-8 -*-
"""
SoftDes 2016 Mini Project 1: Gene Finder

@author: Oliver Steele

"""

from load import load_nitrogenase_seq, load_metagenome

def find_longest_common_substring_length(s1, s2):
    """Return the length of the longest common substring of `s1` and `s2`.

    Examples:
        >>> find_longest_common_substring_length('abc', 'xaz')
        1
        >>> find_longest_common_substring_length('abc', 'xabz')
        2
        >>> find_longest_common_substring_length('abc', 'xabcz')
        3
        >>> find_longest_common_substring_length('abc', 'xyabcz')
        3
        >>> find_longest_common_substring_length('abc', 'abcz')
        3
        >>> find_longest_common_substring_length('abc', 'xyabc')
        3
        >>> find_longest_common_substring_length('abc', 'xyz')
        0
    """

    longest_len = 0
    for i1 in xrange(len(s1)):
        for i2 in xrange(len(s2)):
            offset = 0
            while i1 + offset < len(s1) and i2 + offset < len(s2) and s1[i1 + offset] == s2[i2 + offset]:
                offset += 1
            longest_len = max(longest_len, offset)
    return longest_len

def find_snippet_with_greatest_overlap(target_sequence, snippets):
    """Return the name of the snippet whose sequence has the greatest overlap with `target_sequence`.

    Args:
        snippets (list): a list of tuples of (snippet_name, sequence)
        target_sequence (str): the nucleotide sequence

    Examples:
        >>> snippets = [('1', 'AG'), ('2', 'AGAGAG'), ('3', 'ATAGA')]
        >>> find_snippet_with_greatest_overlap('TAG', snippets)
        '3'
        >>> find_snippet_with_greatest_overlap('AGAG', snippets)
        '2'
    """
    snippet_name, _ = max(snippets, key=lambda snippet: find_longest_common_substring_length(snippet[1], target_sequence))
    return snippet_name


def main():
    nitrogenase = load_nitrogenase_seq()
    metagenome = load_metagenome()
    print find_snippet_with_greatest_overlap(nitrogenase, metagenome)


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == 'test':
        import doctest
        doctest.testmod()
    else:
        main()
