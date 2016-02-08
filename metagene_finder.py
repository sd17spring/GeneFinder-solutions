# -*- coding: utf-8 -*-
"""
SoftDes 2016 Mini Project 1: Gene Finder

@author: Oliver Steele

"""

from load import load_nitrogenase_seq, load_metagenome

def find_longest_common_substring_length(s1, s2):
    """Return the length of the longest substring of `string` in `substring_set_array`,
    where `substring_set_array` is a value returned by `make_substring_set_array`.

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

    L = {}
    for i, c1 in enumerate(s1):
        for j, c2 in enumerate(s2):
            if c1 == c2:
                L[i, j] = L.get((i - 1, j - 1), 0) + 1
    return max(L.values() or [0])

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
