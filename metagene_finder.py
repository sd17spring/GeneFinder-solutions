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
            # Test candidate substrings of different lengths. Don't bother testing substrings
            # with length less than or equal `longest_len`, since finding one won't increase
            # the length of the longest substring. The endpoint of the range is chosen to avoid
            # running off the end of either s1 or s2.
            #
            # If we started with the *longest* possible substring length, we could break as soon as
            # we find a match. But it's slower to compare longer strings, and in general strings
            # will mismatch in the first few characters, so this approach would require both more
            # iterations and more time per average iteration step, in general.
            for substring_len in xrange(longest_len + 1, min(len(s1) + 1 - i1, len(s2) + 1 - i2)):
                if s1[i1:i1 + substring_len] == s2[i2:i2 + substring_len]:
                    longest_len = substring_len
                else:
                    # if the shorter substring doesn't match, making it longer won't going to help
                    break
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
