# -*- coding: utf-8 -*-
"""
SoftDes 2016 Mini Project 1: Gene Finder

@author: Oliver Steele

"""

from load import load_nitrogenase_seq, load_metagenome

def make_substring_set_array(string):
    """Return an array of sets of substrings of `string`, indexed by length.
    This array is suitable as an input to `find_longest_substring_length`

    Examples:
        >>> ar = make_substring_set_array('abc')
        >>> len(ar)
        4
        >>> ar[0] == set([''])
        True
        >>> ar[1] == set(['a', 'b', 'c'])
        True
        >>> ar[2] == set(['ab', 'bc'])
        True
        >>> ar[3] == set(['abc'])
        True
    """
    return [set(string[i:i + n]
            for i in xrange(len(string) + 1 - n))
            for n in xrange(1 + len(string))]

def find_longest_substring_length(string, substring_set_array):
    """Return the length of the longest substring of `string` in `substring_set_array`,
    where `substring_set_array` is a value returned by `make_substring_set_array`.

    Examples:
         >>> ar = make_substring_set_array('abc')
         >>> find_longest_substring_length('xaz', ar)
         1
         >>> find_longest_substring_length('xabz', ar)
         2
         >>> find_longest_substring_length('xabcz', ar)
         3
         >>> find_longest_substring_length('xyabcz', ar)
         3
         >>> find_longest_substring_length('abcz', ar)
         3
         >>> find_longest_substring_length('xyabc', ar)
         3
         >>> find_longest_substring_length('xyz', ar)
         0
    """

    # iterate over the possible substring lengths, in order from longest (length of the max
    # of the string length and the longest item in the substring set array).
    return next((n
                 for n in xrange(min(len(string), len(substring_set_array) - 1), 0, -1)
                 for i in xrange(len(string) + 1 - n)
                 if string[i:i + n] in substring_set_array[n]),
                0)

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
    target_array = make_substring_set_array(target_sequence)
    snippet_name, _ = max(snippets, key=lambda snippet: find_longest_substring_length(snippet[1], target_array))
    return snippet_name


def main():
    nitrogenase = load_nitrogenase_seq()
    metagenome = load_metagenome()
    print 'len', len(metagenome)
    print find_snippet_with_greatest_overlap(nitrogenase, metagenome)


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == 'test':
        import doctest
        doctest.testmod()
    else:
        main()
