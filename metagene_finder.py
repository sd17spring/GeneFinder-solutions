# -*- coding: utf-8 -*-
"""
SoftDes 2016 Mini Project 1: Gene Finder

@author: Oliver Steele

"""

from load import load_nitrogenase_seq, load_metagenome
from suffix_tree import generalized_suffix_tree, depth_first

def common_substrings(s1, s2):
    """
    Examples:
        >>> sorted(common_substrings("ABAB", "BABA"))
        ['A', 'AB', 'ABA', 'B', 'BA', 'BAB']
        >>> sorted(common_substrings("ABAB", "ABBA"))
        ['A', 'AB', 'B', 'BA']
        >>> sorted(common_substrings("BABA", "ABBA"))
        ['A', 'AB', 'B', 'BA']
    """
    prefixes = []

    def leaffn(path, leaf):
        return set([leaf[0]])

    def nodefn(path, sets):
        union = set(i for s in sets for i in s)
        if len(union) == 2:
            prefixes.append(path)
        return union

    depth_first(generalized_suffix_tree([s1, s2]), nodefn, leaffn)
    return filter(None, prefixes)


def longest_common_substring(s1, s2):
    """
    Examples:
        >>> longest_common_substring("ABAB", "BABA")
        'ABA'
        >>> longest_common_substring("ABAB", "ABBA")
        'AB'
        >>> longest_common_substring("BABA", "ABBA")
        'AB'
    """
    # `sorted` makes deterministic for testing
    return max(sorted(common_substrings(s1, s2)), key=len)

def main():
    nitrogenase = load_nitrogenase_seq()
    metagenome = load_metagenome()
    snippet, _ = max(metagenome, key=lambda pair: longest_common_substring(nitrogenase, pair[1]))
    print snippet


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == 'test':
        import doctest
        doctest.testmod()
    else:
        main()
