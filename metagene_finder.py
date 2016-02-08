# -*- coding: utf-8 -*-
"""
SoftDes 2016 Mini Project 1: Gene Finder

@author: Oliver Steele

"""

from load import load_nitrogenase_seq, load_metagenome

def make_substring_tree(string):
    """Return a tree of substrings of string. The tree is represented as a rooted graph
    whose nodes are dicts and whose keys are letters of the string.

    Examples:
        >>> make_substring_tree('abc')
        {'a': {'b': {'c': {}}}, 'c': {}, 'b': {'c': {}}}
    """
    def extend_graph(node, string):
        for edge in string:
            child = node.get(edge)
            if not child:
                node[edge] = child = {}
            node = child

    root = {}
    for i in xrange(len(string)):
        extend_graph(root, string[i:])
    return root

def find_longest_substring_length(string, substring_tree):
    """Return the length of the longest substring of `string` in `substring_set_array`,
    where `substring_set_array` is a value returned by `make_substring_set_array`.

    Examples:
        >>> t = make_substring_tree('abc')
        >>> find_longest_substring_length('xaz', t)
        1
        >>> find_longest_substring_length('xabz', t)
        2
        >>> find_longest_substring_length('xabcz', t)
        3
        >>> find_longest_substring_length('xyabcz', t)
        3
        >>> find_longest_substring_length('abcz', t)
        3
        >>> find_longest_substring_length('xyabc', t)
        3
        >>> find_longest_substring_length('xyz', t)
        0
    """

    longest = 0
    for i in xrange(len(string)):
        if longest > len(string) - i + 1:
            break
        node = substring_tree
        for j in xrange(i, len(string)):
            node = node.get(string[j])
            if node is None:
                break
            longest = max(longest, j + 1 - i)
    return longest

def find_snippet_with_greatest_overlap(target_sequence, snippets):
    """Return the name of the snippet whose sequence has the greatest overlap with `target_sequence`.

    Args:
        snippets (list): a list of tuples of (snippet_name, sequence)
        target_sequence (str): the nucleotide sequence

    Examples:
        >>> t = make_substring_tree('TAG')
        >>> snippets = [('1', 'AG'), ('2', 'AGAGAG'), ('3', 'ATAGA')]
        >>> find_snippet_with_greatest_overlap('TAG', snippets)
        '3'
        >>> find_snippet_with_greatest_overlap('AGAG', snippets)
        '2'
    """
    substring_tree = make_substring_tree(target_sequence)
    snippet_name, _ = max(snippets, key=lambda snippet: find_longest_substring_length(snippet[1], substring_tree))
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
