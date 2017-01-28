# -*- coding: utf-8 -*-
"""
SoftDes 2016 Mini Project 1: Gene Finder

@author: Oliver Steele

"""

from multiprocessing import Pool

from load import load_metagenome, load_nitrogenase_seq


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
    for i in range(len(string)):
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
    for i in range(len(string)):
        if longest > len(string) - i + 1:
            break
        node = substring_tree
        for j in range(i, len(string)):
            node = node.get(string[j])
            if node is None:
                break
            longest = max(longest, j + 1 - i)
    return longest


def p_find_longest_substring_length(sequence):
    """This is an adaptor to make `find_longest_substring_length` callable
    by Pool().map, even though the substring tree is too deep to serialize.
    """
    global g_substring_tree
    return find_longest_substring_length(sequence, g_substring_tree)


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
    # Set a global variable to the target sequence *before* calling Pool().
    # This insures that the value is present within each process that is forked from this one.
    # This wouldn't work if this function were ever called with a different
    # target_sequence, since then the child processes would already have the previous value.
    global g_substring_tree
    g_substring_tree = make_substring_tree(target_sequence)

    p = Pool(12)
    longest_substring_lengths = p.map(p_find_longest_substring_length, (snippet[1] for snippet in snippets))
    max_index = longest_substring_lengths.index(max(longest_substring_lengths))
    max_snippet = snippets[max_index]
    return max_snippet[0]


def main():
    nitrogenase = load_nitrogenase_seq()
    metagenome = load_metagenome()
    print(find_snippet_with_greatest_overlap(nitrogenase, metagenome))


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == 'test':
        import doctest
        doctest.testmod()
    else:
        main()
