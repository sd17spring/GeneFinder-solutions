# -*- coding: utf-8 -*-
"""
Functions to create and traverse suffix trees.

This is the minimal implementation necessary for the "Going Beyond" section of the
`Gene Finder mini-project <https://sites.google.com/site/sd16spring/home/assignments-and-mini-projects/gene-finder>`.

This implementation studiosly avoids classes and tuples because they haven't been covered yet.

Strings are terminated with an digit character e.g. '0', instead of a string e.g. '$0'.
This simplifies the implementation and is sufficient for this problem domain.

References:
    * https://en.wikipedia.org/wiki/Generalized_suffix_tree
"""

# An internal node of a suffix tree is represented as a dictionary {letter: node}.
# A leaf node of a suffix tree is represented as a list of [string number, offset].

def empty_suffix_tree():
    """Create a new empty suffix tree."""
    return {}

def is_leaf(node):
    """Returns true iff a node of a suffix tree is a leaf (terminal node).

    Examples:
        >>> is_leaf(empty_suffix_tree())
        False
    """
    return not isinstance(node, dict)

def is_nonterminal(node):
    """Returns true iff a node of a suffix tree is a non-terminal node.

    Examples:
        >>> is_nonterminal(empty_suffix_tree())
        True
    """
    return isinstance(node, dict)

def get_node_at_substring(suffix_tree, substring):
    """Returns the node or metadata at `substring`.

    Examples:
        >>> bool(get_node_at_substring(empty_suffix_tree(), 'A'))
        False
    """
    node = suffix_tree
    for c in substring:
        node = node.get(c, None)
        if not node:
            break
    return node

def has_suffix(suffix_tree, suffix):
    """Returns true iff `suffix_tree` contains `suffix`.

    Examples:
        >>> has_suffix(empty_suffix_tree(), 'A')
        False
    """
    node = get_node_at_substring(suffix_tree, suffix)
    return bool(node and any(key.isdigit() for key in node.keys()))

def add_terminated_suffix(suffix_tree, suffix, metadata):
    """Modify `suffix_tree` to include `suffix`, terminated by `metadata`.

    Examples:
        >>> st = empty_suffix_tree()
        >>> add_terminated_suffix(st, 'A0', 'd0')
        >>> has_suffix(st, 'A')
        True
        >>> has_suffix(st, 'B')
        False

        >>> add_terminated_suffix(st, 'AB0', 'd1')
        >>> has_suffix(st, 'B')
        True
        >>> has_suffix(st, 'AB')
        True

        >>> print st
        {'A': {'0': 'd0', 'B': {'0': 'd1'}}}

        >>> add_terminated_suffix(st, 'A1', 'd2')
        >>> print st
        {'A': {'1': 'd2', '0': 'd0', 'B': {'0': 'd1'}}}
    """
    node = suffix_tree
    for i, c in enumerate(suffix):
        # print 'add', node, c
        child = node.get(c)
        if not child:
            child = empty_suffix_tree() if i + 1 < len(suffix) else metadata
            node[c] = child
        node = child

def add_string_suffixes(suffix_tree, string, string_number):
    """Modify `suffix_tree` to include all the suffixes of `string` (including `string` in its entirety).
    The metadata for each suffix is a list `[string_number, suffix_offset]`, where `suffix_offset` is
    the start position of the suffix within `string`.

    Examples:
        >>> st = empty_suffix_tree()
        >>> add_string_suffixes(st, 'ABA', 0)
        >>> print st
        {'A': {'0': [0, 2], 'B': {'A': {'0': [0, 0]}}}, 'B': {'A': {'0': [0, 1]}}}

        >>> st = empty_suffix_tree()
        >>> add_string_suffixes(st, 'ABAB', 0)
        >>> print st
        {'A': {'B': {'A': {'B': {'0': [0, 0]}}, '0': [0, 2]}}, 'B': {'A': {'B': {'0': [0, 1]}}, '0': [0, 3]}}
        >>> add_string_suffixes(st, 'BABA1', 1)
    """
    assert 0 <= string_number < 10
    terminated_string = string + str(string_number)
    for suffix_offset in xrange(len(string)):
        add_terminated_suffix(suffix_tree, terminated_string[suffix_offset:], [string_number, suffix_offset])

def add_strings_to_suffix_tree(suffix_tree, strings):
    """Add each of the strings in the iterable `strings` to `suffix_tree`.

        >>> st = empty_suffix_tree()
        >>> add_strings_to_suffix_tree(st, ['ABAB', 'BABA', 'ABBA'])
    """
    for string_number, string in enumerate(strings):
        add_string_suffixes(suffix_tree, string, string_number)

def generalized_suffix_tree(strings):
    st = empty_suffix_tree()
    add_strings_to_suffix_tree(st, strings)
    return st

def depth_first(suffix_tree, nodefn, leaffn):
    def traverse(node, path):
        if is_leaf(node):
            return leaffn(path, node)
        else:
            return nodefn(path, [traverse(child, path + c) for c, child in node.items()])
    return traverse(suffix_tree, '')


if __name__ == "__main__":
    import doctest
    doctest.testmod()
    st = empty_suffix_tree()
