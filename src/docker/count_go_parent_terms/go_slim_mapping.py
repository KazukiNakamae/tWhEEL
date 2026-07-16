"""Deterministic assignment of GO terms to one GO Slim ancestor."""


def _term_id(term):
    """Return the canonical identifier exposed by a goatools GO term."""
    return getattr(term, 'item_id', None) or getattr(term, 'id', None)


def map_to_parent(go_id, slim_ids, dag):
    """Return one deterministic nearest ``is_a`` GO Slim ancestor.

    The input term itself has distance zero. Otherwise, ancestors are searched
    breadth-first. Candidates at the first matching distance are restricted to
    the input term's namespace and resolved by ascending GO identifier.
    """
    if go_id in slim_ids:
        return go_id

    node = dag.get(go_id)
    if node is None:
        return None

    namespace = getattr(node, 'namespace', None)
    visited = {_term_id(node) or go_id}
    frontier = [node]

    while frontier:
        next_nodes = {}
        candidates = set()

        for current in frontier:
            for parent in getattr(current, 'parents', ()):
                parent_id = _term_id(parent)
                if not parent_id or parent_id in visited:
                    continue
                visited.add(parent_id)

                parent_namespace = getattr(parent, 'namespace', namespace)
                if namespace is not None and parent_namespace != namespace:
                    continue

                next_nodes[parent_id] = parent
                if parent_id in slim_ids:
                    candidates.add(parent_id)

        if candidates:
            return min(candidates)

        frontier = [next_nodes[parent_id] for parent_id in sorted(next_nodes)]

    return None


def sort_count_rows(rows):
    """Sort count rows by descending count and ascending GO identifier."""
    return sorted(rows, key=lambda row: (-row[2], row[0]))
