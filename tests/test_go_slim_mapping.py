import unittest

from src.docker.count_go_parent_terms.go_slim_mapping import (
    map_to_parent,
    sort_count_rows,
)


class FakeGoTerm:
    def __init__(self, go_id, namespace='biological_process'):
        self.item_id = go_id
        self.namespace = namespace
        self.parents = []


class GoSlimMappingTests(unittest.TestCase):
    @staticmethod
    def _term(go_id, *parents, namespace='biological_process'):
        term = FakeGoTerm(go_id, namespace=namespace)
        term.parents = list(parents)
        return term

    def test_slim_term_maps_to_itself(self):
        term = self._term('GO:0000005')

        self.assertEqual(
            map_to_parent(term.item_id, {term.item_id}, {term.item_id: term}),
            term.item_id,
        )

    def test_nearest_slim_ancestor_wins_over_lower_go_id(self):
        far_slim = self._term('GO:0000001')
        intermediate = self._term('GO:0000009', far_slim)
        near_slim = self._term('GO:0000008')
        leaf = self._term('GO:0000010', intermediate, near_slim)
        dag = {term.item_id: term for term in (
            leaf, intermediate, near_slim, far_slim
        )}

        self.assertEqual(
            map_to_parent(
                leaf.item_id,
                {far_slim.item_id, near_slim.item_id},
                dag,
            ),
            near_slim.item_id,
        )

    def test_equal_distance_tie_uses_go_id_for_both_parent_orders(self):
        lower_id = self._term('GO:0000002')
        higher_id = self._term('GO:0000007')
        slim_ids = {higher_id.item_id, lower_id.item_id}

        for parent_order in (
            [higher_id, lower_id],
            [lower_id, higher_id],
        ):
            with self.subTest(parent_order=[p.item_id for p in parent_order]):
                leaf = self._term('GO:0000010', *parent_order)
                dag = {term.item_id: term for term in (
                    leaf, lower_id, higher_id
                )}
                self.assertEqual(
                    map_to_parent(leaf.item_id, slim_ids, dag),
                    lower_id.item_id,
                )

    def test_cross_namespace_slim_is_not_selected(self):
        wrong_namespace = self._term(
            'GO:0000001', namespace='molecular_function'
        )
        valid_slim = self._term('GO:0000004')
        intermediate = self._term('GO:0000008', valid_slim)
        leaf = self._term('GO:0000010', wrong_namespace, intermediate)
        dag = {term.item_id: term for term in (
            leaf, wrong_namespace, intermediate, valid_slim
        )}

        self.assertEqual(
            map_to_parent(
                leaf.item_id,
                {wrong_namespace.item_id, valid_slim.item_id},
                dag,
            ),
            valid_slim.item_id,
        )

    def test_unknown_or_unmapped_term_returns_none(self):
        root = self._term('GO:0000001')
        leaf = self._term('GO:0000010', root)

        self.assertIsNone(map_to_parent('GO:9999999', set(), {}))
        self.assertIsNone(
            map_to_parent(leaf.item_id, {'GO:0000002'}, {
                leaf.item_id: leaf,
                root.item_id: root,
            })
        )

    def test_count_rows_have_deterministic_tie_order(self):
        rows = [
            ('GO:0000007', 'seven', 2),
            ('GO:0000009', 'nine', 3),
            ('GO:0000002', 'two', 2),
        ]

        self.assertEqual(sort_count_rows(rows), [
            ('GO:0000009', 'nine', 3),
            ('GO:0000002', 'two', 2),
            ('GO:0000007', 'seven', 2),
        ])


if __name__ == '__main__':
    unittest.main()
