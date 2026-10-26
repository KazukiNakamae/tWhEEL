"""Exercise exact-count parsing and the JavaScript emitted by the generator."""

import copy
import json
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

from src.modules.report_generator import ReportGenerator


# S, M, P, and N from the retained v1.3.4 publication benchmark statistics.
BENCHMARK_COUNTS = (
    (2052, 3570, 32, 1),
    (4892, 10648, 155, 3),
    (4416, 12804, 693, 45),
    (4743, 5429, 2, 0),
    (4791, 11500, 1334, 39),
    (2215, 3094, 377, 198),
    (3401, 2505, 9030, 2180),
    (3884, 3261, 4970, 1112),
    (2663, 1571, 7671, 2241),
    (5190, 11855, 803, 10),
    (3627, 2636, 7403, 1711),
)


def statistics_text(counts):
    """Produce the statistics-file format, including rounded percentages."""
    simply, multiply, partially, not_targetable = counts
    total = sum(counts)
    fully = simply + multiply
    percent = lambda count: 100 * count / total if total else 0
    return (
        f'|F:{percent(fully):.1f}%[S:{percent(simply):.1f}%,'
        f'M:{percent(multiply):.1f}%],P:{percent(partially):.1f}%,'
        f'N:{percent(not_targetable):.1f}%,n:{total}|\n'
        f'|{fully}\tFully targetable GOs (F)\t|\n'
        f'|{simply}\tFully and simply targetable GOs (S)\t|\n'
        f'|{multiply}\tFully and multiply targetable GOs (M)\t|\n'
        f'|{partially}\tPartially targetable GOs (P)\t|\n'
        f'|{not_targetable}\tNot targetable GOs (N)\t|\n'
    )


class GOChartTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.generator = ReportGenerator({
            'output': {
                'base_dir': str(self.root),
                'intermediate': str(self.root / 'intermediate'),
            },
        })

    def parse(self, text):
        source = self.root / 'go_kfi_stat.csv'
        source.write_text(text, encoding='utf-8')
        return self.generator._parse_go_kfi_stats(source)

    def render(self, cases):
        """Run the real generated methods against a minimal SVG DOM harness."""
        if not shutil.which('node'):
            self.skipTest('Node.js is required for generated-JavaScript tests')
        self.generator._create_js_files()
        javascript = (self.root / 'static/js/main.js').read_text()
        start = javascript.index('    drawGOKFIChart() {')
        end = javascript.index('    loadGOParentImages() {', start)
        source = r'''
class Element {
    constructor(tag) {
        this.tag = tag;
        this.attrs = {};
        this.children = [];
        this.style = {};
        this.textContent = '';
        this.innerHTML = '';
    }
    setAttribute(key, value) { this.attrs[key] = String(value); }
    appendChild(child) { this.children.push(child); }
}
let container;
const document = {
    createElementNS: (_, tag) => new Element(tag),
    getElementById: () => container,
};
console.log = () => {};
console.error = () => {};
'''
        source += '\nclass Chart {\n' + javascript[start:end] + '\n}\n'
        source += 'const cases = ' + json.dumps(cases) + ';\n'
        source += r'''
const results = cases.map(item => {
    container = new Element('div');
    const chart = new Chart();
    chart.data = {statistics: {go_analysis: item.stats}};
    if (item.direct) {
        try {
            container.appendChild(chart.createGOKFIBarChart(item.stats.chart_data));
        } catch (error) { container.innerHTML = error.message; }
    } else { chart.drawGOKFIChart(); }
    return container;
});
process.stdout.write(JSON.stringify(results));
'''
        result = subprocess.run(
            [shutil.which('node'), '-'], input=source,
            capture_output=True, text=True, timeout=30,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        return json.loads(result.stdout)

    def assert_chart(self, container, counts):
        self.assertEqual(container['innerHTML'], '')
        self.assertEqual(len(container['children']), 1)
        svg = container['children'][0]
        texts = [el['textContent'] for el in svg['children'] if el['tag'] == 'text']
        simply, multiply, partially, not_targetable = counts
        total = sum(counts)
        self.assertIn(
            f'F:{simply + multiply} [S:{simply}, M:{multiply}, '
            f'P:{partially}, N:{not_targetable}, n:{total}]', texts,
        )
        self.assertIn('GO terms', texts)
        rectangles = [
            el for el in svg['children']
            if el['tag'] == 'rect' and el['attrs'].get('height') == '60'
        ]
        self.assertEqual(len(rectangles), sum(value > 0 for value in counts))
        offset = 80.0
        for rect, count in zip(rectangles, (value for value in counts if value)):
            width = float(rect['attrs']['width'])
            self.assertGreater(width, 0)
            self.assertAlmostEqual(width, 770 * count / total)
            self.assertAlmostEqual(float(rect['attrs']['x']), offset)
            offset += width
        self.assertAlmostEqual(offset, 850.0 if total else 80.0)
        serialized = json.dumps(svg)
        self.assertNotIn('NaN', serialized)
        self.assertNotIn('Infinity', serialized)

    def test_parser_preserves_exact_counts(self):
        stats = self.parse(statistics_text(BENCHMARK_COUNTS[0]))
        self.assertEqual(stats['chart_data']['fully_targetable'], {
            'count': 5622, 'simply_count': 2052, 'multiply_count': 3570,
            'percent': 99.4, 'simply_percent': 36.3, 'multiply_percent': 63.1,
        })
        self.assertEqual(stats['chart_data']['not_targetable']['count'], 1)
        self.assertEqual(stats['chart_data']['total'], 5655)

    def test_all_benchmarks_render_exact_labels_and_geometry(self):
        cases = [
            {'stats': self.parse(statistics_text(counts))}
            for counts in BENCHMARK_COUNTS
        ]
        for result, counts in zip(self.render(cases), BENCHMARK_COUNTS):
            with self.subTest(counts=counts):
                self.assert_chart(result, counts)

    def test_direct_chart_and_zero_categories(self):
        counts_cases = [(0, 0, 0, 0), (1, 0, 0, 0), (0, 0, 2, 1), (0, 2, 0, 0)]
        cases = [
            {'stats': self.parse(statistics_text(counts)), 'direct': True}
            for counts in counts_cases
        ]
        for result, counts in zip(self.render(cases), counts_cases):
            with self.subTest(counts=counts):
                self.assert_chart(result, counts)

    def test_legacy_exact_counts_work_without_new_chart_fields(self):
        stats = self.parse(statistics_text(BENCHMARK_COUNTS[0]))
        del stats['chart_data']['fully_targetable']['simply_count']
        del stats['chart_data']['fully_targetable']['multiply_count']
        self.assert_chart(self.render([{'stats': stats}])[0], BENCHMARK_COUNTS[0])

    def test_geometry_does_not_depend_on_percentages(self):
        stats = self.parse(statistics_text(BENCHMARK_COUNTS[0]))
        fully = stats['chart_data']['fully_targetable']
        fully.update(percent=0, simply_percent=100, multiply_percent=100)
        stats['chart_data']['partially_targetable']['percent'] = 100
        stats['chart_data']['not_targetable']['percent'] = 100
        self.assert_chart(self.render([{'stats': stats}])[0], BENCHMARK_COUNTS[0])

    def test_percentage_labels_are_recomputed_when_counts_are_complete(self):
        text = statistics_text(BENCHMARK_COUNTS[0]).replace('S:36.3%', 'S:0.0%')
        self.assertEqual(self.parse(text)['simply_targetable_percent'], 36.3)

    def test_missing_and_invalid_count_rows_are_not_converted_to_zero(self):
        base = statistics_text(BENCHMARK_COUNTS[0])
        row = '|2052\tFully and simply targetable GOs (S)\t|\n'
        texts = [base.replace(row, ''), base + row]
        texts.extend(base.replace('|2052\t', f'|{value}\t') for value in ['-1', '2.5', 'x'])
        stats_cases = [self.parse(text) for text in texts]
        for stats in stats_cases:
            self.assertIsNone(stats['chart_data']['fully_targetable']['simply_count'])
        for result in self.render([{'stats': stats} for stats in stats_cases]):
            self.assertIn('Exact counts are required', result['innerHTML'])
            self.assertEqual(result['children'], [])

    def test_invalid_or_conflicting_payload_counts_show_an_explicit_error(self):
        original = self.parse(statistics_text(BENCHMARK_COUNTS[0]))
        cases = []
        for value in [None, -1, 1.5, '2052', True, float('nan'), float('inf'), 9007199254740992]:
            stats = copy.deepcopy(original)
            stats['chart_data']['fully_targetable']['simply_count'] = value
            stats['simply_targetable_count'] = value
            cases.append({'stats': stats})
        for field in ['count', 'simply_count', 'multiply_count']:
            stats = copy.deepcopy(original)
            stats['chart_data']['fully_targetable'][field] += 1
            cases.append({'stats': stats})
        for field in ['simply_count', 'multiply_count']:
            stats = copy.deepcopy(original)
            del stats['chart_data']['fully_targetable'][field]
            del stats[field.replace('_count', '_targetable_count')]
            cases.append({'stats': stats})
        stats = copy.deepcopy(original)
        stats['total_gos'] += 1
        stats['chart_data']['total'] += 1
        cases.append({'stats': stats})
        stats = copy.deepcopy(original)
        stats['fully_targetable_count'] += 1
        stats['chart_data']['fully_targetable']['count'] += 1
        cases.append({'stats': stats})
        cases.extend([{'stats': {}}, {'stats': None}])
        for result in self.render(cases):
            self.assertIn('error', result['innerHTML'])
            self.assertEqual(result['children'], [])


if __name__ == '__main__':
    unittest.main()
