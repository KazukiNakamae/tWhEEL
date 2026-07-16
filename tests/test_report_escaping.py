import csv
import io
import json
import re
import shutil
import subprocess
import tempfile
import unittest
import xml.etree.ElementTree as ET
from html.parser import HTMLParser
from pathlib import Path

from src.modules.pathway_report_generator import PathwayReportGenerator
from src.modules.report_generator import (
    ReportGenerator,
    serialize_json_for_html_script,
)


DANGEROUS_TEXT = (
    '</script><img src=x onerror="owned"> & "\' \u2028\u2029'
)


class ReportHTMLInspector(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.tags = []
        self.event_attributes = []
        self.text = []
        self.titles = []

    def handle_starttag(self, tag, attrs):
        self.tags.append(tag)
        for name, value in attrs:
            if name.lower().startswith('on'):
                self.event_attributes.append((tag, name, value))
            if name.lower() == 'title':
                self.titles.append(value)

    def handle_data(self, data):
        self.text.append(data)


class ReportEscapingTests(unittest.TestCase):
    @staticmethod
    def config(root, dangerous=DANGEROUS_TEXT):
        output = root / 'output'
        intermediate = output / 'intermediate'
        intermediate.mkdir(parents=True)
        return {
            'input': {
                'genome_assembly': dangerous,
                'transcriptome': dangerous,
                'protein_db_for_goanalysis': dangerous,
                'uniprot_db': dangerous,
            },
            'crispr': {
                'pam_pattern': dangerous,
                'max_mismatch': dangerous,
                'threads': dangerous,
            },
            'busco': {
                'lineage': dangerous,
                'cpu': dangerous,
            },
            'database': {'force_rebuild': False},
            'report': {
                'generate_summary': False,
                'generate_html': True,
            },
            'output': {
                'base_dir': str(output),
                'intermediate': str(intermediate),
            },
            'pathway_analysis': {},
            'pathway_report': {},
        }

    @staticmethod
    def standard_results(config, final_path):
        return {
            'timestamp': DANGEROUS_TEXT,
            'statistics': {
                'gene_prediction': {'protein_count': 0},
                'annotation': {
                    'multi_genes': 0,
                    'single_genes': 1,
                    'total_genes': 1,
                },
                'guide_design': {},
                'go_analysis': {},
            },
            'config': config,
            'data': {'final_annotation_guides': str(final_path)},
        }

    @staticmethod
    def extract_json_script(document, element_id):
        marker = f'<script type="application/json" id="{element_id}">'
        return json.loads(document.split(marker, 1)[1].split('</script>', 1)[0])

    @staticmethod
    def extract_javascript_variable(source, variable_name):
        marker = f'var {variable_name} = '
        value, _ = json.JSONDecoder().raw_decode(
            source[source.index(marker) + len(marker):]
        )
        return value

    def test_html_script_json_serializer_round_trips_dangerous_text(self):
        serialized = serialize_json_for_html_script({
            'value': DANGEROUS_TEXT,
        })

        self.assertNotIn('<', serialized)
        self.assertNotIn('>', serialized)
        self.assertNotIn('&', serialized)
        self.assertNotIn('\u2028', serialized.replace('\\u2028', ''))
        self.assertNotIn('\u2029', serialized.replace('\\u2029', ''))
        self.assertEqual(json.loads(serialized), {'value': DANGEROUS_TEXT})

    def test_standard_report_escapes_json_html_and_static_table(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            root = Path(temp_dir_name)
            config = self.config(root)
            final_path = root / 'annotation_with_guides_final.tsv'
            final_path.write_text(
                '\t'.join([
                    'query_id', 'subject_id', 'identity',
                    'alignment_length', 'evalue', 'bit_score',
                    'description', 'gene_name', 'classification',
                    'unique_id', 'guide_targets',
                ])
                + '\n'
                + '\t'.join([
                    DANGEROUS_TEXT, DANGEROUS_TEXT, '99', '10', '0', '1',
                    DANGEROUS_TEXT, DANGEROUS_TEXT, 'Single',
                    DANGEROUS_TEXT, '',
                ])
                + '\n',
                encoding='utf-8',
            )
            generator = ReportGenerator(config)
            results = self.standard_results(config, final_path)

            published = generator.publish_browser_payload(results)
            interactive_html = generator._generate_html_content(
                results,
                payload_script_src=published['core_relative_path'],
            )
            payload_source = Path(published['core_path']).read_text(
                encoding='utf-8'
            )
            payload = self.extract_javascript_variable(
                payload_source, 'payload'
            )['data']
            static_html = generator._generate_static_html_content(results)

        self.assertEqual(
            payload['config']['input']['genome_assembly'], DANGEROUS_TEXT
        )
        self.assertEqual(
            payload['gene_annotation'][0]['gene_name'], DANGEROUS_TEXT
        )
        self.assertNotIn('<img src=x onerror="owned">', interactive_html)
        self.assertNotIn('<img src=x onerror="owned">', static_html)
        self.assertNotIn('twheel-data', interactive_html)
        self.assertIn('\\u003c/script\\u003e', payload_source)
        self.assertIn('&lt;/script&gt;', static_html)
        self.assertIn('&quot;', static_html)
        self.assertIn('&#x27;', static_html)

        inspector = ReportHTMLInspector()
        inspector.feed(static_html)
        self.assertEqual(inspector.event_attributes, [])
        self.assertNotIn('img', inspector.tags)
        self.assertIn(DANGEROUS_TEXT, inspector.titles)
        self.assertIn(DANGEROUS_TEXT, ''.join(inspector.text))

    def test_standard_javascript_escapes_dynamic_report_values(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            root = Path(temp_dir_name)
            generator = ReportGenerator(self.config(root, 'safe'))
            generator._create_js_files()
            javascript = (
                Path(generator.output_dir) / 'static' / 'js' / 'main.js'
            ).read_text(encoding='utf-8')

        for unsafe_expression in (
            '${gene.query_id}',
            '${gene.subject_id}',
            '${gene.gene_name}',
            '${guide.target_sequence}',
            '${go.go_name || \'Loading...\'}',
            '${protein.id}',
        ):
            self.assertNotIn(unsafe_expression, javascript)
        self.assertIn('escapeHTML(value)', javascript)
        self.assertIn('data-go-download-id=', javascript)
        self.assertIn('${this.escapeHTML(parent.description)}', javascript)
        self.assertLess(
            javascript.index(
                'delete window.__TWHEEL_REPORT_GUIDE_CHUNKS__[registryKey]'
            ),
            javascript.index(
                'if (!Array.isArray(rows) || rows.length !== chunk.row_count)'
            ),
        )
        self.assertNotIn(
            'onclick="window.tWheelReport.downloadGOGuideRNAs', javascript
        )
        self.assertIn('encodeURIComponent(go.go_id)', javascript)
        self.assertIn("document.createTextNode(` ${String(message ?? '')}`)", javascript)

    @unittest.skipUnless(shutil.which('node'), 'Node.js is not available.')
    def test_standard_csv_downloads_escape_quotes_and_commas(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            root = Path(temp_dir_name)
            generator = ReportGenerator(self.config(root, 'safe'))
            generator._create_js_files()
            javascript = (
                Path(generator.output_dir) / 'static' / 'js' / 'main.js'
            ).read_text(encoding='utf-8')

        helper_start = javascript.index('    escapeCSVField(value) {')
        helper_end = javascript.index(
            '    // Guide RNA CSV download functionality', helper_start
        )
        guide_start = javascript.index(
            '    generateGOGuideRNACSV(guides, goId, goName) {'
        )
        guide_end = javascript.index(
            '    downloadCSV(csvContent, filename) {', guide_start
        )
        parent_start = javascript.index(
            '    downloadGOParentTermsCSV() {'
        )
        parent_end = javascript.index(
            '    // Copy to clipboard functionality', parent_start
        )

        go_name = 'Alpha "quoted", term'
        query_id = 'query "quoted", 1'
        description = 'regulation of "alpha, beta" process'
        guide = {
            'matched_gene': 'Single_0000000001',
            'query_id': query_id,
            'reference': 'Single_0000000001',
            'start': 1,
            'end': 23,
            'target_sequence': 'ACGT',
            'strand': '+',
            'mm0': 0,
            'mm1_2': 0,
            'off_target_score': 100,
        }
        node_source = (
            'const outputs = [];\n'
            'class CSVExporter {\n'
            + javascript[helper_start:helper_end]
            + javascript[guide_start:guide_end]
            + javascript[parent_start:parent_end]
            + '    downloadCSV(content) { outputs.push(content); }\n'
            + '    showSuccess() {}\n'
            + '    showError(message) { throw new Error(message); }\n'
            + '}\n'
            + 'console.log = () => {};\n'
            + 'const exporter = new CSVExporter();\n'
            + 'outputs.push(exporter.generateGOGuideRNACSV('
            + json.dumps([guide])
            + ", 'GO:0000001', "
            + json.dumps(go_name)
            + '));\n'
            + 'exporter.data = {go_parent_counts: ['
            + json.dumps({
                'parent_go': 'GO:0000002',
                'description': description,
                'count': 1,
            })
            + ']};\n'
            + 'exporter.downloadGOParentTermsCSV();\n'
            + 'process.stdout.write(JSON.stringify(outputs));\n'
        )
        result = subprocess.run(
            [shutil.which('node'), '-e', node_source],
            check=True,
            capture_output=True,
            text=True,
        )
        guide_csv, parent_csv = json.loads(result.stdout)
        guide_rows = list(csv.reader(
            io.StringIO(guide_csv, newline=''), strict=True
        ))
        parent_rows = list(csv.reader(
            io.StringIO(parent_csv, newline=''), strict=True
        ))

        self.assertEqual(len(guide_rows[0]), 12)
        self.assertEqual(len(guide_rows[1]), 12)
        self.assertEqual(guide_rows[1][1], go_name)
        self.assertEqual(guide_rows[1][3], query_id)
        self.assertEqual(len(parent_rows[0]), 3)
        self.assertEqual(len(parent_rows[1]), 3)
        self.assertEqual(parent_rows[1][1], description)

    def test_pathway_json_and_dynamic_rendering_contract_are_safe(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            root = Path(temp_dir_name)
            config = self.config(root)
            generator = PathwayReportGenerator(config)
            summary = {
                'timestamp': DANGEROUS_TEXT,
                'statistics': {
                    'pathway_analysis': {
                        'total_pathways': 1,
                        'fully_targetable_pathways': 0,
                        'partially_targetable_pathways': 0,
                        'not_targetable_pathways': 1,
                        'chart_data': {
                            'targetability': {
                                'data': [0, 0, 1],
                                'percentages': [0, 0, 100],
                            },
                        },
                    },
                },
                'pathway_guide_data': {
                    DANGEROUS_TEXT: {
                        'guides': [{
                            'guide_sequence': DANGEROUS_TEXT,
                            'description': DANGEROUS_TEXT,
                        }],
                    },
                },
                'pathway_svgs': {},
            }
            pathways = [{
                'pathway_id': DANGEROUS_TEXT,
                'pathway_name': DANGEROUS_TEXT,
                'total_genes': 1,
                'total_gene_nodes': 1,
                'evaluable_gene_node_count': 1,
                'unique_xref_count': 1,
                'matched_genes': '[]',
                'genes_with_guides': 0,
                'pw_kfi_score': 0,
                'targetable_type': 'Not targetable',
                'coverage_rate': 0,
            }]

            final_path = root / 'empty_annotation_with_guides.tsv'
            final_path.write_text(
                '\t'.join([
                    'query_id', 'subject_id', 'identity',
                    'alignment_length', 'evalue', 'bit_score',
                    'description', 'gene_name', 'classification',
                    'unique_id', 'guide_targets',
                ]) + '\n',
                encoding='utf-8',
            )
            common_generator = ReportGenerator(config)
            common_results = self.standard_results(config, final_path)
            pathway_payload = generator._build_pathway_browser_data(
                pathways, summary
            )
            published = common_generator.publish_browser_payload(
                common_results,
                namespace='pathway',
                extra_globals={
                    '__TWHEEL_PATHWAY_PAYLOAD__': pathway_payload,
                },
            )
            document = generator.create_pathway_html_template(
                pathways,
                summary,
                {},
                common_results=common_results,
                common_report_generator=common_generator,
                payload_script_src=published['core_relative_path'],
            )
            payload_source = Path(published['core_path']).read_text(
                encoding='utf-8'
            )
            payload = self.extract_javascript_variable(
                payload_source, 'extraGlobals'
            )['__TWHEEL_PATHWAY_PAYLOAD__']
            common_payload = self.extract_javascript_variable(
                payload_source, 'payload'
            )['data']

        self.assertEqual(payload['pathways'][0]['pathway_id'], DANGEROUS_TEXT)
        self.assertEqual(
            payload['summary']['timestamp'], DANGEROUS_TEXT
        )
        self.assertEqual(
            common_payload['config']['input']['genome_assembly'],
            DANGEROUS_TEXT,
        )
        self.assertNotIn('<img src=x onerror="owned">', document)
        self.assertNotIn('pathway-data', document)
        self.assertNotIn('twheel-data', document)
        self.assertIn('\\u003c/script\\u003e', payload_source)
        self.assertNotIn('onclick="showPathwayGenes(\'', document)
        self.assertNotIn('onclick="openSgdLink(\'', document)
        self.assertNotIn('function generateSampleSVG', document)
        self.assertIn('validatedSVGMarkup(pathwaySvgData.svg_content)', document)
        self.assertIn("button.addEventListener(", document)

    def test_pathway_svg_escapes_every_external_node_field(self):
        label = '" onmouseover="owned & <tag>'
        database = 'DB" onload="owned'
        identifier = 'ID&<"'
        pathway_data = {
            'board_width': 400,
            'board_height': 300,
            'edges': [],
            'nodes': [
                {
                    'type': 'GeneProduct',
                    'label': label,
                    'database': database,
                    'xref_id': identifier,
                    'x': 50,
                    'y': 50,
                    'width': 80,
                    'height': 20,
                },
                {
                    'type': 'Metabolite',
                    'label': label,
                    'database': '',
                    'xref_id': '',
                    'x': 150,
                    'y': 50,
                    'width': 80,
                    'height': 20,
                },
                {
                    'type': 'Unknown',
                    'label': label,
                    'database': '',
                    'xref_id': '',
                    'x': 250,
                    'y': 50,
                    'width': 80,
                    'height': 20,
                },
            ],
        }
        with tempfile.TemporaryDirectory() as temp_dir_name:
            generator = PathwayReportGenerator(
                self.config(Path(temp_dir_name), 'safe')
            )
            svg = generator.generate_pathway_svg(pathway_data, {})

        root = ET.fromstring(svg)
        allowed_attributes = {
            'width', 'height', 'viewBox', 'd', 'class', 'cx', 'cy', 'rx',
            'ry', 'x', 'y', 'fill', 'stroke', 'stroke-width', 'opacity',
            'id', 'font-weight', 'font-size', 'data-gene-database',
            'data-gene-id', 'data-gene-key', 'data-sgd-id',
            'data-gene-name', 'title',
        }
        for element in root.iter():
            self.assertFalse(any(
                name.lower().startswith('on') for name in element.attrib
            ))
            self.assertTrue(set(element.attrib) <= allowed_attributes)
            for paint_name in ('fill', 'stroke'):
                if paint_name in element.attrib:
                    self.assertRegex(
                        element.attrib[paint_name],
                        re.compile(r'^(?:#[0-9a-f]{3,8}|none|white)$', re.I),
                    )
        ellipse = root.find('{http://www.w3.org/2000/svg}ellipse')
        rectangles = root.findall('{http://www.w3.org/2000/svg}rect')
        self.assertEqual(ellipse.attrib['data-gene-name'], label)
        self.assertEqual(ellipse.attrib['data-gene-database'], database)
        self.assertEqual(ellipse.attrib['data-gene-id'], identifier)
        self.assertEqual(rectangles[0].attrib['title'], label)
        self.assertEqual(rectangles[1].attrib['title'], label)


if __name__ == '__main__':
    unittest.main()
