import contextlib
import io
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest.mock import Mock, patch


try:
    import yaml  # noqa: F401
except ModuleNotFoundError:
    yaml_stub = types.ModuleType('yaml')
    yaml_stub.safe_load = lambda stream: {}
    sys.modules['yaml'] = yaml_stub

try:
    import requests  # noqa: F401
except ModuleNotFoundError:
    sys.modules['requests'] = types.ModuleType('requests')

SRC_DIR = Path(__file__).resolve().parents[1] / 'src'
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

import twheel


@unittest.skipUnless(
    hasattr(sys.modules['yaml'], 'SafeLoader'),
    'PyYAML is required for configuration loader tests',
)
class StrictYamlConfigTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir_context = tempfile.TemporaryDirectory()
        self.temp_dir = Path(self.temp_dir_context.name)
        self.pipeline = object.__new__(twheel.tWhEELPipeline)

    def tearDown(self):
        self.temp_dir_context.cleanup()

    def _assert_duplicate_is_rejected(self, yaml_text, duplicate_key):
        config_path = self.temp_dir / 'duplicate.yaml'
        config_path.write_text(yaml_text, encoding='utf-8')
        output = io.StringIO()

        with contextlib.redirect_stdout(output):
            with self.assertRaises(SystemExit) as error_context:
                self.pipeline._load_config(config_path)

        self.assertEqual(error_context.exception.code, 1)
        self.assertIn(
            f"found duplicate key ('{duplicate_key}')",
            output.getvalue(),
        )

    def test_duplicate_database_section_is_rejected(self):
        self._assert_duplicate_is_rejected(
            'database:\n'
            '  force_rebuild: true\n'
            'database:\n'
            '  force_rebuild: false\n',
            'database',
        )

    def test_duplicate_nested_force_rebuild_key_is_rejected(self):
        self._assert_duplicate_is_rejected(
            'database:\n'
            '  force_rebuild: true\n'
            '  force_rebuild: false\n',
            'force_rebuild',
        )

    def test_denovo_config_keeps_force_rebuild_enabled(self):
        config_path = (
            Path(__file__).resolve().parents[1]
            / 'config'
            / 'input_for_denovotranscriptomeassembly_GOAnalysis_20250718.yaml'
        )

        config = self.pipeline._load_config(config_path)

        self.assertIs(config['database']['force_rebuild'], True)

    def test_repository_configs_have_no_duplicate_keys(self):
        config_dir = Path(__file__).resolve().parents[1] / 'config'

        for config_path in sorted(config_dir.glob('*.yaml')):
            with self.subTest(config_path=config_path.name):
                config = self.pipeline._load_config(config_path)
                self.assertIsInstance(config, dict)


class OutputDirectoryConfigTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir_context = tempfile.TemporaryDirectory()
        self.temp_dir = Path(self.temp_dir_context.name)
        self.pipeline = object.__new__(twheel.tWhEELPipeline)

    def tearDown(self):
        self.temp_dir_context.cleanup()

    def test_configured_base_directory_is_used_as_output_root(self):
        output_root = self.temp_dir / 'custom-output-root'
        self.pipeline.config = {
            'output': {
                'base_dir': str(output_root),
                'intermediate': 'ignored-legacy-intermediate',
            }
        }

        self.pipeline.create_output_directory(
            'tWhEEL_on_20260716_custom_analysis'
        )
        self.pipeline.setup_directories()

        expected_output = (
            output_root / 'tWhEEL_on_20260716_custom_analysis'
        )
        self.assertEqual(
            Path(self.pipeline.config['output']['base_dir']),
            expected_output,
        )
        self.assertEqual(
            Path(self.pipeline.config['output']['intermediate']),
            expected_output / 'intermediate',
        )
        self.assertTrue((expected_output / 'intermediate').is_dir())
        self.assertTrue((expected_output / 'web_data').is_dir())

    def test_complete_output_path_is_not_nested_twice(self):
        run_name = 'tWhEEL_on_20260716_existing_analysis'
        configured_run_path = self.temp_dir / run_name
        self.pipeline.config = {
            'output': {'base_dir': str(configured_run_path)}
        }

        self.pipeline.create_output_directory(run_name)

        self.assertEqual(
            Path(self.pipeline.config['output']['base_dir']),
            configured_run_path,
        )

    def test_empty_base_directory_is_rejected(self):
        self.pipeline.config = {'output': {'base_dir': '  '}}

        with self.assertRaisesRegex(
            ValueError, 'output.base_dir must be a non-empty path'
        ):
            self.pipeline.create_output_directory('analysis')


class OptionalGenomeAssemblyConfigTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir_context = tempfile.TemporaryDirectory()
        self.temp_dir = Path(self.temp_dir_context.name)
        self.transcriptome_path = self.temp_dir / 'transcriptome.fa'
        self.transcriptome_path.write_text(
            '>tx1\nACGT\n', encoding='utf-8'
        )
        self.pipeline = object.__new__(twheel.tWhEELPipeline)
        self.pipeline.config = {
            'input': {'transcriptome': str(self.transcriptome_path)},
            'execution': {'steps': []},
        }
        self.pipeline.logger = Mock()

    def tearDown(self):
        self.temp_dir_context.cleanup()

    def test_validate_inputs_accepts_omitted_genome_assembly(self):
        self.assertTrue(self.pipeline.validate_inputs())
        self.pipeline.logger.warning.assert_called_once_with(
            "Genome assembly not provided. Off-target search will use "
            "transcriptome (less comprehensive)"
        )


class PathwayAnalysisHandoffTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir_context = tempfile.TemporaryDirectory()
        self.temp_dir = Path(self.temp_dir_context.name)
        self.pipeline = object.__new__(twheel.tWhEELPipeline)
        self.pipeline.config = {
            'output': {
                'base_dir': str(self.temp_dir),
                'intermediate': str(self.temp_dir),
            }
        }
        self.pipeline.logger = Mock()

    def tearDown(self):
        self.temp_dir_context.cleanup()

    def artifact_paths(self):
        filenames = {
            'kfi_table_path': 'pathway_kfi_table.csv',
            'fully_targetable_path': 'pathway_kfi_fully_targetable.csv',
            'mapping_path': 'pathway_gene_mapping.tsv',
            'gene_match_audit_path': 'pathway_gene_match_audit.tsv',
            'annotation_identifier_audit_path': (
                'pathway_annotation_identifier_audit.tsv'
            ),
        }
        paths = {}
        for artifact_name, filename in filenames.items():
            path = self.temp_dir / filename
            path.write_text('header\n', encoding='utf-8')
            paths[artifact_name] = str(path)
        return paths

    def test_none_result_is_not_reported_as_success(self):
        with patch.object(twheel, 'PathwayAnalysis') as analysis_class:
            analysis_class.return_value.run.return_value = None

            success = self.pipeline.run_step('pathway_analysis')

        self.assertFalse(success)

    def test_zero_pathways_with_complete_empty_artifacts_is_successful(self):
        result = {
            'processed_pathways': [],
            'statistics': {'processable_pathways': 0},
            'output_files': self.artifact_paths(),
        }
        with patch.object(twheel, 'PathwayAnalysis') as analysis_class:
            analysis_class.return_value.run.return_value = result

            success = self.pipeline.run_step('pathway_analysis')

        self.assertTrue(success)
        self.pipeline.logger.warning.assert_called_once_with(
            "Pathway analysis completed with 0 reportable pathways; "
            "empty report artifacts were created"
        )

    def test_missing_empty_artifact_stops_pipeline_before_report(self):
        artifact_paths = self.artifact_paths()
        Path(artifact_paths['mapping_path']).unlink()
        result = {
            'processed_pathways': [],
            'statistics': {'processable_pathways': 0},
            'output_files': artifact_paths,
        }
        with patch.object(twheel, 'PathwayAnalysis') as analysis_class:
            analysis_class.return_value.run.return_value = result

            success = self.pipeline.run_step('pathway_analysis')

        self.assertFalse(success)


if __name__ == '__main__':
    unittest.main()
