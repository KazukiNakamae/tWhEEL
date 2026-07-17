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
import twheel_metadata


class CliMetadataTests(unittest.TestCase):
    def test_version_output_uses_canonical_metadata(self):
        output = io.StringIO()

        with contextlib.redirect_stdout(output):
            with self.assertRaises(SystemExit) as error_context:
                twheel.build_parser().parse_args(['--version'])

        self.assertEqual(error_context.exception.code, 0)
        self.assertEqual(twheel.__version__, twheel_metadata.__version__)
        self.assertEqual(
            output.getvalue().strip(),
            f'tWhEEL version {twheel_metadata.__version__}',
        )

    def test_help_and_step_choices_cover_all_implemented_steps(self):
        expected_steps = (
            'gene_prediction',
            'annotation',
            'sgd_annotation',
            'guide_design',
            'go_analysis',
            'pathway_analysis',
            'report_generation',
            'pathway_report_generation',
        )
        parser = twheel.build_parser()
        step_action = next(
            action for action in parser._actions if action.dest == 'step'
        )

        self.assertEqual(twheel_metadata.STEP_NAMES, expected_steps)
        self.assertEqual(tuple(step_action.choices), expected_steps)
        help_text = parser.format_help()
        for step_name in expected_steps:
            self.assertIn(step_name, help_text)
        self.assertIn(
            'config/pw_kfi_existing_protein_genome.yaml', help_text
        )
        self.assertIn('does not run or validate prerequisites', help_text)
        self.assertIn('preliminary knockout-feasibility', help_text)
        self.assertIn('maximum mismatch is fixed at 2', help_text)
        self.assertIn(
            'crispr.max_mismatch is not a YAML option', help_text
        )
        self.assertIn('MM0=0', help_text)
        self.assertIn('intended gene', help_text)
        self.assertIn('GFF/GTF exon structure', help_text)

    def test_list_steps_uses_shared_descriptions_without_loading_config(self):
        output = io.StringIO()

        with contextlib.redirect_stdout(output):
            return_code = twheel.main(['--list-steps'])

        self.assertEqual(return_code, 0)
        listing = output.getvalue()
        self.assertIn(f'tWhEEL {twheel_metadata.__version__}', listing)
        for step in twheel_metadata.STEP_DEFINITIONS:
            self.assertIn(step.name, listing)
            self.assertIn(step.workflow, listing)
            self.assertIn(step.summary, listing)
        self.assertIn('does not execute or validate prerequisites', listing)

    def test_help_exits_without_initializing_pipeline(self):
        output = io.StringIO()

        with patch.object(twheel, 'tWhEELPipeline') as pipeline_class:
            with contextlib.redirect_stdout(output):
                with self.assertRaises(SystemExit) as error_context:
                    twheel.main(['--help'])

        self.assertEqual(error_context.exception.code, 0)
        pipeline_class.assert_not_called()
        self.assertIn('Implemented steps:', output.getvalue())

    def test_unknown_step_is_rejected_by_argument_parser(self):
        error_output = io.StringIO()

        with contextlib.redirect_stderr(error_output):
            with self.assertRaises(SystemExit) as error_context:
                twheel.build_parser().parse_args(['--step', 'unknown_step'])

        self.assertEqual(error_context.exception.code, 2)
        self.assertIn('invalid choice', error_output.getvalue())


@unittest.skipUnless(
    hasattr(sys.modules['yaml'], 'SafeLoader'),
    'PyYAML is required for configuration loader tests',
)
class StrictYamlConfigTests(unittest.TestCase):
    EXPECTED_CONFIG_FILES = {
        'default.yaml',
        'go_kfi_existing_protein.yaml',
        'pw_kfi_existing_protein.yaml',
        'go_kfi_existing_protein_genome.yaml',
        'pw_kfi_existing_protein_genome.yaml',
    }
    GO_STEPS = [
        'gene_prediction',
        'annotation',
        'guide_design',
        'go_analysis',
        'report_generation',
    ]
    PATHWAY_STEPS = [
        'gene_prediction',
        'sgd_annotation',
        'guide_design',
        'pathway_analysis',
        'pathway_report_generation',
    ]

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

    def test_repository_configs_have_no_duplicate_keys(self):
        config_dir = Path(__file__).resolve().parents[1] / 'config'
        config_paths = sorted(config_dir.glob('*.yaml'))

        self.assertEqual(
            {path.name for path in config_paths},
            self.EXPECTED_CONFIG_FILES,
        )

        for config_path in config_paths:
            with self.subTest(config_path=config_path.name):
                config = self.pipeline._load_config(config_path)
                self.assertIsInstance(config, dict)

    def test_removed_max_mismatch_option_is_rejected(self):
        config_path = self.temp_dir / 'removed-option.yaml'
        config_path.write_text(
            'crispr:\n  max_mismatch: 3\n', encoding='utf-8'
        )
        output = io.StringIO()

        with contextlib.redirect_stdout(output):
            with self.assertRaises(SystemExit) as error_context:
                self.pipeline._load_config(config_path)

        self.assertEqual(error_context.exception.code, 1)
        self.assertIn(
            'crispr.max_mismatch has been removed', output.getvalue()
        )
        self.assertIn('fixed at 2', output.getvalue())

    def test_repository_configs_match_the_five_supported_profiles(self):
        config_dir = Path(__file__).resolve().parents[1] / 'config'
        scenarios = {
            'default.yaml': (False, False, self.GO_STEPS),
            'go_kfi_existing_protein.yaml': (
                True, False, self.GO_STEPS,
            ),
            'pw_kfi_existing_protein.yaml': (
                True, False, self.PATHWAY_STEPS,
            ),
            'go_kfi_existing_protein_genome.yaml': (
                True, True, self.GO_STEPS,
            ),
            'pw_kfi_existing_protein_genome.yaml': (
                True, True, self.PATHWAY_STEPS,
            ),
        }

        for filename, (has_protein, has_genome, steps) in scenarios.items():
            with self.subTest(config=filename):
                config = self.pipeline._load_config(config_dir / filename)
                inputs = config['input']

                self.assertEqual(
                    bool(inputs.get('protein_sequences')), has_protein
                )
                self.assertEqual(
                    bool(inputs.get('genome_assembly')), has_genome
                )
                self.assertEqual(config['execution']['steps'], steps)
                self.assertNotIn('intermediate', config['output'])
                self.assertNotIn('docker', config)
                self.assertNotIn('parallel', config['execution'])
                self.assertNotIn('LOCAL_DATA', repr(config))
                self.assertNotIn('max_mismatch', config['crispr'])

                if has_protein:
                    self.assertTrue(inputs['protein_transcript_mapping'])
                    self.assertIs(
                        config['transdecoder']['skip_if_protein_provided'],
                        True,
                    )
                    self.assertNotIn('docker_image', config['transdecoder'])
                else:
                    self.assertNotIn('protein_transcript_mapping', inputs)
                    self.assertTrue(config['transdecoder']['docker_image'])

                if steps == self.GO_STEPS:
                    self.assertTrue(inputs['protein_db_for_goanalysis'])
                    self.assertIn('go_analysis', config)
                    self.assertIn('report', config)
                    self.assertNotIn('pathway_analysis', config)
                else:
                    self.assertTrue(inputs['sgd_db'])
                    self.assertIn('pathway_analysis', config)
                    self.assertNotIn('go_analysis', config)
                    self.assertNotIn('report', config)


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

    def test_validate_inputs_accepts_preliminary_no_genome_profile(self):
        self.assertTrue(self.pipeline.validate_inputs())
        self.pipeline.logger.warning.assert_called_once_with(
            "Genome assembly not provided. Off-target search will use the "
            "transcriptome; this mode is for preliminary knockout-"
            "feasibility screening only"
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
