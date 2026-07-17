#!/usr/bin/env python3
"""
tWhEEL - Transcriptome Workflow for Hierarchical Exploration and Expression Levels

This script runs the GO and SGD/WikiPathways analysis profiles for protein
prediction, functional annotation, CRISPR guide design, feasibility scoring,
and static report generation.
"""

import argparse
import sys
import os
import logging
from pathlib import Path
from datetime import datetime

from modules.input_mode import (
    EXISTING_PROTEIN_MODE,
    determine_protein_input_mode,
    get_configured_protein_path,
)
from modules.guide_design_contract import (
    CRISFLASH_MAX_MISMATCH,
    SITE_SPECIFIC_FILTER_LABEL,
)
from modules.multi_gene_alignment import annotation_workspace_lock
from twheel_metadata import (
    STEP_NAMES,
    __version__,
    format_step_listing,
)


# Runtime-only imports are loaded after argument parsing so metadata commands
# remain available before the analysis environment is installed.
yaml = None
_UniqueKeySafeLoader = None
GenePrediction = None
Annotation = None
Sgd_annotation = None
GuideDesign = None
GOAnalysis = None
PathwayAnalysis = None
PathwayReportGenerator = None
ReportGenerator = None


def _ensure_yaml_loader():
    """Return a safe YAML loader that rejects duplicate mapping keys."""
    global yaml, _UniqueKeySafeLoader

    if yaml is None:
        try:
            import yaml as yaml_module
        except ModuleNotFoundError as error:
            raise RuntimeError(
                'PyYAML is required to load a pipeline configuration. '
                'Install the dependencies from requirements.txt.'
            ) from error
        yaml = yaml_module

    if _UniqueKeySafeLoader is not None:
        return _UniqueKeySafeLoader
    if not hasattr(yaml, 'SafeLoader'):
        raise RuntimeError('PyYAML SafeLoader is unavailable')

    class UniqueKeySafeLoader(yaml.SafeLoader):
        """Safe YAML loader that rejects duplicate mapping keys."""

    def construct_unique_mapping(loader, node, deep=False):
        """Construct a mapping without silently replacing earlier keys."""
        mapping = {}
        for key_node, value_node in node.value:
            key = loader.construct_object(key_node, deep=deep)
            try:
                hash(key)
            except TypeError as error:
                raise yaml.constructor.ConstructorError(
                    'while constructing a mapping',
                    node.start_mark,
                    'found an unhashable key',
                    key_node.start_mark,
                ) from error
            if key in mapping:
                raise yaml.constructor.ConstructorError(
                    'while constructing a mapping',
                    node.start_mark,
                    f'found duplicate key ({key!r})',
                    key_node.start_mark,
                )
            mapping[key] = loader.construct_object(value_node, deep=deep)
        return mapping

    UniqueKeySafeLoader.add_constructor(
        yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
        construct_unique_mapping,
    )
    _UniqueKeySafeLoader = UniqueKeySafeLoader
    return _UniqueKeySafeLoader


def _load_runtime_modules():
    """Load analysis classes only when a pipeline run is requested."""
    global GenePrediction, Annotation, Sgd_annotation, GuideDesign
    global GOAnalysis, PathwayAnalysis, PathwayReportGenerator, ReportGenerator

    if GenePrediction is None:
        from modules.gene_prediction import GenePrediction as implementation
        GenePrediction = implementation
    if Annotation is None:
        from modules.annotation import Annotation as implementation
        Annotation = implementation
    if Sgd_annotation is None:
        from modules.sgd_annotation import Sgd_annotation as implementation
        Sgd_annotation = implementation
    if GuideDesign is None:
        from modules.guide_design import GuideDesign as implementation
        GuideDesign = implementation
    if GOAnalysis is None:
        from modules.go_analysis import GOAnalysis as implementation
        GOAnalysis = implementation
    if PathwayAnalysis is None:
        from modules.pathway_analysis import PathwayAnalysis as implementation
        PathwayAnalysis = implementation
    if PathwayReportGenerator is None:
        from modules.pathway_report_generator import (
            PathwayReportGenerator as implementation,
        )
        PathwayReportGenerator = implementation
    if ReportGenerator is None:
        from modules.report_generator import ReportGenerator as implementation
        ReportGenerator = implementation


def _validate_config_contract(config):
    """Reject removed options instead of silently ignoring user input."""
    if not isinstance(config, dict):
        raise ValueError("The YAML document must contain a mapping")

    crispr = config.get('crispr')
    if isinstance(crispr, dict) and 'max_mismatch' in crispr:
        raise ValueError(
            'crispr.max_mismatch has been removed. CrisFlash maximum '
            f'mismatch is fixed at {CRISFLASH_MAX_MISMATCH}; remove this '
            'key from the YAML configuration.'
        )

    return config


class tWhEELPipeline:
    """Main class for tWhEEL pipeline"""
    
    def __init__(self, config_path, analysis_name=None):
        """
        Initialize the pipeline
        
        Args:
            config_path (str): Path to configuration file
            analysis_name (str, optional): Analysis name for output directory
        """
        self.config = self._load_config(config_path)
        _load_runtime_modules()
        self.create_output_directory(analysis_name)
        self.setup_directories()
        self.setup_logging()
        
    def _load_config(self, config_path):
        """Load configuration file"""
        try:
            loader = _ensure_yaml_loader()
            with open(config_path, 'r', encoding='utf-8') as f:
                config = yaml.load(f, Loader=loader)
            return _validate_config_contract(config)
        except Exception as e:
            print(f"Failed to load configuration file: {e}")
            sys.exit(1)
    
    def create_output_directory(self, analysis_name=None):
        """Create output directory with timestamp and analysis name"""
        output_config = self.config.get('output')
        if not isinstance(output_config, dict):
            raise ValueError("The output configuration must be a mapping")
        configured_base_dir = output_config.get('base_dir')
        if (
            not isinstance(configured_base_dir, str)
            or not configured_base_dir.strip()
        ):
            raise ValueError("output.base_dir must be a non-empty path")
        output_root = os.path.expanduser(configured_base_dir.strip())

        if analysis_name:
            # Check if a complete directory name is provided
            if analysis_name.startswith('tWhEEL_on_'):
                # Use the complete directory name as provided
                dir_name = analysis_name
                print(f"Using existing directory name: {analysis_name}")
            else:
                # Generate new directory name with timestamp
                timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
                dir_name = f"tWhEEL_on_{timestamp}_{analysis_name}"
                print(f"Creating new directory with timestamp")
        else:
            # No analysis name provided, use timestamp only
            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
            dir_name = f"tWhEEL_on_{timestamp}"
            print(f"Creating new directory with timestamp only")
        
        # Reusing a complete run name must not duplicate the same path segment.
        normalized_output_root = os.path.normpath(output_root)
        if os.path.basename(normalized_output_root) == dir_name:
            output_dir = normalized_output_root
        else:
            output_dir = os.path.join(output_root, dir_name)
        output_config['base_dir'] = output_dir
        output_config['intermediate'] = os.path.join(
            output_dir, 'intermediate'
        )
        
        print(f"Output directory: {output_dir}")
    
    def setup_logging(self):
        """Setup logging configuration"""
        log_level = logging.DEBUG if self.config['execution']['verbose'] else logging.INFO
        
        # Log format
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )
        
        # Console handler
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(formatter)
        
        # File handler
        log_file = os.path.join(self.config['output']['base_dir'], 'twheel.log')
        file_handler = logging.FileHandler(log_file, encoding='utf-8')
        file_handler.setFormatter(formatter)
        
        # Root logger setup
        logging.basicConfig(
            level=log_level,
            handlers=[console_handler, file_handler]
        )
        
        self.logger = logging.getLogger(__name__)
        self.logger.info("Starting tWhEEL pipeline")
    
    def setup_directories(self):
        """Create output directories"""
        directories = [
            self.config['output']['base_dir'],
            self.config['output']['intermediate'],
            os.path.join(self.config['output']['base_dir'], 'web_data')
        ]
        
        for directory in directories:
            Path(directory).mkdir(parents=True, exist_ok=True)
            print(f"Created directory: {directory}")
    
    def validate_inputs(self):
        """Validate input files"""
        # Transcriptome is always required (tWhEEL is a transcriptome-based pipeline)
        transcriptome = self.config['input'].get('transcriptome', '')
        if not transcriptome:
            self.logger.error("Required parameter not set: transcriptome")
            self.logger.error("tWhEEL is a transcriptome-based pipeline - transcriptome file is always required")
            return False
        if not os.path.exists(transcriptome):
            self.logger.error(f"Transcriptome file not found: {transcriptome}")
            return False
        
        # Check if any steps requiring protein_db_for_goanalysis are planned to run
        steps = self.config['execution'].get('steps', [])
        requires_protein_db = any(step in steps for step in ['annotation', 'go_analysis', 'report_generation'])
        
        # Protein database for GO analysis is conditionally required
        protein_db_for_goanalysis = self.config['input'].get('protein_db_for_goanalysis', '')
        if requires_protein_db:
            if not protein_db_for_goanalysis:
                self.logger.error("Required parameter not set: protein_db_for_goanalysis")
                self.logger.error("protein_db_for_goanalysis is required for steps: annotation, go_analysis, report_generation")
                return False
            if not os.path.exists(protein_db_for_goanalysis):
                self.logger.error(f"Protein database file not found: {protein_db_for_goanalysis}")
                return False
        else:
            if protein_db_for_goanalysis:
                self.logger.info(f"protein_db_for_goanalysis specified but not used in current pipeline steps")
            else:
                self.logger.info("protein_db_for_goanalysis not required for current pipeline steps (pathway analysis only)")
        
        # sgd_annotation step: validate sgd_db (S. cerevisiae-specific)
        # Each organism-specific annotation step validates its own database here.
        # Future organism support (e.g., wormbase_annotation) should add a similar block below.
        if 'sgd_annotation' in steps:
            sgd_db = self.config['input'].get('sgd_db', '')
            if not sgd_db:
                self.logger.error("Required parameter not set: sgd_db")
                self.logger.error("sgd_db is required for sgd_annotation step (S. cerevisiae pathway analysis)")
                return False
            if not os.path.exists(sgd_db):
                self.logger.error(f"SGD database file not found: {sgd_db}")
                return False
            self.logger.info(f"SGD database validated: {sgd_db}")
        
        # pathway_analysis / pathway_report_generation steps: validate WikiPathways GPML directory
        pathway_steps_present = any(step in steps for step in ['pathway_analysis', 'pathway_report_generation'])
        if pathway_steps_present:
            pathway_config = self.config.get('pathway_analysis', {})
            gpml_dir = pathway_config.get('wikipathways_gpml_dir', '')
            if not gpml_dir:
                self.logger.error("Required parameter not set: pathway_analysis.wikipathways_gpml_dir")
                self.logger.error("wikipathways_gpml_dir is required for pathway_analysis and pathway_report_generation steps")
                return False
            if not os.path.exists(gpml_dir):
                self.logger.error(f"WikiPathways GPML directory not found: {gpml_dir}")
                return False
            self.logger.info(f"WikiPathways GPML directory validated: {gpml_dir}")
        
        # Check processing mode: existing protein sequences (optional) + transcriptome (required)
        protein_sequences = get_configured_protein_path(self.config)
        input_mode = determine_protein_input_mode(self.config)

        if input_mode == EXISTING_PROTEIN_MODE:
            # Existing protein sequences mode - still requires transcriptome for downstream analysis
            self.logger.info("Existing protein sequences mode: Using pre-computed proteins + transcriptome data")
            self.logger.info(f"Protein sequences file: {protein_sequences}")
            self.logger.info(f"Transcriptome file: {transcriptome}")
            self.logger.info("TransDecoder will be skipped, but transcriptome data will be used for guide design and other analyses")
            
        else:
            # Standard TransDecoder mode
            self.logger.info("TransDecoder mode: Using transcriptome for ORF prediction")
            self.logger.info(f"Transcriptome file: {transcriptome}")
            
            if protein_sequences and not os.path.exists(protein_sequences):
                self.logger.warning(f"Specified protein sequences file not found: {protein_sequences}")
                self.logger.warning("Falling back to TransDecoder mode")
            elif protein_sequences:
                self.logger.info(
                    "Existing protein sequences were provided, but "
                    "skip_if_protein_provided is disabled"
                )
        
        # Optional genome assembly validation
        genome_assembly = self.config['input'].get('genome_assembly', '')
        if genome_assembly:
            if not os.path.exists(genome_assembly):
                self.logger.error(f"Genome assembly file not found: {genome_assembly}")
                return False
            self.logger.info(
                "Genome assembly provided for off-target searching; guide "
                "candidate sequences remain transcriptome-derived"
            )
        else:
            self.logger.warning(
                "Genome assembly not provided. Off-target search will use the "
                "transcriptome; this mode is for preliminary knockout-"
                "feasibility screening only"
            )
                
        self.logger.info("Input file validation completed")
        return True

    def _validate_pathway_analysis_result(self, result):
        """Validate the pathway-analysis handoff before report generation."""
        if not isinstance(result, dict):
            self.logger.error(
                "Pathway analysis did not return a structured result"
            )
            return False

        processed_pathways = result.get('processed_pathways')
        if not isinstance(processed_pathways, list):
            self.logger.error(
                "Pathway analysis result is missing processed_pathways"
            )
            return False

        output_files = result.get('output_files')
        if not isinstance(output_files, dict):
            self.logger.error(
                "Pathway analysis result is missing output_files"
            )
            return False

        required_artifacts = {
            'kfi_table_path',
            'fully_targetable_path',
            'mapping_path',
            'gene_match_audit_path',
            'annotation_identifier_audit_path',
        }
        invalid_artifacts = []
        for artifact_name in sorted(required_artifacts):
            artifact_path = output_files.get(artifact_name, '')
            if (
                not artifact_path
                or not os.path.isfile(artifact_path)
                or os.path.getsize(artifact_path) == 0
            ):
                invalid_artifacts.append(artifact_name)

        if invalid_artifacts:
            self.logger.error(
                "Pathway analysis did not produce required artifacts: "
                + ', '.join(invalid_artifacts)
            )
            return False
        return True
    
    def run_step(self, step_name):
        """Execute one step while holding the shared workspace generation lock."""
        with annotation_workspace_lock(self.config['output']['intermediate']):
            return self._run_step_locked(step_name)

    def _run_step_locked(self, step_name):
        """Execute individual step"""
        self.logger.info(f"Starting step: {step_name}")
        start_time = datetime.now()
        
        try:
            if step_name == "gene_prediction":
                predictor = GenePrediction(self.config)
                predictor.run()
                
            elif step_name == "annotation":
                annotator = Annotation(self.config)
                annotator.run()

            elif step_name == "sgd_annotation":
                sgd_annotator = Sgd_annotation(self.config)
                sgd_annotator.run()
                
            elif step_name == "guide_design":
                guide_designer = GuideDesign(self.config)
                guide_designer.run()
                
            elif step_name == "go_analysis":
                go_analyzer = GOAnalysis(self.config)
                go_analyzer.run()
                
            elif step_name == "pathway_analysis":
                pathway_analyzer = PathwayAnalysis(self.config)
                pathway_result = pathway_analyzer.run()
                if not self._validate_pathway_analysis_result(pathway_result):
                    return False
                if not pathway_result['processed_pathways']:
                    self.logger.warning(
                        "Pathway analysis completed with 0 reportable pathways; "
                        "empty report artifacts were created"
                    )
                
            elif step_name == "pathway_report_generation":
                pathway_report_gen = PathwayReportGenerator(self.config)
                pathway_report_gen.run()
                
            elif step_name == "report_generation":
                report_gen = ReportGenerator(self.config)
                report_gen.run()
                
            else:
                self.logger.error(f"Unknown step name: {step_name}")
                return False
                
        except Exception as e:
            self.logger.error(f"Error occurred in step {step_name}: {e}")
            return False
        
        end_time = datetime.now()
        duration = end_time - start_time
        self.logger.info(f"Step {step_name} completed (execution time: {duration})")
        return True
    
    def run(self):
        """Execute the complete pipeline as one workspace generation."""
        with annotation_workspace_lock(self.config['output']['intermediate']):
            return self._run_locked()

    def _run_locked(self):
        """Execute entire pipeline"""
        if not self.validate_inputs():
            return False
            
        start_time = datetime.now()
        self.logger.info("Starting tWhEEL pipeline execution")
        
        steps = self.config['execution']['steps']
        
        for step in steps:
            if not self.run_step(step):
                self.logger.error(f"Pipeline stopped at step {step}")
                return False
        
        end_time = datetime.now()
        total_duration = end_time - start_time
        self.logger.info(f"tWhEEL pipeline completed successfully (total execution time: {total_duration})")
        
        return True


def build_parser():
    """Build the public command-line parser from canonical CLI metadata."""
    step_listing = format_step_listing()
    parser = argparse.ArgumentParser(
        description=(
            'tWhEEL - Transcriptome Workflow for Hierarchical Exploration '
            'and Expression Levels'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"""
Examples:
  # Run every step after editing the GO-profile template.
  python src/twheel.py --config config/default.yaml --name my_analysis

  # Run every step after replacing paths in the recommended Pathway profile.
  python src/twheel.py \\
    --config config/pw_kfi_existing_protein_genome.yaml \\
    --name my_pathway

  # Rerun exactly one step in an existing run directory.
  python src/twheel.py --config config/default.yaml \\
    --step go_analysis \\
    --name tWhEEL_on_YYYYMMDD_HHMMSS_my_analysis

  python src/twheel.py --list-steps
  python src/twheel.py --version

Implemented steps:
{step_listing}

Execution behavior:
  Without --step, entries in execution.steps run sequentially in config order.
  --step runs only the named step; it does not run or validate prerequisites.
  A simple --name creates a timestamped run directory. A complete name starting
  with tWhEEL_on_ reuses that run directory but does not detect completed steps.

Guide-design scope:
  A no-genome profile is only for preliminary knockout-feasibility screening.
  CrisFlash maximum mismatch is fixed at {CRISFLASH_MAX_MISMATCH};
  crispr.max_mismatch is not a YAML option.
  With a supplied genome, the normal {SITE_SPECIFIC_FILTER_LABEL} filter excludes MM0=0
  candidates, including junction-spanning candidates absent from the genome.
  tWhEEL does not confirm that the retained perfect match is the intended gene
  locus or validate it against GFF/GTF exon structure.
        """
    )
    
    parser.add_argument(
        '--config', '-c',
        default='config/default.yaml',
        help=(
            'YAML configuration file for the selected GO or Pathway profile '
            '(default: config/default.yaml)'
        ),
    )
    
    parser.add_argument(
        '--step', '-s',
        choices=STEP_NAMES,
        metavar='STEP',
        help=(
            'run exactly one implemented step without running prerequisites; '
            'omit to run execution.steps from the configuration'
        ),
    )
    
    parser.add_argument(
        '--name', '-n',
        help=(
            'simple run label, or a complete tWhEEL_on_... directory name to '
            'reuse an existing run'
        ),
    )
    
    parser.add_argument(
        '--list-steps',
        action='store_true',
        help='list implemented steps with workflow and purpose, then exit',
    )
    
    parser.add_argument(
        '--version', '-v',
        action='version',
        version=f'tWhEEL version {__version__}',
        help='show the tWhEEL release version and exit',
    )

    return parser


def main(argv=None):
    """Run the public CLI."""
    parser = build_parser()
    args = parser.parse_args(argv)
    
    if args.list_steps:
        print(f"tWhEEL {__version__} implemented steps:")
        print(format_step_listing())
        print(
            "\nA single --step run does not execute or validate prerequisites."
        )
        return 0
    
    # Execute pipeline
    try:
        pipeline = tWhEELPipeline(args.config, args.name)
        
        if args.step:
            # Execute specific step only
            success = pipeline.run_step(args.step)
        else:
            # Execute all steps
            success = pipeline.run()
            
        return 0 if success else 1
        
    except KeyboardInterrupt:
        print("\nExecution interrupted by user")
        return 1
    except Exception as e:
        print(f"Unexpected error occurred: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
