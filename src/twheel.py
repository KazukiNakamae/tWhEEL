#!/usr/bin/env python3
"""
tWhEEL - Transcriptome Workflow for Hierarchical Exploration and Expression Levels

This script performs gene prediction, functional annotation, CRISPR guide design, 
and GO functional analysis from transcriptome assembly data to identify target genes 
for knockout experiments.
"""

import argparse
import sys
import os
import yaml
import logging
from pathlib import Path
from datetime import datetime

# Module imports
from modules.gene_prediction import GenePrediction
from modules.annotation import Annotation
from modules.sgd_annotation import Sgd_annotation
from modules.guide_design import GuideDesign
from modules.go_analysis import GOAnalysis
from modules.pathway_analysis import PathwayAnalysis
from modules.pathway_report_generator import PathwayReportGenerator
from modules.report_generator import ReportGenerator
from modules.input_mode import (
    EXISTING_PROTEIN_MODE,
    determine_protein_input_mode,
    get_configured_protein_path,
)


if hasattr(yaml, 'SafeLoader'):
    class _UniqueKeySafeLoader(yaml.SafeLoader):
        """Safe YAML loader that rejects duplicate mapping keys."""


    def _construct_unique_mapping(loader, node, deep=False):
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


    _UniqueKeySafeLoader.add_constructor(
        yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
        _construct_unique_mapping,
    )
else:
    _UniqueKeySafeLoader = None


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
        self.create_output_directory(analysis_name)
        self.setup_directories()
        self.setup_logging()
        
    def _load_config(self, config_path):
        """Load configuration file"""
        try:
            if _UniqueKeySafeLoader is None:
                raise RuntimeError('PyYAML SafeLoader is unavailable')
            with open(config_path, 'r', encoding='utf-8') as f:
                config = yaml.load(f, Loader=_UniqueKeySafeLoader)
            return config
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
            self.logger.info("Genome assembly provided - comprehensive off-target search will be performed")
        else:
            self.logger.warning("Genome assembly not provided. Off-target search will use transcriptome (less comprehensive)")
                
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


def main():
    """Main function"""
    parser = argparse.ArgumentParser(
        description='tWhEEL - Transcriptome Workflow for Hierarchical Exploration and Expression Levels',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python twheel.py --config config/default.yaml
  python twheel.py --config my_config.yaml --name my_analysis
  python twheel.py --config my_config.yaml --step annotation --name test_run
  python twheel.py --config my_config.yaml --step annotation --name tWhEEL_on_20250719_082938_E_okinawaensis_denovoanalysis
  python twheel.py --list-steps

Note:
  --name can specify either:
    - A simple analysis name (creates new timestamped directory)
    - A complete directory name starting with 'tWhEEL_on_' (reuses existing directory)
        """
    )
    
    parser.add_argument(
        '--config', '-c',
        default='config/default.yaml',
        help='Configuration file path (default: config/default.yaml)'
    )
    
    parser.add_argument(
        '--step', '-s',
        help='Execute specific step (omit to run all steps)'
    )
    
    parser.add_argument(
        '--name', '-n',
        help='Analysis name for output directory'
    )
    
    parser.add_argument(
        '--list-steps',
        action='store_true',
        help='List available steps'
    )
    
    parser.add_argument(
        '--version', '-v',
        action='version',
        version='tWhEEL version 1.0.0'
    )
    
    args = parser.parse_args()
    
    if args.list_steps:
        print("Available steps:")
        steps = ["gene_prediction", "annotation", "guide_design", "go_analysis", "report_generation"]
        for step in steps:
            print(f"  - {step}")
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
