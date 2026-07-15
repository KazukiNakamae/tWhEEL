# tWhEEL Software Specification

## Table of Contents
1. [Overview](#overview)
2. [Architecture](#architecture)
3. [Module Specifications](#module-specifications)
4. [Data Flow](#data-flow)
5. [API Documentation](#api-documentation)
6. [Development Guidelines](#development-guidelines)
7. [Extension Guide](#extension-guide)
8. [Code Quality Standards](#code-quality-standards)

---

## Overview

**tWhEEL** (Transcriptome Workflow for Hierarchical Exploration and Expression Levels) is a comprehensive bioinformatics pipeline designed for:
- Transcriptome analysis and ORF prediction
- Functional annotation of genes
- CRISPR guide RNA design with specificity analysis
- Gene Ontology-based knockout feasibility assessment
- Interactive reporting and data visualization

### Design Principles
- **Modularity**: Each analysis step is a self-contained module
- **Transcriptome-centric**: Transcriptome data is the foundation for all analyses
- **Flexibility**: Support for both TransDecoder mode and existing protein sequences
- **Scalability**: Designed to handle large transcriptomes (100k+ sequences)
- **Traceability**: Comprehensive logging and intermediate file preservation

---

## Architecture

### High-Level System Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                    User Configuration                        │
│              (config/*.yaml files)                           │
└────────────────────────┬────────────────────────────────────┘
                         │
┌────────────────────────▼────────────────────────────────────┐
│                tWhEELPipeline (Main Orchestrator)          │
│              (src/twheel.py)                                │
└────────────────────────┬────────────────────────────────────┘
                         │
        ┌────────────────┼────────────────┐
        │                │                │
┌───────▼───────┐ ┌─────▼──────┐ ┌──────▼───────┐
│ Gene Pred.    │ │ Annotation │ │ Guide Design │
│ (TransDecoder)│ │ (FASTA)    │ │ (CRISPRFlash)│
└────────┬──────┘ └──────┬─────┘ └───────┬──────┘
         │                │              │
         └────────────────┼──────────────┘
                          │
         ┌────────────────┼────────────────┐
         │                │                │
    ┌────▼──────┐  ┌──────▼────┐  ┌───────▼────┐
    │GO Analysis│  │Pathway    │  │Report      │
    │(GO-KFI)   │  │Analysis   │  │Generation  │
    └────┬──────┘  └──────┬────┘  └───────┬────┘
         │                │               │
         └────────────────┼───────────────┘
                          │
         ┌────────────────▼────────────────┐
         │    Output Generation           │
         │  (HTML, CSV, JSON, SQLite)     │
         └──────────────────────────────────┘
```

### Module Dependency Graph

```
twheel.py (main entry point)
    ├── GenePrediction
    │   └── docker (TransDecoder)
    ├── Annotation
    │   ├── docker (FASTA3)
    │   └── Sgd_annotation
    ├── GuideDesign
    │   └── docker (CRISPRFlash)
    ├── GOAnalysis
    │   ├── docker (add-go-columns)
    │   ├── docker (compute_go_kfi)
    │   └── docker (count_go_parent_terms)
    ├── PathwayAnalysis
    │   └── sqlite3
    ├── PathwayReportGenerator
    │   └── sqlite3
    └── ReportGenerator
        ├── sqlite3
        └── requests (QuickGO API)
```

---

## Module Specifications

### 1. GenePrediction Module
**File:** `src/modules/gene_prediction.py`

#### Purpose
ORF prediction from transcriptome sequences using TransDecoder and quality assessment using BUSCO.

#### Key Classes
- `GenePrediction`: Main class for gene prediction workflow

#### Configuration Parameters
```yaml
transdecoder:
  docker_image: "quay.io/biocontainers/transdecoder:5.7.1--pl5321hdfd78af_0"
  min_protein_len: 85
  skip_if_protein_provided: true

busco:
  docker_image: "ezlabgva/busco:v5.2.2_cv1"
  lineage: "metazoa_odb10"
  cpu: 8
```

#### Input/Output
- **Input:** Transcriptome FASTA file
- **Output:** 
  - `query.pep` - Predicted protein sequences
  - `transdecoder_output/` - TransDecoder results
  - `busco_output/` - BUSCO assessment

#### Key Methods
```python
def __init__(self, config):
    """Initialize with configuration"""

def validate_inputs(self):
    """Validate input transcriptome file"""

def setup_docker(self):
    """Configure Docker containers"""

def run_transdecoder(self):
    """Execute TransDecoder prediction"""

def run_busco(self):
    """Execute BUSCO quality assessment"""

def run(self):
    """Execute complete gene prediction workflow"""
```

---

### 2. Annotation Module
**File:** `src/modules/annotation.py`

#### Purpose
Functional annotation via FASTA search against UniProt database and gene classification.

#### Key Classes
- `Annotation`: Main annotation workflow class
- `ProteinClassifier`: Classifies proteins as Single or Multi (redundant)

#### Configuration Parameters
```yaml
fasta_search:
  docker_image: "quay.io/biocontainers/fasta3:36.3.8--h470a237_4"
  evalue: 0.1
  threads: 8
```

#### Gene Classification Logic
- **Single**: Unique gene (non-redundant)
- **Multi**: Multiple transcripts mapping to same gene (redundant)

#### Input/Output
- **Input:** `query.pep` from GenePrediction
- **Output:**
  - `annotation_classified.txt` - Classified annotations
  - `multi_genes.fa` / `single_genes.fa` - Classified sequences

---

### 3. Sgd_Annotation Module
**File:** `src/modules/sgd_annotation.py`

#### Purpose
Saccharomyces cerevisiae-specific annotation using SGD database.

#### Current Limitations
- **Yeast-specific**: Only works with S. cerevisiae
- **Known Issue**: Needs generalization for multi-organism support

#### Planned Improvements
- Support for other model organisms
- Configurable database mappings
- Cross-organism gene mapping

#### Key Methods
```python
def __init__(self, config):
    """Initialize with configuration"""

def load_sgd_database(self):
    """Load SGD reference database"""

def map_genes_to_sgd(self, genes):
    """Map gene IDs to SGD identifiers"""

def run(self):
    """Execute SGD annotation workflow"""
```

---

### 4. GuideDesign Module
**File:** `src/modules/guide_design.py`

#### Purpose
CRISPR guide RNA design using CRISPRFlash with off-target specificity analysis.

#### Key Features
- Guide RNA sequence design
- Off-target analysis
- Specificity scoring (MM0, MM1-2)
- Site-specific guide identification

#### Configuration Parameters
```yaml
crispr:
  docker_image: "kazukinakamae/crisflash:1.2.0"
  pam_pattern: "NGG"           # SpCas9 PAM
  max_mismatch: 2
  threads: 10
  memory: "40g"
```

#### Guide RNA Scoring
- **MM0**: Perfect matches in off-target regions
- **MM1-2**: Sequences with 1-2 mismatches
- **Off-target score**: (0-1 scale) - Specificity assessment
- **Site-specific**: MM0=1 AND MM1-2=0 → High specificity

#### Input/Output
- **Input:** Gene sequences (single/multi)
- **Output:**
  - `*_targets.bed` - Guide target positions
  - `*_targets_specific_clean.bed` - Filtered high-specificity guides

---

### 5. GOAnalysis Module
**File:** `src/modules/go_analysis.py`

#### Purpose
Gene Ontology analysis with GO-KFI (Knock-out Feasibility Index) calculation.

#### GO-KFI Scoring
```
GO-KFI Score = (Number of genes with guides) / (Total genes in GO term)

Classification:
- Fully targetable: 100% of genes have guides (KFI = 1.0)
- Partially targetable: Some genes have guides (0 < KFI < 1.0)
- Not targetable: No genes have guides (KFI = 0.0)
```

#### GO Namespaces
- **BP**: Biological Process
- **MF**: Molecular Function
- **CC**: Cellular Component

#### Configuration Parameters
```yaml
go_analysis:
  add_go_image: "kazukinakamae/add-go-columns:v1.0"
  compute_kfi_image: "kazukinakamae/compute_go_kfi:v1.0"
  count_parent_image: "kazukinakamae/count_go_parent_terms:v1.0"
```

#### Docker Containers
1. **add-go-columns**: Adds GO annotations to gene list
2. **compute_go_kfi**: Calculates KFI scores
3. **count_go_parent_terms**: Analyzes GO hierarchy

---

### 6. PathwayAnalysis Module
**File:** `src/modules/pathway_analysis.py`

#### Purpose
Pathway-level analysis using WikiPathways GPML files and PW-KFI scoring.

#### PW-KFI Calculation
```
PW-KFI Score = (Number of genes with guides) / (Total genes in pathway)

Identical logic to GO-KFI but applied at pathway level
```

#### GPML Parsing
- Extracts DataNode elements (genes, metabolites)
- Parses Xref information (SGD IDs)
- Reads Graphics coordinates for visualization
- Identifies Interaction elements (edges)

#### Key Methods
```python
def parse_gpml_file(self, gpml_path):
    """Parse GPML XML file"""
    # Returns: nodes, edges, pathway metadata

def extract_pathway_metadata(self, root):
    """Extract pathway name, organism, description"""

def extract_pathway_genes(self, root):
    """Extract GeneProduct DataNodes"""

def compute_pathway_kfi(self, genes_with_guides, total_genes):
    """Calculate PW-KFI score"""
```

#### Output Files
- `pathway_kfi_table.csv` - KFI scores for all pathways
- `pathway_kfi_fully_targetable.csv` - Fully targetable pathways
- `pathway_gene_mapping.tsv` - Pathway-gene relationships

---

### 7. Report Generation Modules

#### PathwayReportGenerator
**File:** `src/modules/pathway_report_generator.py`

- Generates pathway-specific HTML reports
- Creates SVG visualizations of pathway diagrams
- Color-codes genes based on guide availability
- Includes interactive gene information panels

#### ReportGenerator
**File:** `src/modules/report_generator.py`

- Generates comprehensive HTML dashboard
- Creates interactive data tables
- Implements CSV export functionality
- Generates SQLite database for data persistence
- Fetches GO term information from QuickGO API

---

## Data Flow

### Complete Pipeline Data Flow

```
Transcriptome (input)
        ↓
[GenePrediction]
        ↓
Protein sequences (query.pep)
        ├─→ [Annotation]
        │        ↓
        │   Gene classifications
        │        └─→ [GuideDesign]
        │              ↓
        │         Guide RNAs (annotation with guides)
        │              │
        │    ┌─────────┴─────────┐
        │    │ (Two parallel paths: GO OR Pathway analysis)
        │    │
        │    ├─→ [GOAnalysis]    └─→ [PathwayAnalysis]
        │         ↓                    ↓
        │    GO-KFI scores        PW-KFI scores
        │         ↓                    ↓
        │  [Report Generation]  [Pathway Report Generation]
        │         ↓                    ↓
        │   HTML Dashboard      Pathway HTML Reports
        │         ↓                    ↓
        └─────────┬────────────────────┘
                  ↓
     Combined Output (HTML, CSV, JSON, SQLite)
```

**Key Points:**
- GO Analysis and Pathway Analysis are independent parallel processes
- Each requires separate YAML configuration files
- ReportGenerator handles GO-KFI based reports
- PathwayReportGenerator handles Pathway-KFI based reports
- Both can run independently or together

### Pipeline Processing Flow Details

#### Stage 1: Common Pre-analysis Steps
1. **GenePrediction**
   - Input: Transcriptome FASTA file
   - Process: TransDecoder-based ORF prediction (default) or provided protein sequences (existing protein mode)
   - Output: Predicted protein sequences (query.pep)
   - Note: TransDecoder can be skipped if protein sequences are provided via config

2. **Annotation**
   - Input: Protein sequences from GenePrediction
   - Process: FASTA-based similarity search against Protein Database
   - Output: Annotated gene sequences with classifications (Single/Multi)
   - Database dependency: Protein DB for subsequent analysis type

3. **GuideDesign**
   - Input: Annotated gene sequences
   - Process: CRISPR guide RNA design using CRISPRFlash
   - Output: Guide RNA sequences with specificity scoring
   - Off-target search methods:
     - Transcriptome assembly-based (default)
     - Genome assembly-based (recommended): Provides more comprehensive off-target analysis
   - Note: Genome assembly automatically used when provided in config

#### Stage 2: Parallel Analysis Paths (Choose One or Both)

**Path A: GO-KFI Analysis**
1. GOAnalysis processes annotated sequences with Protein DB for GO Analysis
2. Calculates GO-KFI scores for Gene Ontology terms
3. ReportGenerator creates GO-based HTML dashboard
4. Output: Interactive HTML report with GO term visualization

**Path B: Pathway-KFI Analysis**
1. PathwayAnalysis processes annotated sequences with Protein DB for Pathway Analysis
2. Calculates PW-KFI scores for pathway components
3. PathwayReportGenerator creates pathway-specific HTML reports
4. Output: Pathway diagrams with color-coded guide availability

#### Annotation Database Configuration

The Annotation process requires different Protein Database references depending on the analysis path:

- **GO Analysis Path**: Uses `protein_db_for_goanalysis`
  - Database can be any organism-relevant protein database (e.g., UniProt)
  - Provides annotations for Gene Ontology term mapping
  
- **Pathway Analysis Path**: Uses `sgd_db`
  - Database should contain genes present in target pathway databases
  - Currently optimized for Saccharomyces cerevisiae (SGD)
  - Future versions will support additional organisms

**Important**: Select appropriate protein database for your analysis type and organism.

#### Gene Prediction Flexibility

GenePrediction module supports two operational modes:

1. **TransDecoder Mode (Default)**
   - Performs de novo ORF prediction from transcriptome
   - Recommended for novel transcriptomes without existing annotations
   - Configuration: Leave `protein_sequences` empty or undefined

2. **Existing Protein Mode**
   - Skips TransDecoder; uses provided protein sequences
   - Requires: `protein_sequences` file path and `protein_transcript_mapping` file
   - Mapping file: Tab-separated (TSV) format linking transcript IDs to protein IDs
   - Useful when working with previously predicted or known protein sequences
   - Configuration: Specify both `protein_sequences` and `protein_transcript_mapping`

#### Guide Design Off-Target Search Methods

GuideDesign module provides flexible off-target analysis strategies:

1. **Transcriptome Assembly-Based (Default)**
   - Off-target search limited to transcriptome sequences
   - Lower computational cost
   - Suitable when genome assembly unavailable
   - Configuration: `genome_assembly` not specified

2. **Genome Assembly-Based (Recommended)**
   - Comprehensive off-target search across entire genome
   - Identifies off-target sites in non-transcribed regions
   - Provides more reliable specificity assessment
   - Higher computational cost (multi-gigabyte genome indexing)
   - Configuration: Specify `genome_assembly` file path
   - Automatically used when genome assembly provided
   - Recommended for CRISPR experiments requiring high specificity

### Database Schema

#### genes table
```sql
CREATE TABLE genes (
    id INTEGER PRIMARY KEY,
    query_id TEXT,              -- Original transcript ID
    subject_id TEXT,            -- UniProt accession
    gene_name TEXT,
    classification TEXT,        -- Single/Multi
    unique_id TEXT,
    identity REAL,
    evalue REAL,
    description TEXT,
    guide_count INTEGER
);
```

#### guide_rnas table
```sql
CREATE TABLE guide_rnas (
    id INTEGER PRIMARY KEY,
    gene_unique_id TEXT,
    chromosome TEXT,
    start INTEGER,
    end INTEGER,
    sequence TEXT,              -- 20bp guide sequence
    strand TEXT,
    score REAL,
    pam_sequence TEXT,
    gc_content REAL
);
```

#### go_terms table
```sql
CREATE TABLE go_terms (
    id INTEGER PRIMARY KEY,
    go_id TEXT,                 -- GO:xxxxxxx
    go_name TEXT,
    go_definition TEXT,
    go_namespace TEXT,          -- BP/MF/CC
    total_genes INTEGER,
    genes_with_guides INTEGER,
    kfi_score REAL,
    targetable_type TEXT,       -- Fully/Partially/Not
    annotated_genes TEXT,
    annotated_genes_full TEXT
);
```

---

## API Documentation

### Main Pipeline Class

```python
class tWhEELPipeline:
    """Main orchestrator for tWhEEL pipeline"""
    
    def __init__(self, config_path: str, analysis_name: str = None):
        """
        Initialize pipeline with configuration
        
        Args:
            config_path: YAML configuration file path
            analysis_name: Optional analysis name or directory name
        """
    
    def validate_inputs(self) -> bool:
        """Validate all input files and configuration"""
    
    def run(self) -> bool:
        """Execute entire pipeline"""
    
    def run_step(self, step_name: str) -> bool:
        """
        Execute specific pipeline step
        
        Args:
            step_name: One of:
                - gene_prediction
                - annotation
                - guide_design
                - go_analysis
                - pathway_analysis
                - pathway_report_generation
                - report_generation
        """
```

### Module Base Class

```python
class BaseModule:
    """Base class for all pipeline modules"""
    
    def __init__(self, config: dict):
        """Initialize module with configuration"""
        self.config = config
        self.logger = logging.getLogger(__name__)
    
    def run(self) -> bool:
        """Execute module - must be implemented by subclass"""
        raise NotImplementedError
```

### Configuration Structure

```python
config = {
    'input': {
        'transcriptome': str,           # Always required
        'uniprot_db': str,              # Conditional: required for annotation/GO analysis, NOT used for pathway analysis
        'sgd_db': str,                  # Conditional: required only for pathway analysis (yeast-specific)
        'protein_sequences': str,       # Optional: skip TransDecoder when provided
        'protein_transcript_mapping': str,  # Optional: required when protein_sequences is provided
        'genome_assembly': str          # Optional: improves off-target analysis
    },
    'output': {
        'base_dir': str,
        'intermediate': str
    },
    'execution': {
        'steps': List[str],
        'parallel': bool,
        'verbose': bool
    },
    'transdecoder': {...},
    'busco': {...},
    'fasta_search': {...},
    'crispr': {...},
    'go_analysis': {...},
    'report': {...},
    'database': {...},
    'docker': {...}
}
```

---

## Development Guidelines

### Code Style
- **Language**: Python 3.8+
- **Format**: PEP 8 compliant
- **Docstrings**: Google-style docstrings required
- **Type hints**: Recommended for all functions
- **Logging**: Use `logging` module, not `print()`

### Comment and Documentation Standards
- All functions must have docstrings
- Comments should explain WHY, not WHAT
- Comments must match actual implementation
- Keep docstrings and code synchronized

### Error Handling
- Use try-except blocks appropriately
- Log all errors with full context
- Provide meaningful error messages
- Use custom exceptions where needed

### Testing
- Unit tests for critical functions
- Integration tests for module workflows
- Test data should be representative
- Use pytest framework

### Logging Standards
- **INFO**: Pipeline progress milestones
- **DEBUG**: Detailed execution information
- **WARNING**: Potential issues or fallbacks
- **ERROR**: Critical failures that stop execution

---

## Extension Guide

### Adding a New Analysis Step

1. **Create Module File**
```python
# src/modules/new_analysis.py

from modules.base import BaseModule

class NewAnalysis(BaseModule):
    """New analysis module"""
    
    def __init__(self, config):
        super().__init__(config)
        self.logger.info("Initializing NewAnalysis")
    
    def run(self):
        """Execute analysis"""
        try:
            # Implementation
            self.logger.info("NewAnalysis completed successfully")
            return True
        except Exception as e:
            self.logger.error(f"Error in NewAnalysis: {e}")
            return False
```

2. **Add to twheel.py**
```python
# In twheel.py

from modules.new_analysis import NewAnalysis

# In run_step method:
elif step_name == "new_analysis":
    analyzer = NewAnalysis(self.config)
    analyzer.run()

# In available steps:
steps = ["gene_prediction", ..., "new_analysis", ...]
```

3. **Update Configuration**
```yaml
# config/default.yaml

new_analysis:
  parameter1: value1
  parameter2: value2

execution:
  steps:
    - "gene_prediction"
    - "new_analysis"
```

### Docker Integration

When using Docker containers:

1. **Create Dockerfile**
```dockerfile
FROM ubuntu:20.04
WORKDIR /app
COPY . .
RUN apt-get update && apt-get install -y ...
ENTRYPOINT ["python", "main.py"]
```

2. **Reference in Config**
```yaml
new_analysis:
  docker_image: "myregistry/new_analysis:v1.0"
  parameter1: value1
```

3. **Use in Module**
```python
def run_docker(self):
    """Run analysis in Docker container"""
    docker_image = self.config['new_analysis']['docker_image']
    subprocess.run([
        'docker', 'run', '--rm',
        '-v', f'{input_dir}:/input',
        '-v', f'{output_dir}:/output',
        docker_image,
        '--input', '/input/data.txt',
        '--output', '/output/results.txt'
    ])
```

---

## Code Quality Standards

### Required Documentation
- Module docstrings (top of each file)
- Class docstrings
- Function/method docstrings with:
  - Description
  - Args
  - Returns
  - Raises (if applicable)
- Complex logic comments

### Comment Requirements
- Update comments when code changes
- Comments should align with implementation
- Avoid outdated or misleading comments
- Use TODO/FIXME for known issues

### Performance Considerations
- Avoid unnecessary data copying
- Use generators for large datasets
- Profile code for bottlenecks
- Optimize Docker container sizes
- Cache expensive computations

### Testing Requirements
- >80% code coverage for new code
- Unit tests for utility functions
- Integration tests for workflows
- Edge case testing

### Version Control
- Descriptive commit messages
- Logical commit grouping
- Reference issues in commits
- Keep commits focused and atomic

---

## Known Limitations and Future Improvements

### Current Limitations
1. SGD annotation limited to Saccharomyces cerevisiae
2. Japanese text in some output messages and code comments
3. Unused data generation in report scripts
4. Some script comments don't match current implementation

### Planned Improvements
- Multi-organism support (Phase 2)
- English language standardization (Phase 1.1)
- Report optimization (Phase 1.3)
- Code comment and documentation sync (Phase 1.2)
- Advanced visualization (Phase 3)
- Parallel processing (Phase 4)

---

## References

### External Tools Documentation
- [TransDecoder](https://github.com/TransDecoder/TransDecoder)
- [BUSCO](https://busco.ezlab.org/)
- [FASTA/SSEARCH](https://fasta.bioch.virginia.edu/)
- [CRISPRFlash](https://github.com/tsuryuchka/CRISPRFlash)
- [Gene Ontology](http://geneontology.org/)
- [UniProt](https://www.uniprot.org/)
- [WikiPathways](https://www.wikipathways.org/)

### Python Libraries
- `pyyaml` - Configuration file parsing
- `pandas` - Data manipulation
- `sqlalchemy` - Database ORM
- `requests` - HTTP requests
- `docker` - Docker API client

---

## Document Information

- **Version**: 1.0.0
- **Last Updated**: November 2025
- **Maintainer**: tWhEEL Development Team
- **Status**: Active Development
