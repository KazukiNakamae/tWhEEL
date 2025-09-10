# tWhEEL - Transcriptome Workflow for Hierarchical Exploration and Expression Levels

**A comprehensive pipeline for transcriptome functional analysis and CRISPR target identification**

tWhEEL is an integrated bioinformatics pipeline that performs gene prediction, functional annotation, CRISPR guide RNA design, and Gene Ontology (GO) functional analysis from transcriptome assembly data. The pipeline calculates the GO Knock-out Feasibility Index (GO-KFI) to systematically identify target genes for CRISPR-based knockout experiments, enabling researchers to prioritize functional studies based on targetability assessments.

## Features

### Core Pipeline Components
- **Gene Prediction**: ORF prediction using TransDecoder with BUSCO quality assessment
- **Functional Annotation**: FASTA search against UniProt database with gene classification
- **CRISPR Guide Design**: SpCas9 guide RNA design using CRISPRFlash with off-target analysis
- **GO Functional Analysis**: GO Knock-out Feasibility Index (GO-KFI) calculation and classification
- **Interactive Reports**: Comprehensive HTML-based web reports with data visualization

### Advanced Features
- **GO-KFI Assessment**: Quantitative measurement of gene knockout feasibility based on GO term targetability
- **Multi-level Classification**: Genes classified as Multi (redundant) or Single (non-redundant) based on annotation
- **Site-specific Targeting**: Identification of high-specificity guide RNAs with minimal off-target effects
- **Interactive Visualization**: Dynamic charts, tables, and data export functionality
- **Database Integration**: SQLite database with comprehensive data storage and retrieval
- **Flexible Configuration**: YAML-based configuration with extensive customization options

## Architecture

### Pipeline Overview
```
Transcriptome Assembly → Gene Prediction → Functional Annotation → Guide Design → GO Analysis → Report Generation
                     ↓                  ↓                      ↓             ↓              ↓
                 Protein Sequences   Gene Classification   Guide RNAs    GO-KFI Scores  HTML Report
```

### Module Structure
- **`twheel.py`**: Main pipeline orchestrator with step execution control
- **`gene_prediction.py`**: TransDecoder and BUSCO integration for ORF prediction
- **`annotation.py`**: UniProt FASTA search and gene classification system
- **`guide_design.py`**: CRISPRFlash integration with target specificity analysis
- **`go_analysis.py`**: GO annotation retrieval and KFI calculation
- **`report_generator.py`**: HTML report generation with interactive features

### Docker Components
- **`add-go-columns`**: GO term annotation from UniProt IDs
- **`compute_go_kfi`**: GO-KFI score calculation and classification
- **`count_go_parent_terms`**: GO hierarchy analysis and visualization
- **`run_alignment_and_mask`**: Sequence alignment and masking for multi-gene analysis

### Data Flow
1. **Input Processing**: Transcriptome and reference database preparation
2. **Gene Prediction**: ORF identification and protein sequence generation
3. **Annotation Pipeline**: Functional annotation and gene classification
4. **Guide Design**: CRISPR target identification and specificity scoring
5. **GO Analysis**: Functional categorization and targetability assessment
6. **Report Generation**: Data integration and visualization

## Requirements

### System Requirements
- **Operating System**: Linux/macOS with Docker support
- **Memory**: Minimum 16GB RAM (32GB+ recommended for large datasets)
- **Storage**: 50GB+ free space for intermediate files and databases
- **CPU**: Multi-core processor (8+ cores recommended)

### Software Dependencies
- **Python**: 3.8 or higher
- **Docker**: Latest version with daemon running
- **Conda**: For environment management (recommended)

## Installation

### 1. Repository Setup
```bash
# Clone the repository
git clone https://github.com/KazukiNakamae/tWhEEL.git
cd tWhEEL

# Verify directory structure
ls -la src/modules/  # Should show 5 Python modules
ls -la src/docker/   # Should show 4 Docker directories
```

### 2. Environment Setup
```bash
# Create and activate conda environment
conda create -n twheel_env python=3.12 -y
conda activate twheel_env

# Install Python dependencies
python -m pip install -r requirements.txt

# Verify installation
python -c "import yaml, pandas, numpy, biopython; print('Dependencies installed successfully')"
```

### 3. Docker Image Building
```bash
# Build all required Docker images
sudo docker build -t kazukinakamae/add-go-columns:v1.0 src/docker/add_go_columns/
sudo docker build -t kazukinakamae/compute_go_kfi:v1.0 src/docker/compute_go_kfi/
sudo docker build -t kazukinakamae/count_go_parent_terms:v1.0 src/docker/count_go_parent_terms/
sudo docker build -t kazukinakamae/run_alignment_and_mask:v1.0 src/docker/run_alignment_and_mask/

# Verify Docker images
docker images | grep kazukinakamae
```

### 4. External Tool Images
```bash
# Pull required biocontainer images
docker pull quay.io/biocontainers/transdecoder:5.7.1--pl5321hdfd78af_0
docker pull ezlabgva/busco:v5.2.2_cv1
docker pull quay.io/biocontainers/fasta3:36.3.8--h470a237_4
docker pull kazukinakamae/crisflash:1.2.0
```

### 5. Installation Verification
```bash
# Test pipeline installation
python src/twheel.py --version
python src/twheel.py --list-steps

# Test Docker integration
docker run --rm kazukinakamae/add-go-columns:v1.0 --help
```

### Common Installation Issues

#### Docker Permission Issues
```bash
# Add user to docker group
sudo usermod -aG docker $USER
# Logout and login again, or restart terminal

# Alternative: Use sudo for docker commands
sudo docker build -t image_name path/
```

#### Memory Issues
```bash
# Check available memory
free -h

# Adjust Docker memory limits
# Edit /etc/docker/daemon.json:
{
  "default-runtime": "runc",
  "default-shm-size": "1G",
  "storage-opts": ["size=50G"]
}
```

#### Python Package Conflicts
```bash
# Create clean environment
conda create -n twheel_clean python=3.12 -y
conda activate twheel_clean
pip install --no-cache-dir -r requirements.txt
```

## Configuration

### Configuration File Structure
The pipeline uses YAML configuration files for all settings. Copy and modify `config/default.yaml`:

```bash
cp config/default.yaml config/my_analysis.yaml
```

### Input File Configuration

#### Required Parameters
```yaml
input:
  transcriptome: "/path/to/transcriptome.fasta"    # FASTA format transcriptome assembly (always required)
  uniprot_db: "/path/to/uniprot_sprot.fasta"       # UniProt Swiss-Prot database (always required)
```

#### Optional Parameters
```yaml
input:
  protein_sequences: "/path/to/proteins.fasta"     # Existing amino acid sequences (skips TransDecoder)
  genome_assembly: "/path/to/genome.fasta"         # Genome assembly for comprehensive off-target analysis
```

#### Processing Modes
The pipeline supports two processing modes, but **transcriptome is always required**:

**1. Standard TransDecoder Mode (Default)**
- Uses transcriptome assembly for ORF prediction
- Requires `transcriptome` parameter
- Runs TransDecoder → BUSCO → Functional annotation
- Transcriptome data used throughout downstream analysis

**2. Existing Protein Mode**
- Uses pre-computed amino acid sequences + transcriptome data
- Requires both `transcriptome` and `protein_sequences` parameters
- Skips TransDecoder → BUSCO → Functional annotation
- Transcriptome data still used for guide design and other analyses
- Useful when ORF prediction has already been performed

**Important**: tWhEEL is a transcriptome-based pipeline. The transcriptome file is always required as it provides essential information for guide RNA design, sequence extraction, and downstream analyses, even when using existing protein sequences.

**Note**: If `genome_assembly` is not provided, the transcriptome will be used for off-target analysis, which may be less comprehensive but still functional.

### Pipeline Step Configuration

#### Gene Prediction Settings
```yaml
transdecoder:
  docker_image: "quay.io/biocontainers/transdecoder:5.7.1--pl5321hdfd78af_0"
  min_protein_len: 85                              # Minimum protein length (amino acids)

busco:
  docker_image: "ezlabgva/busco:v5.2.2_cv1"
  lineage: "metazoa_odb10"                         # BUSCO lineage database
  cpu: 8                                           # CPU cores for BUSCO analysis
```

#### Annotation Settings
```yaml
fasta_search:
  docker_image: "quay.io/biocontainers/fasta3:36.3.8--h470a237_4"
  evalue: 0.1                                      # E-value threshold for FASTA search
  threads: 8                                       # Number of search threads
```

#### CRISPR Guide Design Settings
```yaml
crispr:
  docker_image: "kazukinakamae/crisflash:1.2.0"
  pam_pattern: "NGG"                               # PAM sequence pattern (SpCas9)
  max_mismatch: 2                                  # Maximum allowed mismatches
  threads: 10                                      # Number of processing threads
  memory: "40g"                                    # Memory allocation for CRISPRFlash
```

#### GO Analysis Settings
```yaml
go_analysis:
  add_go_image: "kazukinakamae/add-go-columns:v1.0"
  compute_kfi_image: "kazukinakamae/compute_go_kfi:v1.0"
  count_parent_image: "kazukinakamae/count_go_parent_terms:v1.0"
```

#### Report Generation Settings
```yaml
report:
  generate_html: true                              # Generate HTML interactive report
  generate_summary: true                           # Generate text summary
  web_port: 5000                                   # Port for local web server (if needed)

database:
  force_rebuild: false                             # Force database recreation
```

#### System Resource Settings
```yaml
docker:
  max_memory: "40g"                                # Maximum Docker memory usage
  use_docker: true                                 # Enable Docker execution

execution:
  steps:                                           # Pipeline steps to execute
    - "gene_prediction"
    - "annotation"
    - "guide_design"
    - "go_analysis"
    - "report_generation"
  parallel: false                                  # Parallel step execution (future feature)
  verbose: true                                    # Verbose logging
```

### Performance Tuning Configurations

#### High-Performance Configuration (32GB+ RAM)
```yaml
crispr:
  memory: "40g"
  threads: 16
busco:
  cpu: 12
fasta_search:
  threads: 12
```

#### Memory-Constrained Configuration (16GB RAM)
```yaml
crispr:
  memory: "12g"
  threads: 6
busco:
  cpu: 4
fasta_search:
  threads: 4
```

#### Large Dataset Configuration
```yaml
docker:
  max_memory: "60g"
crispr:
  memory: "50g"
  threads: 20
database:
  force_rebuild: true  # Ensure clean database for large datasets
```

## Usage

### Basic Usage
```bash
# Run complete pipeline
python src/twheel.py --config config/my_analysis.yaml

# Run with custom analysis name
python src/twheel.py --config config/my_analysis.yaml --name my_species_analysis
```

### Step-by-Step Execution
```bash
# Run individual steps
python src/twheel.py --config config/my_analysis.yaml --step gene_prediction
python src/twheel.py --config config/my_analysis.yaml --step annotation
python src/twheel.py --config config/my_analysis.yaml --step guide_design
python src/twheel.py --config config/my_analysis.yaml --step go_analysis
python src/twheel.py --config config/my_analysis.yaml --step report_generation
```

### Advanced Usage
```bash
# Resume analysis in existing directory
python src/twheel.py --config config/my_analysis.yaml --name tWhEEL_on_20250719_082938_E_okinawaensis_denovoanalysis

# List available steps
python src/twheel.py --list-steps

# Check version
python src/twheel.py --version
```

## Output Files

### Directory Structure
```
outputs/
└── tWhEEL_on_[timestamp]_[analysis_name]/
    ├── summary.txt                    # Text summary report
    ├── report.html                    # Interactive HTML report
    ├── twheel.log                     # Pipeline execution log
    ├── intermediate/                  # Intermediate analysis files
    ├── web_data/                      # Database and JSON files
    └── static/                        # Web assets (CSS, JS, images)
```

### Main Output Files

#### Summary Report (`summary.txt`)
- **Content**: Execution statistics, processing summary, configuration settings
- **Format**: Plain text with structured sections
- **Use case**: Quick overview of pipeline results

#### Interactive HTML Report (`report.html`)
- **Content**: Complete interactive analysis dashboard
- **Format**: Self-contained HTML with embedded CSS/JavaScript
- **Use case**: Detailed analysis exploration and data export

#### Execution Log (`twheel.log`)
- **Content**: Detailed pipeline execution information
- **Format**: Timestamped log entries with module information
- **Use case**: Debugging and execution monitoring

### Database Files

#### SQLite Database (`web_data/twheel_data.db`)
**Schema Structure:**
```sql
-- Gene information table
genes (
    id INTEGER PRIMARY KEY,
    query_id TEXT,                     -- Original transcript ID
    subject_id TEXT,                   -- UniProt accession
    gene_name TEXT,                    -- Gene symbol
    classification TEXT,               -- Multi/Single classification
    unique_id TEXT,                    -- Internal gene group ID
    identity REAL,                     -- Sequence identity (%)
    evalue REAL,                       -- FASTA search E-value
    description TEXT,                  -- Protein description
    guide_count INTEGER                -- Number of guide RNAs
);

-- Guide RNA information table
guide_rnas (
    id INTEGER PRIMARY KEY,
    gene_unique_id TEXT,               -- Reference to genes.unique_id
    chromosome TEXT,                   -- Target sequence ID
    start INTEGER,                     -- Guide RNA start position
    end INTEGER,                       -- Guide RNA end position
    sequence TEXT,                     -- Guide RNA sequence (20bp)
    strand TEXT,                       -- Target strand (+/-)
    score REAL,                        -- CRISPRFlash score
    pam_sequence TEXT,                 -- PAM sequence
    gc_content REAL                    -- GC content (%)
);

-- GO term information table
go_terms (
    id INTEGER PRIMARY KEY,
    go_id TEXT,                        -- GO accession (GO:xxxxxxx)
    go_name TEXT,                      -- GO term name
    go_definition TEXT,                -- GO term definition
    go_namespace TEXT,                 -- GO namespace (BP/MF/CC)
    total_genes INTEGER,               -- Total genes with this GO term
    genes_with_guides INTEGER,         -- Genes with guide RNAs
    kfi_score REAL,                    -- GO-KFI score (0-1)
    targetable_type TEXT,              -- Fully/Partially/Not targetable
    annotated_genes TEXT,              -- Gene list (formatted)
    annotated_genes_full TEXT          -- Complete gene list
);

-- Gene-GO relationships table
gene_go_relations (
    id INTEGER PRIMARY KEY,
    gene_unique_id TEXT,               -- Reference to genes.unique_id
    go_id TEXT                         -- Reference to go_terms.go_id
);
```

#### JSON Summary (`web_data/summary.json`)
```json
{
  "timestamp": "2025-01-XX",
  "statistics": {
    "gene_prediction": {
      "transcript_count": 50000,        // Total transcripts
      "protein_count": 25000,           // Predicted proteins
      "protein_stats": { ... }          // Protein statistics
    },
    "annotation": {
      "total_genes": 20000,             // Total annotated genes
      "multi_genes": 15000,             // Multi (redundant) genes
      "single_genes": 5000              // Single (non-redundant) genes
    },
    "guide_design": {
      "total_guides": 180000,           // Total guide RNAs
      "site_specific_guides": 45000     // High-specificity guides
    },
    "go_analysis": {
      "total_gos": 8500,                // Total GO terms
      "fully_targetable_percent": 65.2, // Fully targetable (%)
      "chart_data": { ... }             // Visualization data
    }
  }
}
```

### Intermediate Files

#### Gene Prediction Outputs
```
intermediate/
├── transcriptome_formatted.fa        # Formatted input transcriptome
├── query.pep                         # Predicted protein sequences
├── transdecoder_output/              # TransDecoder results
└── busco_output/                     # BUSCO quality assessment
```

#### Annotation Outputs
```
intermediate/
├── fasta_search_results.txt          # Raw FASTA search output
├── fasta_search_parsed.txt           # Parsed search results
├── annotation_classified.txt         # Gene classification results
├── multi_genes.fa                    # Multi-gene sequences
├── single_genes.fa                   # Single-gene sequences
└── aligned_multi_genes/              # Aligned multi-gene sequences
```

#### Guide Design Outputs
```
intermediate/
├── multi_genes_targets.bed           # Multi-gene guide targets
├── single_genes_targets.bed          # Single-gene guide targets
├── multi_genes_targets_specific.bed  # Filtered specific targets
├── single_genes_targets_specific.bed # Filtered specific targets
└── annotation_with_guides_final.tsv  # Final integrated results
```

#### GO Analysis Outputs
```
intermediate/
├── annotation_with_guides_go.tsv     # GO-annotated results
├── go_kfi_table.csv                  # GO-KFI calculation results
├── go_kfi_stat.csv                   # GO-KFI statistics
├── go_kfi_table_fully_targetable.csv # Fully targetable GO terms
├── go_reverse_lookup.tsv             # GO-gene mapping
├── go_parent_counts.csv              # GO hierarchy analysis
└── go_parent_counts.png              # GO distribution visualization
```

### Data File Formats

#### TSV Files (Tab-separated values)
- **Encoding**: UTF-8
- **Delimiter**: Tab character (\t)
- **Headers**: First row contains column names
- **Usage**: Gene annotations, guide RNA information

#### CSV Files (Comma-separated values)
- **Encoding**: UTF-8
- **Delimiter**: Comma (,)
- **Quoting**: Double quotes for fields containing commas
- **Usage**: Statistical summaries, GO analysis results

## HTML Report Features

The generated HTML report provides a comprehensive interactive interface for exploring analysis results.

### Dashboard Overview

#### General Statistics Table
- **Total Transcripts**: Number of input transcript sequences
- **Protein Annotated**: Successfully predicted and annotated proteins
- **Gene Classification**: Multi (redundant) vs Single (non-redundant) genes
- **GO Terms**: Total number of associated GO terms
- **Targetability**: Percentage of fully targetable GO terms

#### Navigation Sidebar
- **Responsive design**: Adapts to different screen sizes
- **Section jumping**: Direct navigation to specific analysis sections
- **Progress tracking**: Visual indication of current section
- **Collapsible**: Mobile-friendly sidebar collapse functionality

### Gene Prediction Section

#### Protein Statistics Dashboard
```
┌─────────────────────┬─────────────────────┬─────────────────────┐
│  Total Predicted    │   Average Length    │   Longest Protein   │
│     Proteins        │       (AA)          │       (AA)          │
├─────────────────────┼─────────────────────┼─────────────────────┤
│      25,000         │        342          │       2,847         │
└─────────────────────┴─────────────────────┴─────────────────────┘
```

#### Protein Sequence Viewer
- **Pagination**: Configurable page sizes (10, 50, 100, 200 entries)
- **Search functionality**: Real-time protein header search
- **Length filtering**: Filter by protein length categories
  - Short: ≤100 amino acids
  - Medium: 101-999 amino acids  
  - Long: ≥1000 amino acids
- **Expandable entries**: Click to view detailed protein information
- **Color coding**: Length-based visual indicators

### Gene Annotation Section

#### Interactive Gene Table
```
Search: [________________] Classification: [All ▼] Guide Status: [All ▼] [📥 Download CSV]

┌──────────────┬──────────────┬──────────────┬──────────────────────┬──────────────┬──────────────┐
│   Query ID   │ UniProt Name │  Gene Name   │     Description      │Classification│# of Targets  │
├──────────────┼──────────────┼──────────────┼──────────────────────┼──────────────┼──────────────┤
│ TRINITY_001  │   P12345     │   ACTB       │ Actin beta protein   │    Multi     │      15      │
│ TRINITY_002  │   Q67890     │   GAPDH      │ Glyceraldehyde-3-... │   Single     │       8      │
└──────────────┴──────────────┴──────────────┴──────────────────────┴──────────────┴──────────────┘

Showing 1-25 of 20,000 entries                                    Page 1 of 800 [Previous] [Next]
```

#### Search and Filtering Features
- **Real-time search**: Query ID, gene name, description, classification
- **Multi-level filtering**: 
  - Classification filter: Multi/Single genes
  - Guide RNA status: With/Without guide RNAs
- **Sortable columns**: Click column headers to sort
- **Export functionality**: Download filtered results as CSV

### Guide RNA Design Section

#### Guide RNA Analysis Table
```
Search: [________________] MM0: [All ▼] MM1-2: [All ▼] [📥 Download CSV]

┌──────────────┬──────────────┬─────────┬─────────┬─────────────────────┬─────────┬─────┬───────┬─────────────┐
│   Query ID   │  Reference   │  Start  │   End   │  Target Sequence    │ Strand  │ MM0 │ MM1-2 │Off-target   │
├──────────────┼──────────────┼─────────┼─────────┼─────────────────────┼─────────┼─────┼───────┼─────────────┤
│ TRINITY_001  │ Multi_ACTB   │  1250   │  1269   │ AGCTGACGTGACGTGACGTG│    +    │  1  │   0   │  1.000000   │
│ TRINITY_002  │Single_GAPDH  │   890   │   909   │ TGACGTGACGTGACGTGAC │    -    │  0  │   1   │  0.875000   │
└──────────────┴──────────────┴─────────┴─────────┴─────────────────────┴─────────┴─────┴───────┴─────────────┘
```

#### Guide RNA Features
- **Sequence display**: Monospace font for guide RNA sequences
- **Specificity indicators**: MM0 (perfect matches) and MM1-2 (1-2 mismatches) counts
- **Off-target scoring**: Quantitative specificity assessment
- **Strand indication**: Visual distinction between +/- strands
- **Position mapping**: Precise genomic coordinates

#### Filtering Options
- **Mismatch filtering**: Filter by MM0 and MM1-2 values
- **Sequence search**: Find specific guide sequences
- **Specificity thresholds**: Filter by off-target scores

### GO-KFI Assessment Section

#### GO Targetability Statistics
```
┌─────────────────────┬─────────────────────┬─────────────────────┬─────────────────────┬─────────────────────┐
│   Fully Targetable │ Simply Targetable   │Multiply Targetable  │Partially Targetable │  Not Targetable     │
│    GO Terms (%)     │   GO Terms (%)      │   GO Terms (%)      │   GO Terms (%)      │   GO Terms (%)      │
├─────────────────────┼─────────────────────┼─────────────────────┼─────────────────────┼─────────────────────┤
│       65.2%         │       45.8%         │       19.4%         │       23.1%         │       11.7%         │
└─────────────────────┴─────────────────────┴─────────────────────┴─────────────────────┴─────────────────────┘
```

#### GO-KFI Distribution Chart
- **Horizontal bar chart**: Visual representation of GO targetability
- **Color coding**: 
  - Dark green: Fully targetable
  - Light green: Simply targetable  
  - Yellow: Multiply targetable
  - Light yellow: Partially targetable
  - Red: Not targetable
- **Interactive elements**: Hover tooltips with detailed statistics
- **Statistical overlay**: F:X [S:X, M:X, P:X, N:X, n:X] format display

#### GO Terms Analysis Table
```
Search: [________________] GO-KFI Class: [All ▼] [📥 Download CSV]

┌─────────────┬────────────────────┬────────────────────┬─────────┬─────────┬─────────┬─────────┬─────────────┬──────────────────┬─────────────┐
│    GO ID    │     GO Term        │   GO Definition    │QuickGO  │GO-ann.  │Target.  │ GO-KFI  │GO-KFI Class │GO-ann. Gene Grps │Guide RNA CSV│ 
├─────────────┼────────────────────┼────────────────────┼─────────┼─────────┼─────────┼─────────┼─────────────┼──────────────────┼─────────────┤
│ GO:0006412  │ translation        │ The cellular...    │QuickGO  │   245   │   198   │ 0.8082  │   Fully     │ ACTB, GAPDH...📋 │  📥 CSV     │
│ GO:0005515  │ protein binding    │ Interacting...     │QuickGO  │   1205  │   456   │ 0.3785  │ Partially   │ MYH1, TUBB...📋  │  📥 CSV     │
└─────────────┴────────────────────┴────────────────────┴─────────┴─────────┴─────────┴─────────┴─────────────┴──────────────────┴─────────────┘
```

#### GO Analysis Features
- **QuickGO integration**: Direct links to GO term databases
- **KFI scoring**: Quantitative knockout feasibility (0-1 scale)
- **Classification system**: Fully/Partially/Not targetable categories
- **Expandable gene lists**: Click 📋 to view complete gene lists
- **Individual CSV export**: Download guide RNAs for specific GO terms

#### GO Parent Terms Distribution
- **Hierarchy visualization**: GO parent term distribution charts
- **Treemap display**: Proportional representation of GO categories
- **Statistical table**: Parent term counts and descriptions
- **Export functionality**: Download parent term analysis

### Data Export and Download Features

#### CSV Export Options
1. **Gene Annotation CSV**: Complete gene annotation data with classifications
2. **Guide RNA CSV**: All guide RNA sequences with targeting information
3. **GO Analysis CSV**: GO terms with KFI scores and gene associations
4. **GO Parent Terms CSV**: GO hierarchy analysis results
5. **GO-specific Guide RNA CSV**: Guide RNAs for individual GO terms

#### Export Features
- **Filtered exports**: Download only currently filtered data
- **Complete datasets**: Export all data regardless of current filters
- **Standardized format**: Consistent CSV formatting with headers
- **UTF-8 encoding**: Proper character encoding for international characters
- **Automatic filename generation**: Timestamped filenames for organization

### Mobile and Responsive Design

#### Responsive Features
- **Mobile navigation**: Collapsible sidebar for mobile devices
- **Responsive tables**: Horizontal scrolling for large tables
- **Touch-friendly**: Optimized button sizes and spacing
- **Adaptive layouts**: Dynamic column adjustment
- **Performance optimization**: Efficient rendering on mobile devices

#### Cross-browser Compatibility
- **Modern browsers**: Chrome, Firefox, Safari, Edge support
- **JavaScript fallbacks**: Graceful degradation for older browsers
- **CSS compatibility**: Cross-browser styling consistency
- **Print support**: Optimized layouts for printing

## Usage Examples and Best Practices

### Example 1: Small Transcriptome Analysis (< 50,000 transcripts)
```yaml
# config/small_dataset.yaml
crispr:
  memory: "16g"
  threads: 8
busco:
  cpu: 6
fasta_search:
  threads: 6
```

```bash
python src/twheel.py --config config/small_dataset.yaml --name small_species_analysis
```

### Example Usage: Comparing TransDecoder vs Existing Proteins
```bash
# Method 1: Traditional TransDecoder approach
python src/twheel.py --config config/transcriptome_analysis.yaml --name transcriptome_mode

# Method 2: Using existing protein sequences (faster)
python src/twheel.py --config config/existing_protein.yaml --name protein_mode

# Compare results between both approaches
ls -la outputs/tWhEEL_on_*_transcriptome_mode/
ls -la outputs/tWhEEL_on_*_protein_mode/
```

### Example 2: Large Transcriptome Analysis (> 100,000 transcripts)
```yaml
# config/large_dataset.yaml
crispr:
  memory: "60g"
  threads: 20
busco:
  cpu: 16
fasta_search:
  threads: 16
docker:
  max_memory: "64g"
```

### Example 3: Using Existing Protein Sequences
```yaml
# config/existing_protein.yaml
input:
  protein_sequences: "/path/to/existing_proteins.fasta"  # Pre-computed protein sequences
  uniprot_db: "/path/to/uniprot_sprot.fasta"
  genome_assembly: "/path/to/genome.fasta"              # Optional

transdecoder:
  skip_if_protein_provided: true  # Skip TransDecoder when protein_sequences is provided

# Continue with normal CRISPR and GO analysis settings
crispr:
  memory: "40g"
  threads: 10
```

```bash
python src/twheel.py --config config/existing_protein.yaml --name protein_analysis
```

### Example 4: Memory-Constrained Environment
```yaml
# config/low_memory.yaml
crispr:
  memory: "8g"
  threads: 4
busco:
  cpu: 2
fasta_search:
  threads: 2
```

### Best Practices

#### Performance Optimization
1. **Memory allocation**: Set CRISPR memory to 60-80% of available RAM
2. **Thread management**: Use 1-2 threads per CPU core
3. **Storage optimization**: Use SSD storage for intermediate files
4. **Docker optimization**: Allocate sufficient Docker memory limits

#### Large Dataset Handling
1. **Staged execution**: Run pipeline steps individually for large datasets
2. **Disk space monitoring**: Ensure 3x transcriptome size in free space
3. **Database optimization**: Enable force_rebuild for consistent results
4. **Log monitoring**: Check logs for memory or disk space warnings

#### Quality Control
1. **Input validation**: Verify FASTA format and file integrity
2. **BUSCO assessment**: Check assembly completeness before proceeding
3. **Guide RNA evaluation**: Review off-target scores and specificity
4. **GO annotation coverage**: Verify GO term annotation completeness

## Technical Specifications

### Module APIs and Interfaces

#### Pipeline Class Interface
```python
class tWhEELPipeline:
    def __init__(self, config_path: str, analysis_name: str = None)
    def run(self) -> bool
    def run_step(self, step_name: str) -> bool
    def validate_inputs(self) -> bool
```

#### Module Base Interface
```python
class BaseModule:
    def __init__(self, config: dict)
    def run(self) -> bool
    def setup_logging(self) -> None
```

### Data Format Specifications

#### FASTA Files
- **Format**: Standard FASTA with header lines starting with '>'
- **Encoding**: UTF-8 text encoding
- **Line endings**: Unix-style line endings (\n)
- **Sequence format**: One-line or multi-line sequences supported

#### Configuration Files
- **Format**: YAML version 1.2
- **Encoding**: UTF-8 with optional BOM
- **Structure**: Nested dictionaries with type validation
- **Comments**: Hash-prefixed comments supported

### Extension Guidelines

#### Adding New Analysis Steps
1. Create module in `src/modules/` following base class pattern
