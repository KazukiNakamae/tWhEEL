# tWhEEL - Transcriptome Workflow for Hierarchical Exploration and Expression Levels

A comprehensive bioinformatics pipeline for transcriptome functional analysis and CRISPR target identification. tWhEEL integrates gene prediction, functional annotation, guide RNA design, and Gene Ontology analysis to systematically identify genes for CRISPR-based knockout experiments.

## Quick Start

### 1. Installation

```bash
# Clone repository
git clone https://github.com/KazukiNakamae/tWhEEL.git
cd tWhEEL

# Create conda environment
conda create -n twheel python=3.12 -y
conda activate twheel

# Install dependencies
pip install -r requirements.txt

# Build Docker images
docker build -t kazukinakamae/add-go-columns:v1.0 src/docker/add_go_columns/
docker build -t kazukinakamae/compute_go_kfi:v1.0 src/docker/compute_go_kfi/
docker build -t kazukinakamae/count_go_parent_terms:v1.0 src/docker/count_go_parent_terms/
docker build -t kazukinakamae/run_alignment_and_mask:v1.0 src/docker/run_alignment_and_mask/
```

### 2. Prepare Configuration

```bash
# Copy and edit configuration file
cp config/default.yaml config/my_analysis.yaml

# Edit config file with your input files and parameters
# Required: transcriptome and uniprot_db
```

### 3. Run Pipeline

```bash
# Run complete pipeline
python src/twheel.py --config config/my_analysis.yaml --name my_analysis

# Run specific step
python src/twheel.py --config config/my_analysis.yaml --step annotation
```

## Pipeline Overview

```
Transcriptome → Gene Prediction → Annotation → Guide Design → GO Analysis → Report
```

### Core Components

| Component | Function | Output |
|-----------|----------|--------|
| **Gene Prediction** | ORF prediction using TransDecoder | Protein sequences |
| **Annotation** | Functional annotation via FASTA search | Gene classifications, annotations |
| **Guide Design** | CRISPR guide RNA design using CRISPRFlash | Guide sequences with specificity scores |
| **GO Analysis** | Gene Ontology mapping and KFI calculation | GO-KFI scores, targetability classification |
| **Report Generation** | Interactive HTML dashboard with data visualization | Comprehensive analysis report |

## Usage Examples

### Basic Analysis
```bash
python src/twheel.py --config config/default.yaml
```

### Use Existing Protein Sequences
```yaml
# config/protein_mode.yaml
input:
  transcriptome: "/path/to/transcriptome.fasta"
  protein_sequences: "/path/to/proteins.fasta"
  uniprot_db: "/path/to/uniprot_db.fasta"
```

### Resume Previous Analysis
```bash
python src/twheel.py --config config/my_analysis.yaml \
  --name tWhEEL_on_20250918_120000_my_analysis
```

### List Available Steps
```bash
python src/twheel.py --list-steps
```

## Output Files

### Generated Reports
- **report.html** - Interactive web dashboard with tables, charts, and data export
- **summary.txt** - Text summary of analysis results
- **twheel.log** - Detailed execution log

### Data Files
- **intermediate/** - Intermediate analysis files from each pipeline step
- **web_data/report_sources/** - Durable inputs used to regenerate reports after `intermediate/` is removed
- **web_data/report_sources/manifest.json** - Inventory of the staged report inputs
- **web_data/twheel_data.db** - Minimal SQLite cache containing report metadata and GO term labels/definitions
- **web_data/summary.json** - JSON summary data for programmatic access

Large report inputs remain in their native TSV/CSV/FASTA/BED formats. Each run
publishes one complete, versioned source snapshot and atomically switches the
manifest only after every file is ready. tWhEEL uses copy-on-write filesystem
clones where available and falls back to independent copies, so later writes to
`intermediate/` cannot damage a published snapshot. Pathway reports stage their
KFI table, canonical gene mapping, and referenced GPML files in the same source
store; they do not create a separate SQLite result database.

### Key Results
- **annotation_with_guides_final.tsv** - Gene annotations with guide RNAs
- **go_kfi_table.csv** - GO Knock-out Feasibility Index scores
- **pathway_kfi_table.csv** - Pathway-level targetability assessment

## Configuration

### Required Parameters
```yaml
input:
  transcriptome: "path/to/transcriptome.fasta"  # Always required
```

### Conditional Parameters
```yaml
input:
  # For Gene Prediction, Annotation, Guide Design, GO Analysis:
  uniprot_db: "path/to/protein_database.fasta"  # Required for these steps
  
  # For Pathway Analysis only (alternative to uniprot_db):
  sgd_db: "path/to/sgd_reference.fasta"         # Used for pathway analysis
```

### Optional Parameters
```yaml
input:
  protein_sequences: "path/to/proteins.fasta"   # Skip TransDecoder
  genome_assembly: "path/to/genome.fasta"       # Better off-target analysis
  protein_transcript_mapping: "path/to/mapping.tsv"  # Required with protein_sequences
```

### Database Parameter Details
- **uniprot_db**: Any protein database (FASTA format)
  - Not limited to UniProt - can be custom databases
  - Required for: annotation, GO analysis
  - NOT used in: pathway analysis (uses sgd_db instead)
  - Examples: UniProt Swiss-Prot, TrEMBL, or custom organism-specific databases
- **sgd_db**: SGD (Saccharomyces Genome Database) reference
  - Required only for pathway analysis
  - Currently yeast-specific (S. cerevisiae)

### Resource Configuration
```yaml
crispr:
  memory: "40g"        # For large transcriptomes
  threads: 10

busco:
  cpu: 8

fasta_search:
  threads: 8
```

See `config/default.yaml` for all available options.

## Performance Tuning

### For Large Datasets (>100k transcripts)
```yaml
crispr:
  memory: "60g"
  threads: 20
busco:
  cpu: 16
fasta_search:
  threads: 16
```

### For Memory-Constrained Systems (16GB RAM)
```yaml
crispr:
  memory: "12g"
  threads: 6
busco:
  cpu: 4
fasta_search:
  threads: 4
```

## Understanding Results

### GO-KFI Score
The Gene Ontology Knock-out Feasibility Index (GO-KFI) quantifies the proportion of genes in each GO term that have targetable guide RNAs:
- **1.0 (Fully targetable)** - All genes have guide RNAs
- **0.5 (Partially targetable)** - Half of genes have guide RNAs
- **0.0 (Not targetable)** - No genes have guide RNAs

### Guide RNA Specificity
- **MM0 = 1, MM1-2 = 0** - Highly specific (recommended)
- **MM0 > 1** - Non-specific guides
- **Off-target score** - Quantitative specificity (0-1 scale)

### Gene Classification
- **Single** - Non-redundant genes (unique annotation)
- **Multi** - Redundant genes (multiple transcripts per gene)

## Troubleshooting

### Docker Issues
```bash
# Check Docker daemon
docker ps

# Add user to docker group
sudo usermod -aG docker $USER
```

### Memory Issues
```bash
# Check available memory
free -h

# Reduce CRISPR memory allocation in config
crispr:
  memory: "16g"
```

### Input File Format
```bash
# Verify FASTA format
head -2 transcriptome.fasta
# Should show: ">sequence_header"
#             "sequence_data"
```

## Input File Requirements

### Transcriptome
- Format: FASTA (`.fa`, `.fasta`, `.fna`)
- One transcript per sequence
- Valid DNA sequences (ACGTN)
- Recommended: High-quality assembly from assembler like Trinity, SPAdes

### UniProt Database
- Format: FASTA with UniProt headers
- Download from: https://www.uniprot.org/downloads
- Swiss-Prot recommended over TrEMBL for accuracy
- Organism-specific databases acceptable

### Genome Assembly (Optional)
- Format: FASTA
- Required for comprehensive off-target analysis
- Improves guide RNA specificity prediction

## System Requirements

- **Operating System**: Linux/macOS
- **RAM**: 16GB minimum (32GB+ recommended)
- **Storage**: 50GB+ for intermediate files
- **Docker**: Latest version
- **Python**: 3.8+

## Citation

If you use tWhEEL in your research, please cite:

```
Nakamae, K., et al. (Year). tWhEEL: A comprehensive pipeline for 
transcriptome-based CRISPR target prioritization.
```

## Support

For issues, questions, or feature requests:
- Open an issue on GitHub
- Check existing documentation
- Review the twheel.log file for error details

## License

See LICENSE file for details.

## Contributing

Contributions welcome! Please see CONTRIBUTING.md for guidelines.

## Related Resources

- [CRISPR/Cas9 Guide Design Best Practices](https://www.addgene.org/crispr/)
- [Gene Ontology Project](http://geneontology.org/)
- [UniProt Database](https://www.uniprot.org/)
- [TransDecoder](https://github.com/TransDecoder/TransDecoder)
