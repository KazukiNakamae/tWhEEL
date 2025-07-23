# tWhEEL - Transcriptome Workflow for Hierarchical Exploration and Expression Levels

**A comprehensive pipeline for transcriptome functional analysis and CRISPR target identification**

tWhEEL is a comprehensive pipeline that performs gene prediction, functional annotation, CRISPR guide design, and GO functional analysis from transcriptome assembly data to identify target genes for knockout experiments.

## Features

- **Gene Prediction**: ORF prediction using TransDecoder and BUSCO quality assessment
- **Functional Annotation**: FASTA search against UniProt database and gene classification
- **CRISPR Guide Design**: SpCas9 guide RNA design using CrisFlash
- **GO Functional Analysis**: GO Knock-out Feasibility Index (GO-KF) calculation
- **Interactive Reports**: HTML-based web report generation
- **Flexible Configuration**: YAML-based configuration files

## Requirements

- Python 3.8 or higher
- Docker
- Required Python packages (see requirements.txt)

## Installation

```bash
# Clone the repository
git clone https://github.com/KazukiNakamae/tWhEEL.git
cd tWhEEL

# Install required packages
conda create -n twheel_env -y;
conda activate twheel_env;
conda install conda-forge::python=3.12 -y;
python -m pip install -r requirements.txt;

# Build Docker images if needed
sudo docker build -t kazukinakamae/add-go-columns:v1.0 src/docker/add_go_columns/
sudo docker build -t kazukinakamae/compute_go_kfi:v1.0 src/docker/compute_go_kfi/
sudo docker build -t kazukinakamae/count_go_parent_terms:v1.0 src/docker/count_go_parent_terms/
```

## Usage

### 1. Configuration File Setup

Copy `config/default.yaml` and set the input file paths:

```yaml
input:
  transcriptome: "/path/to/transcriptome.fsa_nt"  # Required
  uniprot_db: "/path/to/uniprot_sprot.fasta"  # Required
  genome_assembly: "/path/to/genome.fna"  # Optional - uses transcriptome if not provided
```

**Note**: The `genome_assembly` parameter is optional. If not provided, the pipeline will use the transcriptome for off-target analysis, though this may be less comprehensive than using a genome assembly.

### 2. Execution

```bash
# Run all steps
python src/twheel.py --config config/my_config.yaml

# Run specific step only
python src/twheel.py --config config/my_config.yaml --step annotation

# Check available steps
python src/twheel.py --list-steps
```

### 3. Output Files

After execution, the following files are generated:

- `outputs/summary.txt`: Text summary
- `outputs/report.html`: HTML interactive report
- `outputs/web_data/`: Web application data
- `outputs/intermediate/`: Intermediate files

## Pipeline Structure

### 1. Gene Prediction (gene_prediction)
- **TransDecoder**: ORF prediction and protein sequence generation
- **BUSCO**: Assembly quality assessment
- **seqkit**: Sequence statistics

### 2. Annotation (annotation)
- **FASTA Search**: ggsearch36 against UniProt database
- **Gene Classification**: Multi/Single gene classification
- **Sequence Extraction**: Classified gene sequence extraction

### 3. Guide Design (guide_design)
- **CrisFlash**: SpCas9 guide RNA design
- **Site Specificity**: Specific target extraction and filtering
- **Annotation Integration**: Guide information and annotation integration

### 4. GO Analysis (go_analysis)
- **GO Annotation**: GO term retrieval from UniProt IDs
- **GO-KFI Calculation**: Knockout feasibility index calculation
- **Functional Classification**: Fully/Partially/Not targetable classification

### 5. Report Generation (report_generation)
- **Text Summary**: Execution statistics and summary
- **HTML Report**: Interactive web report
- **Database**: SQLite and JSON format data storage

## Output Files

### Main Outputs
- `summary.txt`: Execution summary
- `report.html`: HTML report
- `web_data/twheel_data.db`: SQLite database
- `web_data/summary.json`: JSON format statistics

### Intermediate Files
- `intermediate/query.pep`: Predicted protein sequences
- `intermediate/annotation_classified.txt`: Classified annotations
- `intermediate/annotation_with_guides_final.tsv`: Guide information integration results
- `intermediate/go_kfi_table.csv`: GO-KFI table
- `intermediate/go_parent_counts.png`: GO hierarchy visualization

## Configuration Options

### Input Settings
```yaml
input:
  genome_assembly: ""  # Genome assembly file
  transcriptome: ""    # Transcriptome file
  uniprot_db: ""       # UniProt database
```

### CRISPR Settings
```yaml
crispr:
  pam_pattern: "NGG"   # PAM pattern
  max_mismatch: 2      # Maximum mismatch count
  threads: 10          # Thread count
```

### GO Analysis Settings
```yaml
go_analysis:
  add_go_image: "kazukinakamae/add-go-columns:v1.0"
  compute_kfi_image: "kazukinakamae/compute_go_kfi:v1.0"
  count_parent_image: "kazukinakamae/count_go_parent_terms:v1.0"
```

## HTML Report Features

The generated HTML report includes the following features:

### Dashboard
- Processing statistics summary
- GO Knock-out Feasibility distribution graph

### Gene Browser
- Multi/Single gene search and filtering
- UniProt ID, identity, E-value display
- Guide RNA count display

### Guide RNA
- List of designed guide RNAs
- Chromosome position, sequence, strand information
- Search and filtering functionality

### GO Analysis
- GO term search and filtering
- KFI score-based classification
- Targetable classification display

## Troubleshooting

### Docker-related Issues
```bash
# If Docker image not found
docker build -t kazukinakamae/add-go-columns:v1.0 src/docker/add_go_columns/

# Permission issues
sudo usermod -aG docker $USER
```

### Memory Issues
```yaml
# Adjust memory usage in configuration file
crispr:
  memory: "20g"  # Reduce memory usage
```

### Performance Tuning
```yaml
# Adjust CPU usage
crispr:
  threads: 4     # Reduce thread count
busco:
  cpu: 4         # Reduce CPU count
```

## References

- **TransDecoder**: Haas et al. (2013) Nature Protocols
- **BUSCO**: Simão et al. (2015) Bioinformatics
- **CrisFlash**: Nakamae et al. (2024) *In preparation*
- **UniProt**: UniProt Consortium (2023) Nucleic Acids Research

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Contributing

We welcome contributions to the project. Please submit bug reports and feature requests through the Issues page, and code improvements through Pull Requests.

## Author

- **Kazuki Nakamae** - [GitHub](https://github.com/KazukiNakamae)

## Changelog

- **v1.0.0** (2025-01-18): Initial release
  - Basic pipeline functionality
  - HTML report generation
  - GO-KFI calculation
