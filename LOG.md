# tWhEEL Development Log

## Version History

### v1.0.0 (Current Release)

**Release Date:** November 2025

#### Major Features
- Complete transcriptome-based pipeline for CRISPR target identification
- Gene prediction with TransDecoder and BUSCO quality assessment
- Functional annotation via FASTA search against UniProt database
- CRISPR guide RNA design using CRISPRFlash with off-target analysis
- Gene Ontology functional analysis with GO-KFI (GO Knock-out Feasibility Index) scoring
- Pathway analysis integration with interactive visualization
- Interactive HTML reports with data export functionality
- SQLite database backend for data storage and retrieval

#### Core Modules
1. **gene_prediction.py** - ORF prediction and protein sequence generation
2. **annotation.py** - Functional annotation and gene classification
3. **sgd_annotation.py** - Saccharomyces cerevisiae-specific annotation (currently yeast-limited)
4. **guide_design.py** - CRISPR guide RNA design and specificity analysis
5. **go_analysis.py** - Gene Ontology analysis and KFI calculation
6. **pathway_analysis.py** - Pathway-level analysis and visualization
7. **pathway_report_generator.py** - Pathway-specific interactive reports
8. **report_generator.py** - Comprehensive HTML report generation

#### Key Features
- Support for both TransDecoder mode and existing protein sequences mode
- Optional genome assembly for improved off-target analysis
- Multi-level gene classification (Single/Multi redundancy)
- Site-specific guide RNA identification with high specificity
- Interactive HTML dashboards with visualization
- Multiple export formats (CSV, TSV, JSON)
- Comprehensive logging and error handling

#### Known Issues
- SGD annotation module is currently limited to Saccharomyces cerevisiae
- Some command-line output and report text still in Japanese (being migrated to English)
- Script comments occasionally do not match current implementation
- HTMLreport generation includes unused data processing steps

#### Bug Fixes
- Fixed PW-KFI score calculation to use total genes in pathway as denominator
- Corrected Pathway Targetability Distribution from pie chart to horizontal bar chart
- Fixed pathway diagram SVG rendering in HTML reports

### v0.9.0 (Previous Development)

**Release Date:** September 2025

#### Features
- Initial pathway analysis implementation
- Pathway report generation framework
- SVG-based pathway visualization
- GO parent term analysis and distribution

#### Improvements
- Enhanced error handling in pipeline execution
- Improved Docker integration
- Better resource management for large datasets

### v0.8.0

- Basic GO analysis functionality
- GO-KFI scoring implementation

### v0.7.0

- CRISPR guide design integration
- CRISPRFlash Docker container support

### v0.6.0

- Annotation pipeline development
- Gene classification system

### v0.5.0

- Gene prediction module
- TransDecoder integration

## Roadmap

### Phase 1: Code Standardization (Next Priority)

#### 1.1 Internationalization - English Language Migration
- **Status:** In Progress
- **Scope:** Convert all user-facing output to English
  - Command-line logging messages
  - HTML report interface text
  - Error messages and warnings
  - Configuration file comments
- **Timeline:** v1.1.0
- **Impact:** Improved accessibility for international users

#### 1.2 Code Quality - Comment and Documentation Sync
- **Status:** Planned
- **Scope:** Align script comments with actual implementation
  - Review all Python module docstrings
  - Update parameter descriptions to match current behavior
  - Fix outdated comments and documentation
  - Add missing implementation details
- **Timeline:** v1.1.0
- **Impact:** Improved code maintainability and developer experience

#### 1.3 Report Optimization - Data Generation Cleanup
- **Status:** Planned
- **Scope:** Remove unused data generation from reporting scripts
  - Identify unused data structures in pathway_report_generator.py
  - Remove unused data processing steps
  - Optimize memory usage
  - Reduce report generation time
- **Timeline:** v1.1.0
- **Impact:** Improved performance and reduced resource consumption

### Phase 2: Organism Support Expansion

#### 2.1 Generalized SGD Annotation Module
- **Status:** Planned
- **Scope:** Extend SGD annotation to support multiple organisms
  - Refactor sgd_annotation.py for general organism support
  - Implement organism-specific database mappings
  - Support for other model organisms:
    - Schizosaccharomyces pombe
    - Drosophila melanogaster
    - Caenorhabditis elegans
    - Custom organism databases
  - Configuration-based organism selection
- **Timeline:** v1.2.0 (Q1 2026)
- **Impact:** Enables pathway analysis for diverse organisms

#### 2.2 Multi-Organism Database Support
- Add support for OrthoDB for cross-organism gene mapping
- Implement model organism-specific databases
- Support for custom genomic annotation sources

### Phase 3: Feature Enhancement

#### 3.1 Advanced Guide RNA Analysis
- Machine learning-based guide efficiency prediction
- Thermodynamic property calculation
- Secondary structure analysis

#### 3.2 Comparative Analysis Tools
- Multi-species comparison functionality
- Ortholog identification and visualization
- Cross-species pathway conservation analysis

#### 3.3 Enhanced Visualization
- 3D pathway visualization
- Interactive gene network diagrams
- Advanced filtering and search capabilities

### Phase 4: Performance Optimization

#### 4.1 Parallel Processing
- Multi-threaded processing for independent steps
- GPU acceleration for computationally intensive tasks
- Distributed processing support for large datasets

#### 4.2 Caching and Indexing
- Result caching for repeated analyses
- Database indexing optimization
- Incremental update capabilities

## Contributing

We welcome contributions! Please see CONTRIBUTING.md for guidelines.

## Support

For bug reports, feature requests, or questions:
- Open an issue on GitHub
- Provide detailed error messages and log files
- Include configuration and input file information

## License

See LICENSE file for details.

## Acknowledgments

tWhEEL integrates several open-source tools:
- TransDecoder (https://github.com/TransDecoder/TransDecoder)
- BUSCO (https://busco.ezlab.org/)
- FASTA/SSEARCH (https://fasta.bioch.virginia.edu/)
- CRISPRFlash (https://github.com/tsuryuchka/CRISPRFlash)
- Gene Ontology (http://geneontology.org/)

## Contact

For inquiries about tWhEEL, please contact the development team or open an issue on GitHub.
