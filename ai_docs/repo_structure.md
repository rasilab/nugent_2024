# Repository Structure and Conventions Guide

## Overview

This repository contains the complete computational analysis pipeline for the Nature Methods publication "Decoding post-transcriptional regulatory networks by RNA-linked CRISPR screening in human cells". The codebase follows best practices for computational reproducibility in bioinformatics research.

## Repository Architecture

### Top-Level Organization

```
nugent_2024/
├── README.md                    # Main documentation with figure-to-code mapping
├── run_everything.sh            # Master orchestration script
├── analysis/                    # All analysis pipelines and scripts
├── data/                        # Raw experimental data
├── source_data/                 # Processed data for manuscript figures
├── specs/                       # Technical specifications and documentation
├── ai_docs/                     # AI-generated documentation
└── Nature_Methods_manuscript.pdf
```

### Analysis Directory Structure

The `analysis/` directory is organized by experimental technique:

```
analysis/
├── submit_cluster.sh           # SLURM cluster submission wrapper
├── submit_local.sh             # Local execution wrapper
├── cluster.yaml                # SLURM resource specifications
├── barcodeseq/                 # CRISPR barcode sequencing analysis
├── flow_cytometry/             # Flow cytometry (FACS) analysis
├── polysome_profiling/         # Polysome profiling experiments
├── qpcr/                       # qPCR validation experiments
├── riboseq/                    # Ribosome profiling (Ribo-seq)
└── rnaseq/                     # RNA sequencing analysis
```

Each subdirectory follows a consistent structure:
```
<technique>/
├── annotations/               # Sample metadata and annotations
├── data/                      # Processed data (intermediate files)
├── figures/                   # Generated figure PDFs
├── scripts/                   # Analysis scripts (R, Python, Snakemake)
└── logs/                      # Workflow execution logs
```

## Technology Stack

### Core Languages
- **R (4.x)**: Primary analysis language for statistical analysis and plotting
- **Python (3.x)**: Data preprocessing, bioinformatics pipelines
- **Bash/AWK**: High-performance text processing and workflow glue

### Workflow Management
- **Snakemake**: Reproducible workflow orchestration
- **SLURM**: High-performance computing cluster job management
- **GNU Make**: Simple task automation

### Containerization
- **Singularity/Apptainer**: HPC-compatible container runtime
- **Docker**: Container images hosted on GitHub Container Registry
- **Container namespace**: `ghcr.io/rasilab/*`

### Key R Libraries
- **tidyverse**: Data manipulation and visualization
- **rasilabRtemplates**: Custom lab plotting themes and utilities
- **ggplot2**: Publication-quality figures
- **DESeq2**: Differential expression analysis
- **flowCore**: Flow cytometry data processing

### Bioinformatics Tools (Containerized)
- **MaGeCK**: CRISPR screen analysis
- **STAR**: RNA-seq spliced alignment
- **bowtie**: Short read alignment
- **cutadapt**: Adapter trimming
- **samtools**: SAM/BAM file manipulation
- **parallel-fastq-dump**: SRA data download
- **entrez-direct**: NCBI database queries

## Coding Conventions

### R Script Conventions

#### Naming Conventions
- **Variables**: snake_case (e.g., `sample_annotations`, `plot_data`)
- **Functions**: snake_case for custom functions
- **Files**: lowercase with underscores (e.g., `plot_grna_fitness_results.R`)
- **Constants**: UPPERCASE or prefixed (e.g., `cbPalette`)

#### Code Structure Pattern
```r
# 1. Load libraries
library(tidyverse)
library(rasilabRtemplates)

# 2. Set theme
theme_set(theme_rasilab())

# 3. Define color palettes
cbPalette <- c("#999999", "#E69F00", ...)

# 4. Read data
data <- read_csv("file.csv", show_col_types = F)

# 5. Process data
processed_data <- data %>%
  filter(...) %>%
  mutate(...) %>%
  group_by(...) %>%
  summarize(...)

# 6. Export source data
write_csv(processed_data, "../../source_data/figure_X.csv")

# 7. Create plot
p <- ggplot(processed_data, aes(...)) +
  geom_point() +
  theme(...)

# 8. Save figure
ggsave("../figures/figure_name.pdf", p, width = 4, height = 3)
```

#### Tidyverse Patterns
- Extensive use of pipe operator `%>%`
- Prefer `dplyr` verbs over base R
- Use `purrr::map()` for iteration
- Chain operations with intermediate `print()` for debugging

### Python Script Conventions

#### Style
- PEP 8 compliance
- Type hints where appropriate
- Comprehensive docstrings
- Command-line interfaces with `argparse`

#### Bioinformatics Patterns
- GTF/GFF parsing with custom functions
- NumPy/Pandas for numerical operations
- Biopython for sequence manipulation

### Snakemake Workflow Conventions

#### Rule Structure
```python
rule rule_name:
    input:
        file1 = "path/to/input1",
        file2 = "path/to/input2"
    output:
        "path/to/output"
    params:
        param1 = "value"
    threads: 4
    container: "docker://ghcr.io/rasilab/tool:version"
    shell:
        """
        command {input.file1} {params.param1} > {output}
        """
```

#### Best Practices
- Explicit container versions
- Thread specifications for parallelization
- Log files for debugging
- Functions for dynamic input generation
- Checkpoints for variable outputs

## Data Conventions

### Sample Naming
- **Sample IDs**: Numeric patterns (e.g., `219p138`, `225p41`)
- **SRR accessions**: For public data (e.g., `SRR12345678`)
- **Illumina IDs**: Separate from analysis IDs

### File Formats
- **Tabular data**: CSV format exclusively
- **Figures**: PDF format for publication quality
- **Sequences**: FASTQ, FASTA
- **Alignments**: BAM (sorted and indexed)
- **Annotations**: GTF/GFF3

### Directory Conventions
- **Raw data**: Original, unmodified files in `data/`
- **Processed data**: Intermediate files in `analysis/*/data/`
- **Annotations**: Metadata in `analysis/*/annotations/`
- **Figures**: PDFs in `analysis/*/figures/`
- **Source data**: Figure data in `source_data/`

## Workflow Patterns

### Master Orchestration
The `run_everything.sh` script:
1. Downloads data from SRA
2. Runs preprocessing pipelines
3. Executes analysis scripts
4. Generates all figures
5. Exports source data

### Individual Workflows
Each analysis typically follows:
1. **Snakemake pipeline**: Data preprocessing
2. **R scripts**: Statistical analysis and plotting
3. **Source data export**: Before figure generation
4. **Figure generation**: Consistent formatting

### Reproducibility Features
- All software containerized with version pins
- Snakemake ensures complete rerun capability
- README maps every figure to its code
- Source data provided for all plots
- Raw data accessible via SRA

## Quality Assurance

### Code Quality
- Suppress warnings appropriately: `options(warn = -1)`
- Clean package loading: `suppressPackageStartupMessages()`
- Explicit NA handling
- Error handling in workflows

### Documentation
- Section headers in scripts
- Figure references in filenames
- README with comprehensive mapping
- Clean commit messages

### Testing
- Visual inspection of figures
- Source data validation
- Workflow dry-runs with Snakemake

## Best Practices Summary

1. **Reproducibility First**
   - Everything in containers
   - Version control for all code
   - Clear data-to-figure pipeline

2. **Consistent Organization**
   - Standard directory structure
   - Predictable file naming
   - Uniform coding style

3. **Clear Documentation**
   - README as central hub
   - Comments explain "why" not "what"
   - Figure-code mapping maintained

4. **Efficient Workflows**
   - Snakemake for automation
   - Parallel execution where possible
   - Modular script design

5. **Publication Ready**
   - PDF figures for print quality
   - Source data for transparency
   - Clean, professional output

## Common Operations

### Running a Specific Analysis
```bash
cd analysis/<technique>/scripts
bash ../submit_local.sh <workflow>.smk
```

### Generating a Figure
```bash
cd analysis/<technique>/scripts
Rscript <plot_script>.R
```

### Adding a New Analysis
1. Create directory structure under `analysis/`
2. Add Snakemake workflow if needed
3. Write R script following conventions
4. Update `run_everything.sh`
5. Add figure mapping to README

## Container Usage

### Local Execution
Containers are automatically pulled and cached by Singularity when running Snakemake workflows.

### Cluster Execution
Mount the lab's shared container directory to avoid redundant downloads:
```bash
export SINGULARITY_CACHEDIR=/shared/containers
```

## Troubleshooting

### Common Issues
1. **Missing containers**: Check internet connection and registry access
2. **Path errors**: Ensure proper bind mounts in submission scripts
3. **Memory errors**: Adjust resource requests in `cluster.yaml`
4. **Missing dependencies**: Verify container versions

### Debug Strategies
- Run Snakemake with `-n` for dry-run
- Check log files in `logs/` directories
- Test individual rules with `--until`
- Verify data paths are absolute where required

## Future Maintenance

### Adding New Data
1. Place raw data in appropriate `data/` subdirectory
2. Update `sample_annotations.csv`
3. Modify Snakemake workflow if needed
4. Run analysis pipeline

### Updating Analyses
1. Test changes locally first
2. Update container versions if needed
3. Verify figure output matches expected
4. Update source data files
5. Commit with clear message

This repository represents a gold standard for reproducible computational biology research, with clear organization, comprehensive documentation, and robust automation.