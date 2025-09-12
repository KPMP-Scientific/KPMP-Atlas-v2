# Single-Nucleus Ligand-Receptor Analysis Pipeline

This Snakemake pipeline reproduces the ligand-receptor (LR) analysis results presented in the KPMP Atlas v2 paper. The workflow performs cell-cell communication analysis on single-nucleus RNA-seq data to identify significant ligand-receptor interactions between different cell types in the kidney.

## Overview

The pipeline consists of two main modules:

1. **Data Conversion Module** (`workflow/rules/conversion.smk`)
   - Converts Seurat objects to H5 and H5AD formats
   - Extracts cell type hierarchy information
   - Prepares data for downstream ligand-receptor analysis

2. **Ligand-Receptor Analysis Module** (`workflow/rules/LR.smk`)
   - Performs cell-state-specific ligand-receptor analysis using LIANA
   - Analyzes specific ligand-receptor pairs (CCL2-CCR2, IL1B-IL1R1)
   - Generates publication-ready plots

## Paper Figures Reproduced

The pipeline generates the following figures from the KPMP Atlas v2 paper:

- **Figure 4i**: CCL2-CCR2 ligand-receptor interaction analysis
- **Figure 4k**: IL1B-IL1R1 ligand-receptor interaction analysis
- **Extended Data Figure 11**: Proximal tubule (PT) cell-state ligand-receptor analysis
- **Supplementary Figure 11**: Thick ascending limb (TAL) cell-state ligand-receptor analysis

## Prerequisites

### 1. Singularity/Apptainer Containers

The pipeline requires specific container images. Make sure to download and convert the Docker images to Singularity format before running the workflow:

```bash
# Navigate to the sings directory and follow the instructions
cd sings/
# See sings/README.md for detailed instructions on downloading:
# - scanpy.0.0.9.sif (Python environment with scanpy)
# - seurat.0.0.3.sif (R environment with Seurat)
```

Refer to [`sings/README.md`](sings/README.md) for detailed instructions on downloading and converting the required Docker images.

### 2. Snakemake Installation

Install Snakemake (version 6.0 or higher):

```bash
# Via conda/mamba (recommended)
conda install -c conda-forge snakemake

# Or via pip
pip install snakemake
```

### 3. Data Requirements

*[PLACEHOLDER: Data download rules will be added in future updates]*

The pipeline expects the following data structure:
```
smLR/
└── data/
    ├── KPMPv2/
    │   └── snRNA/
    │       ├── Kidney_AtlasV2_Seurat_05162024.rds
    │       ├── full_kidney_count_set_0424/
    │       ├── Supplementary Tables Version 1 05-2024-3.xlsx
    │       └── markers/
    │           └── Subclass Level 2 within each Subclass Level 1/
    │               ├── Human_Kidney_AtlasV2_PT_l2_markers.txt
    │               └── Human_Kidney_AtlasV2_TAL_l2_markers.txt
    └── sample_metadata.csv
```

## Usage

### Run the Complete Pipeline

Execute the entire ligand-receptor analysis pipeline:

```bash
# Run all rules to generate paper figures
snakemake --use-singularity --cores 8

# Or run with specific resource constraints
snakemake --use-singularity --cores 8 --resources mem_mb=250000
```

## Configuration

The pipeline behavior is controlled by `config.yaml`, which defines:

- **Cell type mappings**: Source and target cell types for ligand-receptor analysis
- **PT analysis**: Proximal tubule cell states and their target cell types (based on Visium niches)
- **TAL analysis**: Thick ascending limb cell states and their target cell types (based on Visium niches)

## Computational Requirements

- **Memory**: Up to 230 GB RAM for data conversion steps
- **Runtime**: Approximately 1-2 hours for complete pipeline
- **Storage**: Several GB for intermediate and final results
- **CPU**: Multi-threading supported (up to 8 cores recommended)