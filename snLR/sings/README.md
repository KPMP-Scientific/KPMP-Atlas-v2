# Singularity Container Images

This directory contains Singularity (`.sif`) container images used by the Snakemake workflow for single-nucleus ligand-receptor analysis.

## Required Container Images

The workflow requires two container images:
- `scanpy.0.0.9.sif` - Python environment with scanpy for single-cell analysis
- `seurat.0.0.3.sif` - R environment to convert Seurat objects to H5AD format

## Prerequisites

Make sure you have Apptainer (formerly Singularity) installed on your system:

```bash
# Check if apptainer is installed
apptainer --version
```

If not installed, you only require the unprivileged installation to download and convert Docker images to Singularity format. You can find installation instructions [here](https://apptainer.org/docs/admin/main/installation.html#install-unprivileged-from-pre-built-binaries).

## Downloading and Converting Docker Images

### Method 1: Direct conversion from Docker Hub (Recommended)

Convert the Docker images directly to Singularity format:

```bash
# Make sure you are in the snLR/sings directory
cd snLR/sings

# Convert scanpy Docker image to Singularity
apptainer build scanpy.0.0.9.sif docker://falrob/scanpy:0.0.9

# Convert seurat Docker image to Singularity
apptainer build seurat.0.0.3.sif docker://falrob/seurat:0.0.3
```

### Method 2: Using Docker and then converting

If you prefer to pull with Docker first:

```bash
# Pull Docker images
docker pull falrob/scanpy:0.0.9
docker pull falrob/seurat:0.0.3

# Convert to Singularity format
apptainer build scanpy.0.0.9.sif docker-daemon://falrob/scanpy:0.0.9
apptainer build seurat.0.0.3.sif docker-daemon://falrob/seurat:0.0.3
```

## Usage in Snakemake

These container images are automatically used by the Snakemake workflow rules defined in:
- `workflow/rules/conversion.smk`
- `workflow/rules/LR.smk`

Ensure the paths in the `singularity` directive of each rule match the location of your `.sif` files.

## Container Contents

### scanpy.0.0.9.sif
- Python environment with scanpy for single-cell RNA-seq analysis
- Additional packages for ligand-receptor analysis
- Required for Python-based analysis scripts

### seurat.0.0.3.sif
- R environment with Seurat for single-cell RNA-seq analysis
- Required for converting to H5AD format

## Notes

- Container images should be placed in the `sings/` directory relative to the Snakemake workflow
- The workflow expects specific container names as defined in the `.smk` files