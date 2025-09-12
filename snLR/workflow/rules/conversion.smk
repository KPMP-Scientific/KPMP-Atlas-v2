##################################################
# Conversion of the single nuclei atlas
##################################################

rule sn_to_h5:
    input:
        data = 'data/KPMPv2/snRNA/Kidney_AtlasV2_Seurat_05162024.rds',
        data_dir = 'data/KPMPv2/snRNA/full_kidney_count_set_0424'
    output:
        temp('results/KPMPv2/snRNA/atlas_v2.h5')
    params:
        assays = ['RNA'],
        reductions = {'umap': 'umap', 'pca': 'pca'}
    resources:
        mem_mb = 30000
    singularity:
        'sings/seurat.0.0.3.sif'
    script:
        '../scripts/conversion/sn_seurat_to_h5.R'

rule sn_to_h5ad:
    input:
        h5 = 'results/KPMPv2/snRNA/atlas_v2.h5',
        meta = 'data/KPMPv2/snRNA/Supplementary Tables Version 1 05-2024-3.xlsx'
    output:
        adata = 'results/KPMPv2/snRNA/atlas_v2.h5ad'
    params:
        assays = ['RNA'],
        reductions = {'umap': 'umap', 'pca': 'pca'}
    resources:
        mem_mb = 230000,
        runtime = 30
    singularity:
        'sings/scanpy.0.0.9.sif'
    script:
        '../scripts/conversion/sn_h5_to_h5ad.py'

rule get_celltype_hierarchy:
    input:
        adata = 'results/KPMPv2/snRNA/atlas_v2.h5ad'
    output:
        hierarchy = 'results/KPMPv2/snRNA/celltype_hierarchy.csv'
    resources:
        mem_mb = 2000,
        runtime = 5 #minutes
    singularity:
        'sings/scanpy.0.0.9.sif'
    script:
        '../scripts/conversion/get_celltype_hierarchy.py'