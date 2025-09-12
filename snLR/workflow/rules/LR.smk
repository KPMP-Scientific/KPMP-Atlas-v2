
rule cellstate_marker_liana:
    input:
        atlas = 'results/KPMPv2/{omic}/atlas_v2.h5ad',
        markers = 'data/KPMPv2/snRNA/markers/Subclass Level 2 within each Subclass Level 1/Human_Kidney_AtlasV2_{celltype}_l2_markers.txt',
    output:
        'results/LR/{omic}/{celltype}_liana.csv'
    params:
        sources = lambda w: config['LR'][w.omic][w.celltype]['sources'],
        targets = lambda w: config['LR'][w.omic][w.celltype]['targets'],
        cellstate_col = 'celltype_3',
        expr_prop = 0.1,
        top_n = 40,
        max_specificity = 0.01
    resources:
        mem_mb = 120000,
        runtime = 20
    threads: 8
    singularity:
        'sings/scanpy.0.0.9.sif'
    script:
        '../scripts/LR/cellstate_marker_liana.py'


rule plot_cellstate_marker_liana:
    input:
        liana_results = 'results/LR/{omic}/{celltype}_liana.csv'
    output:
        'plots/LR/{omic}/{celltype}_liana_plot.pdf'
    params:
        significance_threshold = 0.01
    resources:
        mem_mb = 5000,
        runtime = 5
    threads: 1
    singularity:
        'sings/scanpy.0.0.9.sif'
    script:
        '../scripts/LR/plot_cellstate_marker_liana.py'


rule single_lr_analysis:
    input:
        atlas = 'results/KPMPv2/{omic}/atlas_v2.h5ad'
    output:
        'results/LR/{omic}/{ligand}_{receptor}_single_lr.csv'
    params:
        celltypes = lambda w: ['PT', 'TAL', 'Lymphoid', 'Myeloid', 'FIB', 'Interstitial'] if w.ligand == 'IL1B' else ['PT', 'Lymphoid', 'Myeloid', 'FIB'],
        cellstate_col = 'celltype_3',
        expr_prop = 0.05
    resources:
        mem_mb = 150000,
        runtime = 15
    threads: 4
    singularity:
        'sings/scanpy.0.0.9.sif'
    script:
        '../scripts/LR/run_single_LRs.py'


rule plot_single_lr:
    input:
        atlas = 'results/KPMPv2/{omic}/atlas_v2.h5ad',
        lr_results = 'results/LR/{omic}/{ligand}_{receptor}_single_lr.csv'
    output:
        'plots/LR/{omic}/{ligand}_{receptor}_single_lr_{greyed}.pdf'
    params:
        cellstate_col = 'celltype_3',
        significance_threshold = 0.1
    resources:
        mem_mb = 15000,
        runtime = 5
    threads: 1
    singularity:
        'sings/scanpy.0.0.9.sif'
    script:
        '../scripts/LR/plot_single_LRs.py'