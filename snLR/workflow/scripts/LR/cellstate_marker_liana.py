# %%
# Author: Robin Fallegger
# Date: 2025-July


# Imports
import os
import pandas as pd

import scanpy as sc
import liana as li

from itertools import product


def in_jupyter_notebook():
    try:
        get_ipython()
        return True
    except NameError:
        return False

if 'snakemake' in locals():
    wildcards = snakemake.wildcards
    # input files
    atlas_path = snakemake.input['atlas']
    marker_path = snakemake.input['markers']
    # lr_path = snakemake.input['lr']
    
    # output
    output_file = snakemake.output[0]
    
    threads = snakemake.threads
    expr_prop = snakemake.params.get('expr_prop', 0.1)
    cellstate_col = snakemake.params.get('cellstate_col', 'celltype_3')
    sources = snakemake.params.get('sources')
    targets = snakemake.params.get('targets')
    top_n = snakemake.params.get('top_n', 40)
    max_specificity = snakemake.params.get('max_specificity', 0.01)

elif in_jupyter_notebook():
    wildcards = {'omic': 'snRNA', 'celltype': 'PT', 'subset': 'xenium'}
    os.chdir('/gpfs/bwfor/work/ws/hd_mt265-work/SpaceKidneys')

    # input files
    atlas_path = 'sds/results/KPMPv2/{omic}/atlas_v2.h5ad'.format(omic=wildcards['omic'])
    marker_path = 'sds/data/KPMPv2/snRNA/markers/Subclass Level 2 within each Subclass Level 1/Human_Kidney_AtlasV2_{celltype}_l2_markers.txt'.format(celltype=wildcards['celltype'])
    lr_path = 'sds/resources/lr.csv'
    xenium_gl_path = 'sds/resources/Xenium_gene_list.csv'
    cytokine_path = 'sds/resources/uniprotkb_Cytokine_2025_05_30.tsv'
    threads = 1
    expr_prop = 0.1
    
    cellstate_col = 'celltype_3' if wildcards['omic'] == 'snRNA' else 'celltype_2'
    
    sources =  ['aPT2','aPT1','aPT-S1_S2','PT-S1','PT-S2','frPT-S1_S2','PT-S3','frPT-S3','cycPT']
    targets =  ["C-FIB", "C_M-FIB", "C-FIB-PATH", "C-FIB-OSMRhi", "C-FIB-OSMRlo", "C-MYOF",  "pvFIB", "pvMYOF", "B", "PL", "Naive Th", "MAIT", "CD8 TEM_TRM", "resMAC-LYVE1", "resMAC-HLAIIhi", "ncMON", "MON", "pDC", "moMAC-HBEGF", "moMAC-CXCL10"]

# %%
# load marker genes of the cell type
marker_genes = sorted(pd.read_csv(marker_path, sep='\t').query('p_val_adj < 0.05')['gene'].unique().tolist())

# %%
lrs = li.rs.select_resource('consensus')
# filter to ligand-receptor pairs that include marker genes
lrs = lrs.loc[lrs['ligand'].isin(marker_genes) | lrs['receptor'].isin(marker_genes),:]


#make a list of tuples from lrs
interactions_list = [tuple(x) for x in lrs.to_numpy()]
print('INFO: Found {} interactions in the resource'.format(len(interactions_list)))

# %%
adata = sc.read_h5ad(atlas_path)
adata = adata[adata.obs[cellstate_col].isin(sources) | adata.obs[cellstate_col].isin(targets)]

if wildcards['omic'] == 'scRNA':
    adata.obs.loc[adata.obs['celltype_2'] == "PL", 'celltype_1'] = "Lymphoid"

sc.pp.filter_genes(adata, min_cells=1)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# %%
pairs = list(product(sources, targets)) + list(product(targets, sources))
pairs = pd.DataFrame(pairs, columns=['source', 'target'])

# %%
li.mt.rank_aggregate(adata, 
                    groupby=cellstate_col, 
                    resource_name='consensus',
                    groupby_pairs=pairs,
                    expr_prop=expr_prop, 
                    use_raw = False,
                    n_jobs = threads,
                    inplace = True,
                    interactions =interactions_list,
                    verbose=True)


# %%
ontology = adata.obs[['celltype_1', cellstate_col]].drop_duplicates()
ontology['celltype_1'] = ontology['celltype_1'].astype('str')
ontology.loc[ontology['celltype_1'].isin(['Lymphoid', 'Myeloid']), 'celltype_1'] = 'Immune'

# %%
def get_top_cell_interactions(results, max_specificity=0.01, from_cs=None, CSOI=None, top_n=10):
    """
    Extract top cell-cell interactions from LIANA results based on specificity and magnitude ranks.
    
    Parameters:
    -----------
    adata : AnnData
        AnnData object containing LIANA results in adata.uns['liana_res']
    inter_pairs : pandas.DataFrame
        DataFrame containing a 'celltype_3' column for filtering source cells
    max_specificity : int, default 50
        Maximum specificity rank threshold for significant interactions
    from_cs : str, default None
        Cell source indicator. If 'PT', statistics are computed differently
    top_n : int, default 10
        Number of top interactions to return
        
    Returns:
    --------
    pandas.DataFrame
        Filtered interactions with magnitude ranks set to NA where 
        specificity_rank exceeds max_specificity
    """
    # results = results.loc[results['source'].isin(inter_pairs['celltype_3']),:]
    
    sig_inter = (
        results
        .query('specificity_rank <= @max_specificity')
        [['ligand_complex', 'receptor_complex']]
        .drop_duplicates()
        .merge(results, on=['ligand_complex', 'receptor_complex'], how='left')
    )

    if from_cs != CSOI:
        stats = sig_inter.groupby(['ligand_complex', 'receptor_complex', 'source']).agg({'magnitude_rank': ['mean', 'std']}).reset_index()
        stats.columns = ['ligand_complex', 'receptor_complex', 'source', 'magnitude_rank_mean', 'magnitude_rank_std']
        stats = stats.groupby(['ligand_complex', 'receptor_complex']).mean().reset_index()
    else:
        stats = sig_inter.groupby(['ligand_complex', 'receptor_complex', 'target']).agg({'magnitude_rank': ['mean', 'std']}).reset_index()
        stats.columns = ['ligand_complex', 'receptor_complex', 'target', 'magnitude_rank_mean', 'magnitude_rank_std']
        stats = stats.groupby(['ligand_complex', 'receptor_complex']).mean().reset_index()
    
    stats = stats.sort_values(by='magnitude_rank_std', ascending=False).head(top_n)
    stats['interaction'] = stats['ligand_complex'] + '_' + stats['receptor_complex']
    order = stats['interaction'].unique()
    
    results['interaction'] = results['ligand_complex'] + '_' + results['receptor_complex']
    results = results.loc[results['interaction'].isin(order), :].copy()
    results['interaction'] = pd.Categorical(results['interaction'], categories=order, ordered=True)
    results = results.sort_values(by='interaction')
    
    
    return results.drop(columns=['interaction'])

# %%
liana_res = (
    adata.uns['liana_res']
    .merge(ontology.rename(columns={'celltype_1': 'source_celltype', cellstate_col: 'source'}), on = 'source', how = 'left')
    .merge(ontology.rename(columns={'celltype_1': 'target_celltype', cellstate_col: 'target'}), on = 'target', how = 'left')
)
liana_res = liana_res.loc[:, liana_res.columns[-2:].tolist() + liana_res.columns[0:-2].tolist()]

# %%
sorted_results = []
for (src, trgt), LR_df in liana_res.groupby(['source_celltype', 'target_celltype'], observed=True):
    print(f'INFO: Processing interactions from {src} to {trgt}')
    sorted_results.append(get_top_cell_interactions(LR_df, max_specificity=max_specificity, from_cs=src, CSOI=wildcards['celltype'], top_n=top_n))
LR_df = pd.concat(sorted_results, ignore_index=True)

# %%
if not in_jupyter_notebook():
    LR_df.to_csv(output_file, index=False)