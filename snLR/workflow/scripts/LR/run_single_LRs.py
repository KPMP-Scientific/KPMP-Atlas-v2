# %%
# Author: Robin Fallegger
# Date: 2025-July


# Imports

import scanpy as sc
import liana as li


def in_jupyter_notebook():
    try:
        get_ipython()
        return True
    except NameError:
        return False

if 'snakemake' in locals():
    wildcards = snakemake.wildcards
    # input files
    input_path = snakemake.input[0]
    # output
    output_file = snakemake.output[0]
    
    ligand = wildcards['ligand']
    receptor = wildcards['receptor']
    
    celltypes = snakemake.params.get('celltypes', ['PT', 'TAL', 'Lymphoid', 'Myeloid', 'FIB', 'Interstitial'])
    cellstate_col = snakemake.params.get('cellstate_col', 'celltype_3')
    
    threads = snakemake.threads
    expr_prop = snakemake.params.get('expr_prop', 0.05)

else:
    raise ValueError('Not running in snakemake, or jupyter')

print('INFO: Using expr_prop:', expr_prop)

# %%
adata = sc.read_h5ad(input_path)
if wildcards['ligand'] == 'CCL2':
    fib_cs = ['C-FIB-OSMRhi', 'C-FIB-OSMRlo']
    adata = adata[adata.obs['celltype_1'].isin(celltypes) | adata.obs[cellstate_col].isin(fib_cs)]
else:
    adata = adata[adata.obs['celltype_1'].isin(celltypes),:]

if wildcards['omic'] == 'snRNA':
    adata.obs['celltype_3'] = adata.obs['celltype_3'].astype(str)
    adata.obs.loc[adata.obs['celltype_3'].isin(['PT-S1', 'PT-S2', 'PT-S3']), 'celltype_3'] = 'PT'
    adata.obs.loc[adata.obs['celltype_3'].isin(['C_M-TAL-A','C-TAL-A','C_M-TAL-B','C-TAL-B', 'M-TAL']), 'celltype_3'] = 'TAL'
    
adata = adata[~adata.obs[cellstate_col].str.startswith('d'),:]
sc.pp.filter_genes(adata, min_cells = 3)
sc.pp.normalize_total(adata, target_sum = 1e5)
sc.pp.log1p(adata)

# %%
interactions_list = [(ligand, receptor)]

# %%
li.mt.rank_aggregate(adata,
                    groupby=cellstate_col,
                    resource_name='consensus',
                    expr_prop=expr_prop,
                    use_raw = False,
                    n_jobs = threads,
                    inplace = True,
                    interactions =interactions_list,
                    verbose=True)

# %%
if not in_jupyter_notebook():
    adata.uns['liana_res'].to_csv(output_file, index=False)