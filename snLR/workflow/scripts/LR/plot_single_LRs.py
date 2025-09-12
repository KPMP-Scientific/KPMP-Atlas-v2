# %%
# Author: Robin Fallegger
# Date: 2025-July


# Imports
import pandas as pd
import numpy as np

import scanpy as sc

import marsilea as ma
import marsilea.plotter as mp

import matplotlib.pyplot as plt

from matplotlib.backends.backend_pdf import PdfPages


if 'snakemake' in locals():
    wildcards = snakemake.wildcards
    ligand = wildcards['ligand']
    receptor = wildcards['receptor']
    
    # input files
    input_path = snakemake.input[0]
    lr_path = snakemake.input[1]
    
    # output
    output_file = snakemake.output[0]
    
    cellstate_col = snakemake.params.get('cellstate_col', 'celltype_3')
    significance_threshold = snakemake.params.get('significance_threshold', 0.1)
else:
    raise ValueError('Not running in snakemake, or jupyter')


# %%
# load results
adata = sc.read_h5ad(input_path, backed = 'r')
liana_res = pd.read_csv(lr_path, sep = ',')

if wildcards['omic'] == 'snRNA':
    adata.obs['celltype_3'] = adata.obs['celltype_3'].astype(str)
    adata.obs.loc[adata.obs['celltype_3'].isin(['PT-S1', 'PT-S2', 'PT-S3']), 'celltype_3'] = 'PT'

    adata.obs.loc[adata.obs['celltype_3'].isin(['C_M-TAL-A','C-TAL-A','C_M-TAL-B','C-TAL-B', 'M-TAL']), 'celltype_3'] = 'TAL'


# %%
# define color maps
vmap = plt.cm.viridis
vmap.set_bad('lightgrey', 1.0)

celltype_map = {
    'PT': '#1f77b4',
    'TAL': '#9467bd',
    'Lymphoid': '#ff7f0e',
    'Myeloid': '#2ca02c'
}
if wildcards['omic'] == 'scRNA':
    celltype_map['Interstitial'] = '#d62728'
elif wildcards['omic'] == 'snRNA':
    celltype_map['FIB'] = '#d62728'
    
color_max = (-np.log10(liana_res['magnitude_rank'])).quantile(0.975)
size_max = (-np.log10(liana_res['specificity_rank'])).quantile(0.975)
    
# prepare row and column annotations
celltypes = adata.obs[['celltype_1', cellstate_col]].drop_duplicates().sort_values(by = 'celltype_1')


# %%
# data transformation
mat_long = liana_res.query('receptor_complex == @receptor and ligand_complex == @ligand')

# wide format and log10 transformation
mat_mag = mat_long.pivot(index = 'source', columns = 'target', values = 'magnitude_rank')
mat_mag = -np.log10(mat_mag)
mat_spec = mat_long.pivot(index = 'source', columns = 'target', values = 'specificity_rank')
mat_spec = -np.log10(mat_spec)

# filter for significance
if wildcards['greyed'] == 'greyed':
    mat_mag[mat_spec < -np.log10(significance_threshold)] = np.nan


mat_mag = mat_mag.loc[celltypes.loc[celltypes[cellstate_col].isin(mat_mag.index), cellstate_col], celltypes.loc[celltypes[cellstate_col].isin(mat_mag.columns), cellstate_col]]
mat_spec = mat_spec.loc[celltypes.loc[celltypes[cellstate_col].isin(mat_spec.index), cellstate_col], celltypes.loc[celltypes[cellstate_col].isin(mat_spec.columns), cellstate_col]]

# keep only columns and rows with at least one non-NaN value and ensure PT and TAL are always included

if wildcards['greyed'] == 'greyed':
    keep_targets = (mat_mag.notna()).sum(axis = 0) > 0
    keep_sources = (mat_mag.notna()).sum(axis = 1) > 0
else:
    keep_targets = (mat_spec > -np.log10(significance_threshold)).sum(axis = 0) > 0
    keep_sources = (mat_spec > -np.log10(significance_threshold)).sum(axis = 1) > 0

keep_targets = keep_targets[keep_targets].index.tolist() + [col for col in mat_mag.columns if 'PT' in col or 'TAL' in col]
keep_sources = keep_sources[keep_sources].index.tolist() + [row for row in mat_mag.index if 'PT' in row or 'TAL' in row]

if ligand == "CCL2" and receptor == "CCR2":
    # add PT and TAL to targets
    if 'PT' not in keep_sources:
        keep_sources.append('PT')

# order rows and columns by celltype_1 and cellstate_col
target_annot = celltypes.loc[celltypes[cellstate_col].isin(keep_targets), :].sort_values(by=['celltype_1', cellstate_col])
source_annot = celltypes.loc[celltypes[cellstate_col].isin(keep_sources), :].sort_values(by=['celltype_1', cellstate_col])

# reorder matrices
mat_mag = mat_mag.reindex(source_annot[cellstate_col], axis=0).reindex(target_annot[cellstate_col], axis=1)
mat_spec = mat_spec.reindex(source_annot[cellstate_col], axis=0, fill_value=0).reindex(target_annot[cellstate_col], axis=1, fill_value=0)

#clip mat_spec at size_max
mat_spec = mat_spec.clip(upper = size_max)


# %%
# plot heatmap
def make_heatmap(mat_spec, mat_mag, celltype_map, color_max, ligand, receptor, row_annot, column_annot, row_title, col_title):
    """
    Create a heatmap with the given matrices and parameters.
    """
    h = ma.SizedHeatmap(mat_spec, 
                        color=mat_mag, 
                        cmap =  vmap, 
                        vmax = color_max, 
                        color_legend_kws = {'title': 'Magnitude\nrank'} , 
                        size_legend_kws = {'title': 'Specificity\nrank'}, 
                        plotnonfinite = True, 
                        width = 0.15 * mat_spec.shape[1] + 1, 
                        height = 0.175 * mat_spec.shape[0])

    h.add_left(mp.Colors(row_annot['celltype_1'], palette = celltype_map, label_props={'rotation': 0}), size = 0.1, pad = 0.05, name =row_title)
    h.add_left(mp.Labels(mat_spec.index, fontsize=9,label=row_title, label_loc='left'))

    if col_title == "Source":
        h.add_top(mp.Colors(column_annot['celltype_1'], palette = celltype_map, label_props={'rotation': 0}), size = 0.1, pad = 0.05, name =col_title)
        h.add_top(mp.Labels(mat_spec.columns, fontsize=9, label=col_title, label_loc='top', label_props={'rotation': 0}))
    elif col_title == "Target":
        h.add_bottom(mp.Colors(column_annot['celltype_1'], palette = celltype_map, label_props={'rotation': 0}), size = 0.1, pad = 0.05, name =col_title)
        h.add_bottom(mp.Labels(mat_spec.columns, fontsize=9, label=col_title, label_loc='bottom', label_props={'rotation': 0}))

    h.add_title(ligand + '-' + receptor, pad = 0.5, fontsize=12)
    h.add_legends()
    h.render()
    plt.close()

    fig = h.figure
    fig.set_dpi(150)
    
    return fig

st_fig = make_heatmap(mat_spec, mat_mag, celltype_map, color_max, ligand, receptor, source_annot, target_annot, 'Source', 'Target')
ts_fig = make_heatmap(mat_spec.T, mat_mag.T, celltype_map, color_max, ligand, receptor, target_annot, source_annot, 'Target', 'Source')

# %%
if not in_jupyter_notebook():
    with PdfPages(output_file) as pdf:
        pdf.savefig(st_fig, bbox_inches='tight')
        pdf.savefig(ts_fig, bbox_inches='tight')
