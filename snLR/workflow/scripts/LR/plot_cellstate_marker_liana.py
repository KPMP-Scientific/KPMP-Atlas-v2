# %%
# Author: Robin Fallegger
# Date: 2025-July


# Imports
import pandas as pd
import numpy as np

import marsilea as ma
import marsilea.plotter as mp
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


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
    
    significance_threshold = snakemake.params.get('significance_threshold', 0.01)
else:
    raise ValueError('Not running in snakemake, or jupyter')


# %%
LRs = pd.read_csv(input_path, sep = ',')

# %%
plots = []
color_max = (-np.log10(LRs['magnitude_rank'])).quantile(0.975)
size_max = (-np.log10(LRs['specificity_rank'])).quantile(0.975)

for (src, trgt), LR_df in LRs.groupby(['source_celltype', 'target_celltype'], observed=True):
    print(f'INFO: Processing interactions from {src} to {trgt}')
    
    LR_df['interaction'] = LR_df.apply(lambda x: f"{x['ligand_complex']} -> {x['receptor_complex']}", axis=1)
    order = LR_df[['ligand_complex', 'receptor_complex']].drop_duplicates()
    order['interaction'] = order.apply(lambda x: f"{x['ligand_complex']} -> {x['receptor_complex']}", axis=1)
    
    mat_mag = LR_df.pivot(index = ['interaction'], columns = ['source', 'target'], values= 'magnitude_rank').loc[order['interaction'], :]
    mat_mag.sort_values(by=list(mat_mag.columns.names), axis=1, inplace=True)
    mat_mag = np.log10(mat_mag).abs()
    mat_spec = LR_df.pivot(index = ['interaction'], columns = ['source', 'target'], values= 'specificity_rank').loc[order['interaction'], :]
    mat_spec = np.log10(mat_spec).abs()
    mat_spec = mat_spec.loc[:, mat_mag.columns]

    # filter for significance
    # mat_mag[mat_spec < -np.log10(significance_threshold)] = np.NAN
    mat_spec.fillna(0, inplace=True)

    mat_spec = mat_spec.clip(upper = size_max)

    col_annot = mat_spec.columns.to_frame()

    vmap = plt.cm.viridis
    vmap.set_bad('lightgrey', 1.0)

    h = ma.SizedHeatmap(mat_spec, 
                        color=mat_mag, 
                        cmap =  vmap, 
                        vmax = color_max, 
                        color_legend_kws = {'title': 'Magnitude rank'} , 
                        size_legend_kws = {'title': 'Specificity rank'}, 
                        plotnonfinite = True, 
                        width = 0.175 * mat_spec.shape[1], 
                        height = 0.165 * mat_spec.shape[0])
    h.add_left(mp.Labels(mat_spec.index, fontsize=8, rotation=0))

    h.add_bottom(mp.Labels(col_annot['target'], fontsize=8, rotation=90))

    h.group_cols(mat_spec.columns.to_frame()['source'],spacing= .0005)
    h.add_top(mp.Chunk(col_annot['source'].unique(), bordercolor="gray"), pad=0.1)
    h.add_dendrogram('top')

    h.add_title('Intercellular interactions from {0} to {1} cells ({2})'.format(src, trgt, wildcards['omic']), fontsize=12)

    h.add_legends()

    h.render()
    plt.close()

    fig = h.figure
    fig.set_dpi(150)
    plots.append(fig)

# %%
if not in_jupyter_notebook():
    with PdfPages(output_file) as pdf:
        for fig in plots:
            pdf.savefig(fig, bbox_inches='tight')