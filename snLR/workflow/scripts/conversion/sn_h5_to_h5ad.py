# %%
from scipy.sparse import csr_matrix
import h5py
import pandas as pd

import platform
import os

import scanpy as sc
import muon as mu
import numpy as np

# %%
def readh5_to_anndata(file, main_assay = 'RNA'):
    """ Reads from h5 file and returns anndata object

    Currently only support a single type of data (i.e. data OR counts OR scaled data from Seurat).
    Args:
        input_file (path): Path to h5 file
        main_assay (string): Name of assay in h5 file to use as X

    Returns:
        Anndata: Anndata with main_assay as X and all other assays as layers.
    """
    with h5py.File(file, "r") as f:

        # check that main assay is in f
        # if main_assay not in list(f['assays'].keys()):
        #     raise ValueError('Main assay %s not in h5 file' % main_assay)
        
        # check that obs is in f
        if 'obs' not in list(f.keys()):
            raise ValueError('Obs not in h5 file')

        # create anndata object
        print('Loading main assay %s' % main_assay)
        X = f['assays'][main_assay]

        mat = csr_matrix((X['data'], X['indices'], X['indptr']), shape=X.attrs['shape'])
        adata = sc.AnnData(X=mat, dtype=np.float64)

        print('INFO: first 20 cells and genes in the X matrix')
        print(adata.X[0:20, 0:20])

        adata.obs_names = f['obs']['index'][()].astype(str)
        adata.var_names = f['var']['index'][()].astype(str)

        # add obs
        print('Loading obs')
        obs = f['obs']
        for key in obs.keys():
            print('INFO: Loading obs: %s' % key)
            # Check if key is a factor and if so convert
            if key + '_levels' in obs.keys():
                idx = obs[key][()]
                adata.obs[key] = obs[key + '_levels'][()][idx].astype(str)#.str.decode("utf-8")
            elif '_levels' not in key and 'colnames' not in key:
                # if bytes in [x for x in obs[key][()].apply(type).unique()]:]
                #     adata.obs[key] = obs[key][()].astype(str)#.str.decode("utf-8")
                # else:
                adata.obs[key] = obs[key][()]
            print('INFO: done loading obs %s' % key)

        print('INFO: head of obs')
        print(adata.obs.head())

        # # load remaining assays as layers
        # layers = [key for key in f['assays'].keys() if 'RNA' not in key]
        # for layer in layers:
        #     print('Loading layer: %s' % layer)
        #     X = f['assays'][layer]

        #     # add only if dims match
        #     if X['dims'] == adata.shape:
        #         adata.layer[layer] = csr_matrix((X['values'], X['indices'], X['indptr']), shape=(X['dims'][0], X['dims'][1]))
        #     else:
        #         raise Warning('Layer %s not added because dims do not match' % layer)
        #     print('INFO: first 20 cells and genes in the %s layer' % layer)
        #     print(adata.layer[layer][0:20, 0:20])

        # load reductions
        if 'reductions' in f.keys():
            print('Loading reductions')
            reductions = f['reductions']
            for red in reductions.keys():
                print('\tLoading reduction: %s' % red)
                embedding = reductions[red]

                adata.obsm['X_' + red] = np.array([embedding[key][()] for key in embedding.keys()], dtype=np.float64).T

    return adata

# %%

if 'snakemake' in locals():
    input_file = snakemake.input['h5']
    meta_file = snakemake.input['meta']
    wildcards = snakemake.wildcards
    
    output_file = snakemake.output[0]

    main_assay = snakemake.params['assays'][0]

else:
    os.chdir('/gpfs/bwfor/work/ws/hd_mt265-work/SpaceKidneys')
    wildcards = {'omic': 'scRNA', 'source': 'GEO'}
    input_file = 'results/KPMPv2/atlas_v2.h5'
    input_file = os.path.join('sds', input_file)
    
    meta_file = 'sds/data/KPMPv2/snRNA/Supplementary Tables Version 1 05-2024-3.xlsx'

    main_assay = 'RNA'

# %%
adata = readh5_to_anndata(input_file, main_assay = main_assay)
# adata = sc.read_h5ad(input_file+'ad', backed = 'r')
adata.obs = adata.obs.rename(columns={'v2.clusters': 'clusters'}).drop(columns=['index', 'n_genes'], errors = 'ignore')

# %%
if 'bpcells_name' in adata.obs.columns:
    adata.obs.drop(columns='bpcells_name', inplace=True)
del adata.obsm

print(adata)

# %%
print('Adding patient metadata')
patient_metadata = pd.read_excel(meta_file,
              sheet_name='Table S1',
              skiprows = 4).iloc[:,1:]
patient_metadata['patient'] = patient_metadata['patient'].astype(str)

# remove trailing whitespace in all character columns
patient_metadata = patient_metadata.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
metadata = adata.obs.reset_index().merge(patient_metadata, on = ['library', 'patient', 'specimen'], how = 'left').rename(columns={'index': 'cell_id'})

# %%
print('Adding cell type annotations')
cluster_metadata = pd.read_excel(meta_file,
                                 sheet_name='Table S - sn annot',
                                    skiprows = 4).iloc[:,1:]

# look for string with 'ï' or + signs in the character columns
cluster_metadata = cluster_metadata.applymap(lambda x: x.replace('ï','i') if isinstance(x, str) else x)
cluster_metadata = cluster_metadata.applymap(lambda x: x.replace('+','') if isinstance(x, str) else x)
cluster_metadata = cluster_metadata.applymap(lambda x: x.replace('/','_') if isinstance(x, str) else x)

#replace all . in the column names with _
cluster_metadata.columns = cluster_metadata.columns.str.replace('.','_')
cluster_metadata = cluster_metadata.filter(regex='^v2_|integration_stable', axis = 1)

# %%
metadata = metadata.merge(cluster_metadata, left_on = 'clusters', right_on='v2_clusters', how = 'left')

# %%
column_mapping = {
    #ids
    'cell_id': 'cell_id',
    'library': 'uid_library',
    'patient': 'uid_donor_id',
    'specimen': 'uid_specimen',
    #conditions
    'condition_level3': 'condition_3',
    'condition_level2': 'condition_2',
    'condition_level1': 'condition_1',
    'condition': 'condition_full',
    #celltypes
    'v2_clusters': 'clusters',
    'integration_stable': 'intergration_stable',
    'v2_subclass_full': 'celltype_full',
    'v2_subclass_l3': 'celltype_3',
    'v2_subclass_l2': 'celltype_2',
    'v2_subclass_sp': 'celltype_sp',
    'v2_subclass_l1': 'celltype_1',
    'v2_state_l2': 'cellstate_2',
    'v2_state_l1': 'cellstate_1',
    'v2_class': 'celltype_class',
    'v2_structure' : 'structure',
    #technical
    'source': 'tech_source',
    'assay': 'tech_assay',
    'experiment': 'tech_experiment',
    'percent_cortex': 'tech_percent_cortex',
    'percent_medulla': 'tech_percent_medulla',
    'region_level3': 'tech_region_3',
    'region_level2': 'tech_region_2',
    'region_level1': 'tech_region_1',
    'protocol': 'tech_protocol',
    'tissue_type_full': 'tech_tissue_type_full',
    'tissue_type': 'tech_tissue_type',
    'atlas_version': 'tech_atlas_version'
}

# %%
metadata = metadata.filter(column_mapping.keys()).rename(columns = column_mapping).set_index('cell_id')
metadata.index.name = None

# %%
adata.obs = metadata.loc[adata.obs.index,:]

# %%
patients_to_remove = ['3535'] #remove patient with two conditions
for patient in patients_to_remove:
    mu.pp.filter_obs(adata, 'uid_donor_id', lambda x: x != patient)

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

print(adata)

adata.write_h5ad(output_file)
