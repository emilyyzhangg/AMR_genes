import anndata as ad
import scanpy as sc
import pandas as pd
import scipy as sp
import glob
import numpy as np
import pickle as pkl


"""
Discovers different cell types, and identifies the representative genes of each type.

Parameters
adata: stores the Cell x Gene matrix and all associated meta information
feature_space: which matrix of observations to use in generating the neighborhood graph (select from adata.obsm)
neighbor_method: method to generate neighborhood graph (only supports 'gauss' and knn)   

Return
Dataframe containing each cluster's genes ranked from most->least indicative of the cell state associated with the cluster 
"""
def mrvi_identify_cell_states(adata: ad.AnnData, feature_space: str, neighbor_method: str):
    
    # Assumes Bacdrop replicate 1 cell read depth of 5,000
    sc.pp.filter_cells(adata, min_genes=35)
    sc.pp.normalize_total(adata, target_sum=5000)
    sc.pp.log1p(adata)
    
    adata.raw = adata

    if (neighbor_method == 'gauss'):
        sc.pp.neighbors(adata, use_rep=feature_space, knn=False, method='gauss')
    else:
        sc.pp.neighbors(adata, use_rep=feature_space, knn=True)

    sc.tl.leiden(adata)
    sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
    df = pd.DataFrame(adata.uns['rank_genes_groups']['names'])

    return df

"""
Removes duplicate genes by creating a new dataframe with columns of length n_genes each 
containing the first n_genes non-null entries from the given dataframe.  
"""
def filter_duplicates(df: pd.DataFrame, n_genes: int):
    
    d = {}
    # 1) Find and filter duplicates 
    for c in df.columns:
        # For the entries in each column
        potential_duplicates = set()
        for i in range(len(df[c])):
            temp = df[c][i].split('-')
            
            if len(temp) > 2:
                # replace with empty string
                if temp[1] in potential_duplicates:
                    df[c][i] = ''
                else:
                    potential_duplicates.add(temp[1])
        
    for c in df.columns:
        # 2) Return the top genes while ignoring the duplicates
        d[c] = []
        max_itr = len(df[c])
        itr = 0
        while (max_itr > itr) and (len(d[c]) < n_genes):
            if df[c][itr] != '':
                d[c].append(df[c][itr])
            
            itr += 1
    res = pd.DataFrame(data=d)
    return res


"""
Save genes mrvi finds to a text file.

Parameters
ranked_groups: a Dataframe output by @mrvi_identify_cell_states
outfile: file path of where to save results to
n_genes: the number of top genes to save from all the clusters (default is all)   

Return
"""
def write_cell_state_genes(ranked_groups: pd.DataFrame, outfile: str, n_genes=0):

    if n_genes == 0:
        ranked_groups.to_csv(outfile)
    else:
        df = ranked_groups.head[n_genes]
        df.to_csv(outfile)




"""
Saves files in any text format to anndata objects in h5ad file format.

Parameters
input_file: path to the file to save
output_file: desired location to save the results
delimiter: seperator used between entries in the text file (eg. ',' for csv)

Return
"""
def text_to_adata(input_file: str, output_file: str, delimiter='\t'):
    adata = sc.read_csv(input_file, delimiter=delimiter)
    adata.X = sp.sparse.csr_matrix(adata.X, dtype=np.float32)
    adata = adata.transpose()
    # TODO: replace duplicate names here, and re-run the entire pipeline ?
    adata.write_h5ad(output_file)


"""
Saves all files with given format in a directory to adata objects in h5ad format to specified directory. 
"""
def directory_to_adata(data_dir: str, output_dir: str):
    for filepath in glob.iglob(data_dir + "*.tsv"):
        outfile = output_dir + (filepath.split('/')[-1]).split('.')[0] + ".h5ad"
        text_to_adata(filepath, outfile)


"""
Convert and save the given object as a bytestream for quick future reloads. 
"""
def save_pickle(object, save_path: str):
    with open(save_path, 'wb') as outstream:
        pkl.dump(obj, outstream)

"""
Load in a pickled object. 
"""
def load_pickle(path: str):
    with open(path, 'rb') as instream:
      object_name = pkl.load(instream)
    return object_name


"""
Combines all adata objects in a given directory.
(merges observations from different conditions/batches into a single sample-aware datastructure) 

Parameters
data_dir: path to directory housing all the .h5ad files to merge.
union: keep all genes (T) or only those found in each sample matrix (F).
unique_names: if true hyphenates each cell's barcode with its corresponding sample label.
obs_label: the name to use for the adata.obs column that will be created in the merged adata object.

Return
The merged adata object.
"""
def concatenate_samples(data_dir: str, union=True, unique_names=False, obs_label='sample'):
    adatas = []
    sample_labels = []
    for filepath in glob.iglob(data_dir + "*.h5ad"):
        filename = filepath.split('/')[-1]
        sample_labels.append(filename.split('.')[0])
        adatas.append(sc.read(data_dir + filename))

    join_val = 'outer' if union else 'inner'
    iu_val = sample_labels if unique_names else None

    adata_temp = ad.concat(adatas, 
                            label=obs_label, 
                            keys=sample_labels,
                            axis=0,
                            join=join_val,
                            index_unique=iu_val, 
                            fill_value=np.float32(0.0))

    adata_temp.var_names_make_unique()
    return adata_temp
