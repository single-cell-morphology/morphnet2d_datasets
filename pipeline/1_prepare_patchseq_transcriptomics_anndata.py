"""
Preprocess Tolias transcript data.

This code depends on patchseq_celltype.py - please run that first.

Steps:
1. Combine exon and intron counts.
2. Make AnnData object from Pandas DataFrame.
3. Preprocess using Scanpy.
"""

import os
import pandas as pd
import anndata
from scipy import sparse
import scanpy as sc

def combine_exon_intron_counts():
    transcript_dir = "data/raw/tolias/rnaseq"
    exon_path = "m1_patchseq_exon_counts.csv.gz"
    intron_path = "m1_patchseq_intron_counts.csv.gz"

    # Read exon and intron counts
    exon_count = pd.read_csv(os.path.join(transcript_dir, exon_path), compression="gzip", index_col=0)
    intron_count = pd.read_csv(os.path.join(transcript_dir, intron_path), compression="gzip", index_col=0)
    count = exon_count.add(intron_count, fill_value=0)
    count = count.T # cell x gene matrix
    return count

def make_adata(raw_count_df, cell_source):
    metadata = pd.read_csv("data/processed/metadata/tolias_allen_celltype_5764cells.csv", index_col=0)
    obs = metadata[metadata["cell_source"] == cell_source]
    obs = obs.reindex(raw_count_df.index) # Make sure cell ids are in same order
    var = pd.DataFrame(index = raw_count_df.columns, columns=["genes"])
    adata = anndata.AnnData(X=sparse.csr_matrix(raw_count_df.values), var=var, obs=obs)

    return adata

def preprocess_adata(adata):
    # Follow tutorial on https://docs.scvi-tools.org/en/stable/user_guide/notebooks/api_overview.html
    sc.pp.filter_genes(adata, min_counts=3)
    adata.layers["counts"] = adata.X.copy() # preserve counts
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata # freeze the state in `.raw`
    """
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes = 10000,
        subset=True,
        layer="counts",
        flavor="seurat_v3"
    )
    """
    return adata

def add_metadata(adata):
    """Add metadata to adata.obs: `has_image`, `use_train`, `to_infer`"""

    obs = adata.obs
    obs["has_image"] = 0
    obs["use_train"] = 0
    obs["to_infer"] = 0

    img_names = os.listdir("data/processed/step2_padded")
    img_cellids = [name.split(".")[0] for name in img_names] # 1135 cellids
    celltypes_to_consider = ["Pvalb", "Sst", "IT", "Sncg", "Vip", "Lamp5"]

    for index, row in obs.iterrows():
        if index in img_cellids:
            obs.at[index, "has_image"] = 1
        if obs.at[index, "has_image"] == 1 and row["celltype"] in celltypes_to_consider:
            obs.at[index, "use_train"] = 1
        if obs.at[index, "celltype"] in celltypes_to_consider:
            obs.at[index, "to_infer"] = 1

    adata.obs = obs
    return adata

def change_allen_id(df):
    allen_metadata_filename = "data/raw/allen/transcripts/20200625_patchseq_metadata_mouse.csv"
    allen_metadata = pd.read_csv(allen_metadata_filename)

    for index, row in df.iterrows():
        cid = allen_metadata[allen_metadata["transcriptomics_sample_id"] == index]["cell_specimen_id"]
        df.rename(index={index:str(cid.values[0])}, inplace=True)

    return df

def main():
    # Tolias
    tolias_raw_count = combine_exon_intron_counts()
    tolias = make_adata(tolias_raw_count, "tolias")
    # tolias = preprocess_adata(tolias)
    """
    AnnData object with n_obs × n_vars = 1329 × 10000 (42466 total)
    obs: 'celltype', 'cell_source', 'transcriptomics_id'
    var: 'genes', 'n_counts', 'highly_variable', 'highly_variable_rank', 'means', 'variances', 'variances_norm'
    uns: 'log1p', 'hvg'
    layers: 'counts'
    """


    # Allen
    allen_transcript_path = "data/raw/allen/transcripts/20200513_Mouse_PatchSeq_Release_count.v2.csv"
    allen_raw_count = pd.read_csv(allen_transcript_path, index_col=0)
    allen_raw_count = allen_raw_count.T
    allen_raw_count = change_allen_id(allen_raw_count)

    allen  = make_adata(allen_raw_count, "allen")
    # allen = preprocess_adata(allen)
    """
    AnnData object with n_obs × n_vars = 4435 × 10000 (45768 total)
    obs: 'celltype', 'cell_source', 'transcriptomics_id'
    var: 'genes', 'n_counts', 'highly_variable', 'highly_variable_rank', 'means', 'variances', 'variances_norm'
    uns: 'log1p', 'hvg'
    layers: 'counts'
    """

    # Combine Raw Counts
    patchseq = anndata.concat([tolias, allen], axis=0) # 5764 x 29158
    patchseq = preprocess_adata(patchseq)
    """
    AnnData object with n_obs × n_vars = 5764 × 3552
    obs: 'celltype', 'cell_source', 'transcriptomics_id'
    layers: 'counts'
    """
    patchseq = add_metadata(patchseq) # 5764 x 28728
    # patchseq_raw_count = pd.concat([tolias_raw_count, allen_raw_count], axis=0, join="inner")
    # patchseq_raw_adata = make_adata(patchseq)
    # Combine Tolias and Allen adata
    patchseq.write_h5ad("data/processed/gene_expression/patchseq_Xgenes.h5ad")
    """
    AnnData object with n_obs × n_vars = 5764 × 28728
    obs: 'celltype', 'cell_source', 'transcriptomics_id', 'has_image', 'use_train', 'to_infer'
    var: 'n_counts'
    uns: 'log1p'
    layers: 'counts'
    """

if __name__ == "__main__":
    main()
