import os
import shutil
import pandas as pd

def celltype_remap(celltype):
    """Re-map celltype in Allen and Tolias to be consistent."""
    cell_map = {
        "ET": "ET",
        "CT": "CT",
        "IT": "IT",
        "Lamp5": "Lamp5",
        "NP": "NP",
        "Pvalb": "Pvalb",
        "Sncg": "Sncg",
        "Sst": "Sst",
        "Vip": "Vip",
        "low quality": "low quality",
        "L2/3_IT": "IT",
        "Interneuron-CGE-Lamp5": "Lamp5", 
        "Interneuron-CGE-Vip": "Vip", 
        "Interneuron-MGE-Pvalb": "Pvalb", 
        "Interneuron-MGE-Sst": "Sst", 
        "Interneuron-MGE-Sst-Chodl": "Sst",
    }
    for i in celltype.index:
        cell_class_orig = celltype.at[i, "celltype"]
        celltype.at[i, "celltype"] = cell_map[cell_class_orig]
    
    return celltype

def prepare_allen_celltype():
    # Load annotated celltype information 
    celltype_filename = "data/raw/allen/celltype/Allen_PatchSeq_annotated_new.csv"
    allen_celltype_df = pd.read_csv(celltype_filename)
    allen_celltype_df.columns=["cell_id", "celltype"]

    transcriptomics_id = allen_celltype_df["cell_id"].tolist()
    transcriptomics_id = [tid.replace('.', '-') for tid in transcriptomics_id]
    allen_celltype_df["transcriptomics_id"] = transcriptomics_id
    
    # Find cell_id corresponding to transcriptomic_id
    allen_celltype_df["cell_id"] = ""
    metadata_filename = "data/raw/allen/transcripts/20200625_patchseq_metadata_mouse.csv"
    metadata = pd.read_csv(metadata_filename)
    

    for index, row in allen_celltype_df.iterrows():
        tid = row["transcriptomics_id"]
        cid = metadata[metadata["transcriptomics_sample_id"] == tid]["cell_specimen_id"]
        row["cell_id"] = str(cid.values[0])

    allen_celltype_df["cell_source"] = "allen"

    return allen_celltype_df

def prepare_tolias_celltype():
    metadata_filename = "data/raw/tolias/metadata/m1_patchseq_meta_data.csv"
    metadata = pd.read_csv(metadata_filename, sep="\t", index_col=0)
    tolias_celltype_df = metadata[["Cell", "RNA family"]]
    tolias_celltype_df.columns = ["cell_id", "celltype"]
    tolias_celltype_df["cell_source"] = "tolias"

    return tolias_celltype_df

def prepare_nuclei_celltype():
    import pdb
    pdb.set_trace()
    annot = pd.read_csv("data/raw/macosko_10x_nuclei_v3/cluster.annotation.csv")
    membership = pd.read_csv("data/raw/macosko_10x_nuclei_v3/cluster.membership.csv")
    membership.columns = ["cell_id", "cluster_id"]
    membership["celltype"] = membership["cluster_id"].map(annot.set_index("cluster_id")["subclass_label"].to_dict())
    membership["class_label"] = membership["cluster_id"].map(annot.set_index("cluster_id")["class_label"].to_dict())

    return membership

def main():
    allen = prepare_allen_celltype()
    tolias = prepare_tolias_celltype()
    # nuclei = prepare_nuclei_celltype()

    # Concatenate information
    celltype = pd.concat([tolias, allen], ignore_index=True)
    celltype["cell_id"] = celltype["cell_id"].astype(str)
    celltype = celltype.set_index("cell_id")

    # Remap celltype labels and save
    celltype = celltype_remap(celltype)
    celltype.to_csv("data/processed/metadata/tolias_allen_celltype_5764cells.csv")
    """
    Sst            1886
    Pvalb          1622
    Vip             873
    Lamp5           488
    IT              369
    Sncg            269
    CT              106
    low quality      97
    ET               49
    NP                5
    """

    """
    def select_cells_with_images(celltype, images_path):
        cellids = []
        for filename in os.listdir(images_path):
            cellid = "_".join(filename.split("_")[:-1])
            cellids.append(cellid)
        
        celltype = celltype.loc[celltype["cell_id"].isin(cellids)]
        return celltype
    # Select Cells that are in images folder
    celltype = select_cells_with_images(celltype, images_path)

    # Save
    if not os.path.exists(savepath):
        print(f"Creating directory {savepath}...")
        os.makedirs(savepath)

    celltype.to_csv(os.path.join(savepath, "celltype.csv"))
    """

if __name__ == "__main__":
    main()

    

   
