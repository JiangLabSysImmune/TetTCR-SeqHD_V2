import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scvi
import os
import sys
from glob import glob

transform_batch = "S-COV_22-07-28" #batch with most abseq data used to impute expression on other cells/batches

if not os.path.exists('./out'):
    os.makedirs('./out')
else:
    sys.exit("Error! Analysis previously run in this directory.\nDelete or move 'out' directory to proceed")

#read in files
h5ad_files = glob("*.h5ad")
adata_list = [ad.read_h5ad(h5ad) for h5ad in h5ad_files]
#concatenate
if len(h5ad_files) > 1:
    adata = ad.concat(adata_list,join="outer")
else:
    adata = adata_list[0]
#replace NaN values with zeroes in protein data
adata.obsm["protein_expression"] = adata.obsm["protein_expression"].replace(np.nan,0)

#detect highly variable genes (may need to optimize n)
sc.pp.highly_variable_genes(
    adata,
    batch_key="experiment",
    flavor="seurat_v3",
    n_top_genes=350, 
    subset=True
)

#set up scvi
scvi.model.TOTALVI.setup_anndata(adata, batch_key="experiment", protein_expression_obsm_key="protein_expression")

#prepare and train model
model = scvi.model.TOTALVI(
    adata, 
    latent_distribution="normal", 
    n_layers_decoder=2
)
model.train()

#determine protein foreground and predicted expression based on input batch
adata.obsm["X_totalVI"] = model.get_latent_representation()
adata.obsm["protein_fg_prob"] = model.get_protein_foreground_probability(transform_batch=transform_batch)
#extract normalized rna + protein expression values
rna, protein = model.get_normalized_expression(transform_batch=transform_batch,n_samples=25, return_mean=True)

#write files
rna.to_csv('./out/RNA_totalvi.csv')
protein.to_csv('./out/Protein_totalvi.csv')
adata.obs.to_csv("./out/metadata_totalvi.csv")
pd.DataFrame(data = adata.X, index = adata.obs.index, columns = adata.var.index).to_csv('./out/mrna_counts.csv')
np.savetxt('./out/latent.csv',adata.obsm["X_totalVI"],delimiter=",")
adata.write('./out/totalvi-out.h5ad')

#END
