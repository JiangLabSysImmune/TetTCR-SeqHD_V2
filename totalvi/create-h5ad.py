import pandas as pd
import numpy as np
import sys
import anndata as ad
import scanpy as sc

#input/output files
rna_file = sys.argv[1] #csv file where row = cells, col = features, data = RNA umi counts
ab_file = sys.argv[2] #same as above, but with protein umi counts
md_file = sys.argv[3] #csv file with all metadata, should have same indexes (cell ids) as rna and ab files
experiment = rna_file.split("/")[0]
out = experiment + ".h5ad"

def readin(f):
    """ read in csv to pandas and change cell label to string """
    """ because anndata can't handle numeric cell ids         """
    df = pd.read_csv(f,index_col=0,comment='#')
    df.index = experiment + df.index.astype(str)

    return(df)

#read in files
rna = readin(rna_file)
ab = readin(ab_file)
md = readin(md_file)

#create anndata
adata = ad.AnnData(rna,dtype=np.float32)
adata.obs = md
adata.obsm["protein_expression"] = ab

#write in h5ad
adata.write(out)

#END
