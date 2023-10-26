# Totalvi model with TetTCR-SeqHD

Use these scripts to merge TetTCR-SeqHD experiments (output from seurat.R) and use totalvi for batch-correction, protein expression background correction, integrating proteinand mRNA expression data, normalization, etc.
For more information on totalvi, including installation and usage, visit: https://docs.scvi-tools.org/en/stable/user_guide/models/totalvi.html.

# Usage:
1. Move mrna.csv, abseq.csv, and meta_data.csv files into a directory named after your experiment. Then use the following command:
   `python create-h5ad.py Experiment_1/mrna.csv Experiment_1/abseq.csv Expriment_1/meta_data.csv`
   to create an h5ad file slotting all of the data. You should run this command for all experiments you wish to integrate.
2. Run the command:
   `python total-vi.py`.
   Note: you will need to edit the `transform_batch` variable in this script to indicate which batch should be used to inform protein expression imputation. To remove this behavior see the docs at the site above. This will output several files in a directory called "out", which can be used for downstream analyses.

# Python Dependencies:
`pandas
numpy
anndata
scanpy
scvi`

#END
