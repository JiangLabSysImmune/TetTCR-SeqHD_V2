# TetTCR-SeqHD_V2
Converts raw sequencing data of pMHC tetramer barcodes to a count matrix
Runs on data generated from BD Rhapsody Enhanced Bead kit (Change file names accordingly)
UPDATED: 6/21/22 by Michael Malone

1.MAKE SURE THE PEPTIDE REFERENCE CSV FILE FORMAT IS CORRECT (and called peptide.csv).
Example:
name,Ntsequence,AAsequence,gene,category
IA2,TTCATTGACTCTTACATCTGCCAGGTT,SLSPLQAEL,IA2,self
Empty,AACCTGGTTCCGATGGTTGCTACCGTT,Empty,Empty,Empty
Note: for the Empty, Capitalize E!

2.copy all the files in this directory to the working directory where all the fastq files are.

3.run the run.sh script
bash run.sh Rhapsody_dbec.csv tet_R1.fastq.gz tet_R2.fastq.gz > log &
Remember to send the stdout to log file, otherwise the last step will be interupted.

Note: you may need to change some of the variables in src/tet.sh (particularly the conda environments)
      you may also need to edit some of the variables in run.sh, depending on needs

Input files:
  1. Rhapsody_dbec.csv = *_DBEC_MolsPerCell.csv (output from SBG)
  2. peptide.csv = peptide reference file (for each experiment)
  3. tet_R1.fastq.gz & tet_R2.fastq.gz = raw reads from tetramer libraries

# It is recommended that you run your data through more stringent and educated single cell workflows.
# Most analyses above should only be used preliminarily

# Dependencies:
  1. umi_tools
  2. seqtk
  3. cutadapt
  4. java
  5. R & Rscript
       - installed in a conda environment ($R_env in src/tet.sh) w/ r-base 4
  6. python
       - base w/ pulp, pandas, numpy installed
       - $conda_env w/ samtools and javac installed
  7. umis

#END
