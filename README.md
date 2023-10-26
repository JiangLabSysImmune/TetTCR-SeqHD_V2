# TetTCR-SeqHD_V2
Converts raw sequencing data of pMHC tetramer barcodes to a count matrix
Runs on data generated from BD Rhapsody Enhanced Bead kit (Change file names accordingly)<br><br>
UPDATED: 6/21/22 by Michael Malone

# Usage:
1. MAKE SURE THE PEPTIDE REFERENCE CSV FILE FORMAT IS CORRECT (and called peptide.csv).<br>
Example:<br>
name,Ntsequence,AAsequence,gene,category
IA2,TTCATTGACTCTTACATCTGCCAGGTT,SLSPLQAEL,IA2,self
Empty,AACCTGGTTCCGATGGTTGCTACCGTT,Empty,Empty,Empty
Note: for the Empty, Capitalize E!

2. Copy all the files in this directory to the working directory where all the fastq files are.

3. Run the run.sh script:<br>
`bash run.sh Rhapsody_dbec.csv tet_R1.fastq.gz tet_R2.fastq.gz > log &`<br>
Remember to send the stdout to log file, otherwise the last step will be interupted.

Note: you may need to change some of the variables in src/tet.sh (particularly the conda environments)
      you may also need to edit some of the variables in run.sh, depending on needs
      specifically update the conda source path on line 6 of run.sh

Input files:
  1. Rhapsody_dbec.csv = *_DBEC_MolsPerCell.csv (output from SBG)
  2. peptide.csv = peptide reference file (for each experiment). This file needs to be in execution directory.
  3. tet_R1.fastq.gz & tet_R2.fastq.gz = raw reads from tetramer libraries

# Dependencies:
UMI tools (https://github.com/CGATOxford/UMI-tools)<br>
seqtk (https://github.com/lh3/seqtk)<br>
cutadapt (https://cutadapt.readthedocs.io/en/stable/)<br>
R & Rscript<br>
      - installed in a conda environment (`$R_env` in `src/tet.sh`) w/ r-base 4<br>
python<br>
      - base w/ pulp, pandas, numpy installed<br>
      - `$conda_env` w/ samtools and javac installed<br>
umis (https://pypi.org/project/umis/)<br>

#END
