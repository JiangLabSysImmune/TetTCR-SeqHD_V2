rhapsody=$1 # BD output (*.csv)
tetR1=$2 # tetramer read1 (*.fq.gzip)
tetR2=$3 # tetramer read2 (*.fq.gzip)
threads=15

source /central_scratch/mjmalone/tools/miniconda3/etc/profile.d/conda.sh

python src/barcode_conversion.py "$rhapsody"
mkdir fq
mkdir tmp
mkdir out
mkdir figures
cp peptide.csv out/
bash src/tet.sh "$tetR1" "$tetR2" $threads peptide.csv
umi_tools count -I tet/tetramer_align.dbec.sorted.bam --extract-umi-method umis --gene-tag=GX --per-cell --method unique --wide-format-cell-counts -S tet/tetramer_dbec.tsv
python src/umitools_index_conversion.py tet/tetramer_dbec.tsv out/tetramer_dbec.csv

#END
