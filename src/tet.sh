read1=$1
read2=$2
core=$3
ref=$4
n_bases=10 #number of n bases after 5' handle in barcode seqs

source /central_scratch/mjmalone/tools/miniconda3/etc/profile.d/conda.sh
conda_env="R"
R_env="R-4"

if [ -f tet_statistics.csv ]; then
        rm -f tet_statistics.csv  
fi                                     

dir=$(pwd)
cnt1=$(zcat $read1|wc -l)
echo "total tet reads,$((cnt1/4))" >> tet_statistics.csv
umis fastqtransform --cores $core --separate_cb ref/transform.json $read1 $read2|gzip > tmp/tet_R2.transform.fastq.gz

umis cb_filter --nedit 2 --cores $core --bc1 ref/bc1.txt --bc2 ref/bc2.txt --bc3 ref/bc3.txt tmp/tet_R2.transform.fastq.gz|gzip > tmp/tet_R2.filter.fastq.gz

zcat tmp/tet_R2.filter.fastq.gz |grep --no-group-separator -A 3 -Ff truecell.txt |gzip > fq/tet_R2.true.fastq.gz

cnt2=$(zcat fq/tet_R2.true.fastq.gz|wc -l)
echo "total tet truecell reads,$((cnt2/4))" >> tet_statistics.csv

cutadapt -a ATGGACGACGACGACAAG...TAACGAAGCACCTCGCT -M 40 -m 10 -e 0.1 -o fq/tet_R2.trim.fastq.gz fq/tet_R2.true.fastq.gz
seqtk seq -a fq/tet_R2.trim.fastq.gz > fq/tet_R2.trim.tmp.fasta
seqtk trimfq -b $n_bases fq/tet_R2.trim.tmp.fasta > fq/tet_R2.trim.fasta

awk -F ',' '{if(NR>1){ printf ">%s\n%s\n",$1,$2 }}' $ref > ref/peptide.fasta

mkdir tet
perl src/fasta-splitter.pl --n-parts $core fq/tet_R2.trim.fasta
cd src/
conda activate $conda_env
for i in $(seq -w 1 $core)
do
	(java -XX:+UseSerialGC alignHam ../tet_R2.trim.part-"$i".fasta $dir/ref/peptide.fasta > ../tetramer_align_"$i") &
done
wait
conda deactivate
cd ../
cat tetramer_align_* > tet/tetramer_align

algn=$(wc -l < tet/tetramer_align)
fasta=$(grep ">" fq/tet_R2.trim.fasta|wc -l)
if [ $algn -ne $fasta ]; then
	echo "alignment error!" >> error.log
fi

rm -f tetramer_align_*
rm -f tet_R2.trim.part-*.fasta

awk -F ',' 'NR==1{print "@HD\tVN:1.0\tSO:coordinate"};NR>1{printf "@SQ\tSN:%s\tLN:50\n",$1}' $4 > ref/samheader

cnt3=$(awk '$2==0' OFS='\t' tet/tetramer_align|wc -l)
echo "total tet aligned reads,$cnt3" >> tet_statistics.csv

cd tet
awk '$2==0' OFS='\t' tetramer_align |cat $dir/ref/samheader - > tetramer_align.sam
conda activate $conda_env
samtools sort -O bam tetramer_align.sam > tetramer_align.sorted.bam
samtools index tetramer_align.sorted.bam
conda deactivate

###cluster UMIs and remove UMIs based on read count distribution
umi_tools group -I tetramer_align.sorted.bam --extract-umi-method umis --group-out=tetramer.group --output-bam --log=group.log -S tetramer.group.bam --gene-tag=GX --per-cell --per-gene
umi_tools dedup -I tetramer_align.sorted.bam --extract-umi-method umis --output-stats tetramer --log=dedup.log -S tetramer.dedup.bam --gene-tag=GX --per-cell --per-gene
conda activate $R_env
R CMD BATCH $dir/src/umi_count_filter.R

rm -f tetramer.dedup.bam
rm -f tetramer.group.bam

cutoff=$(Rscript $dir/src/cutoff.R)
echo $cutoff
awk -v thr="$cutoff" '$6 < thr' tetramer.group|awk '{print $1}' > tetramer.group.filterlist
conda deactivate
conda activate $conda_env
fgrep -v -w -f tetramer.group.filterlist tetramer_align.sam|samtools sort -O bam - > tetramer_align.dbec.sorted.bam
samtools index tetramer_align.dbec.sorted.bam

final=$(samtools view tetramer_align.dbec.sorted.bam|wc -l)          
filter=$(wc -l < tetramer.group.filterlist)
total=$(wc -l < tetramer.group)
if [ $final -ne $(($total - $filter -1)) ]; then               
        echo "tetramer alignment filter error!" >> error.log
fi                                          
conda deactivate

cd ../

#END
