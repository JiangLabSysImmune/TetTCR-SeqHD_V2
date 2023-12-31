tetramer_per_umi_per_position <- read.delim("./tetramer_per_umi_per_position.tsv")
freq <- tetramer_per_umi_per_position$counts * tetramer_per_umi_per_position$instances_post / sum(tetramer_per_umi_per_position$counts * tetramer_per_umi_per_position$instances_post)
plot(x=tetramer_per_umi_per_position$counts,y= freq, log="x",xlab="read counts",ylab="percentage of reads among total reads",pch=19)
d1 <- diff(tetramer_per_umi_per_position$counts * tetramer_per_umi_per_position$instances_post)/diff(tetramer_per_umi_per_position$counts)
cutoffs <- which.min(abs(d1[1:10]))
cat(cutoffs)
