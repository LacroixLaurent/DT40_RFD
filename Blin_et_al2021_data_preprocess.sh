#### Blin et al 2021
#### Data preprocessing
#### LL 20210119

## genome from http://hgdownload.cse.ucsc.edu/goldenPath/galGal4/bigZips/galGal4.fa.gz
## raw data are from the GEO deposit (GSE164252)

adap1=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
adap2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

cutadapt -e 0.1 -q 20 -O 3 -a $adap1 -A $adap2 -o O1-DT40_S3_DMD.R1.cut.fastq -p O1-DT40_S3_DMD.R2.cut.fastq O1-DT40_S3_DMD_R1_001.fastq O1-DT40_S3_DMD_R2_001.fastq
bwa mem -M galGal4.fa.gz O1-DT40_S3_DMD.R1.cut.fastq O1-DT40_S3_DMD.R2.cut.fastq > O1-DT40_S3_DMDgal4.new.sam
samtools view -b -F2304 -q10 -o O1-DT40_S3_DMDgal4.new.bam O1-DT40_S3_DMDgal4.new.sam
samtools sort -o O1-DT40_S3_DMDgal4.new.sorted.bam O1-DT40_S3_DMDgal4.new.bam
samtools index O1-DT40_S3_DMDgal4.new.sorted.bam
rm *.new.sam
rm *.new.bam

cutadapt -e 0.1 -q 20 -O 3 -a $adap1 -A $adap2 -o O1-DT40_S3_CCSER1.R1.cut.fastq -p O1-DT40_S3_CCSER1.R2.cut.fastq O1-DT40_S3_CCSER1_R1_001.fastq O1-DT40_S3_CCSER1_R2_001.fastq
bwa mem -M genomes/galGal4.fa.gz O1-DT40_S3_CCSER1.R1.cut.fastq O1-DT40_S3_CCSER1.R2.cut.fastq > O1-DT40_S3_CCSER1gal4.new.sam
samtools view -b -F2304 -q10 -o O1-DT40_S3_CCSER1gal4.new.bam O1-DT40_S3_CCSER1gal4.new.sam
samtools sort -o O1-DT40_S3_CCSER1gal4.new.sorted.bam O1-DT40_S3_CCSER1gal4.new.bam
samtools index O1-DT40_S3_CCSER1gal4.new.sorted.bam
rm *.new.sam
rm *.new.bam

cutadapt -e 0.1 -q 20 -O 3 -a $adap1 -A $adap2 -o E14-DT40_S1_DMD.R1.cut.fastq -p E14-DT40_S1_DMD.R2.cut.fastq E14-DT40_S1_DMD_R1_001.fastq E14-DT40_S1_DMD_R2_001.fastq
bwa mem -M galGal4.fa.gz E14-DT40_S1_DMD.R1.cut.fastq E14-DT40_S1_DMD.R2.cut.fastq > E14-DT40_S1_DMDgal4.new.sam
samtools view -b -F2304 -q10 -o E14-DT40_S1_DMDgal4.new.bam E14-DT40_S1_DMDgal4.new.sam
samtools sort -o E14-DT40_S1_DMDgal4.new.sorted.bam E14-DT40_S1_DMDgal4.new.bam
samtools index E14-DT40_S1_DMDgal4.new.sorted.bam
rm *.new.sam
rm *.new.bam

cutadapt -e 0.1 -q 20 -O 3 -a $adap1 -A $adap2 -o E14-DT40_S1_CCSER1.R1.cut.fastq -p E14-DT40_S1_CCSER1.R2.cut.fastq E14-DT40_S1_CCSER1_R1_001.fastq E14-DT40_S1_CCSER1_R2_001.fastq
bwa mem -M galGal4.fa.gz E14-DT40_S1_CCSER1.R1.cut.fastq E14-DT40_S1_CCSER1.R2.cut.fastq > E14-DT40_S1_CCSER1gal4.new.sam
samtools view -b -F2304 -q10 -o E14-DT40_S1_CCSER1gal4.new.bam E14-DT40_S1_CCSER1gal4.new.sam
samtools sort -o E14-DT40_S1_CCSER1gal4.new.sorted.bam E14-DT40_S1_CCSER1gal4.new.bam
samtools index E14-DT40_S1_CCSER1gal4.new.sorted.bam
rm *.new.sam
rm *.new.bam

cutadapt -e 0.1 -q 20 -O 3 -a $adap1 -A $adap2 -o E15-DT40_S2_DMD.R1.cut.fastq -p E15-DT40_S2_DMD.R2.cut.fastq E15-DT40_S2_DMD_R1_001.fastq E15-DT40_S2_DMD_R2_001.fastq
bwa mem -M galGal4.fa.gz E15-DT40_S2_DMD.R1.cut.fastq E15-DT40_S2_DMD.R2.cut.fastq > E15-DT40_S2_DMDgal4.new.sam
samtools view -b -F2304 -q10 -o E15-DT40_S2_DMDgal4.new.bam E15-DT40_S2_DMDgal4.new.sam
samtools sort -o E15-DT40_S2_DMDgal4.new.sorted.bam E15-DT40_S2_DMDgal4.new.bam
samtools index E15-DT40_S2_DMDgal4.new.sorted.bam
rm *.new.sam
rm *.new.bam

cutadapt -e 0.1 -q 20 -O 3 -a $adap1 -A $adap2 -o E15-DT40_S2_CCSER1.R1.cut.fastq -p E15-DT40_S2_CCSER1.R2.cut.fastq E15-DT40_S2_CCSER1_R1_001.fastq E15-DT40_S2_CCSER1_R2_001.fastq
bwa mem -M galGal4.fa.gz E15-DT40_S2_CCSER1.R1.cut.fastq E15-DT40_S2_CCSER1.R2.cut.fastq > E15-DT40_S2_CCSER1gal4.new.sam
samtools view -b -F2304 -q10 -o E15-DT40_S2_CCSER1gal4.new.bam E15-DT40_S2_CCSER1gal4.new.sam
samtools sort -o E15-DT40_S2_CCSER1gal4.new.sorted.bam E15-DT40_S2_CCSER1gal4.new.bam
samtools index E15-DT40_S2_CCSER1gal4.new.sorted.bam
rm *.new.sam
rm *.new.bam

### the *.cut.fastq could be deleted or kept to map on another chicken genome is necessary
# rm *.cut.fastq

#### Data are then imported to R as pair end reads and deduplicated

require("BSgenome.Ggallus.UCSC.galGal4")
genome <- BSgenome.Ggallus.UCSC.galGal4
seqinf <- seqinfo(genome)
require(GenomicAlignments)

scanpar2 <- ScanBamParam(flag=scanBamFlag(isDuplicate=FALSE,isSecondaryAlignment=FALSE,isProperPair=TRUE),mapqFilter=10)

# O1 replicate
bf1 <- BamFile("O1-DT40_S3_DMDgal4.new.sorted.bam")
ga1 <- readGAlignmentPairs(bf1,param=scanpar2)
seqlevels(ga1,pruning.mode="coarse") <- seqlevels(seqinf)
seqinfo(ga1) <- seqinf
gra1 <- granges(ga1)
gra1 <- gra1[!duplicated(gra1)]
export(gra1,con="data/O1_galGal4_DMD_OKreads_exp.bed")

bf2 <- BamFile("O1-DT40_S3_CCSER1gal4.new.sorted.bam")
ga2 <- readGAlignmentPairs(bf2,param=scanpar2)
seqlevels(ga2,pruning.mode="coarse") <- seqlevels(seqinf)
seqinfo(ga2) <- seqinf
gra2 <- granges(ga2)
gra2 <- gra2[!duplicated(gra2)]
export(gra2,con="data/O1_galGal4_CCSER1_OKreads_exp.bed")

# E14 replicate
bf1 <- BamFile("E14-DT40_S1_DMDgal4.new.sorted.bam")
ga1 <- readGAlignmentPairs(bf1,param=scanpar2)
seqlevels(ga1,pruning.mode="coarse") <- seqlevels(seqinf)
seqinfo(ga1) <- seqinf
gra1 <- granges(ga1)
gra1 <- gra1[!duplicated(gra1)]
export(gra1,con="data/E14_galGal4_DMD_OKreads_exp.bed")

bf2 <- BamFile("E14-DT40_S1_CCSER1gal4.new.sorted.bam")
ga2 <- readGAlignmentPairs(bf2,param=scanpar2)
seqlevels(ga2,pruning.mode="coarse") <- seqlevels(seqinf)
seqinfo(ga2) <- seqinf
gra2 <- granges(ga2)
gra2 <- gra2[!duplicated(gra2)]
export(gra2,con="data/E14_galGal4_CCSER1_OKreads_exp.bed")

# E15 replicate
bf1 <- BamFile("E15-DT40_S2_DMDgal4.new.sorted.bam")
ga1 <- readGAlignmentPairs(bf1,param=scanpar2)
seqlevels(ga1,pruning.mode="coarse") <- seqlevels(seqinf)
seqinfo(ga1) <- seqinf
gra1 <- granges(ga1)
gra1 <- gra1[!duplicated(gra1)]
export(gra1,con="data/E15_galGal4_DMD_OKreads_exp.bed")

bf2 <- BamFile("E15-DT40_S2_CCSER1gal4.new.sorted.bam")
ga2 <- readGAlignmentPairs(bf2,param=scanpar2)
seqlevels(ga2,pruning.mode="coarse") <- seqlevels(seqinf)
seqinfo(ga2) <- seqinf
gra2 <- granges(ga2)
gra2 <- gra2[!duplicated(gra2)]
export(gra2,con="data/E15_galGal4_CCSER1_OKreads_exp.bed")
