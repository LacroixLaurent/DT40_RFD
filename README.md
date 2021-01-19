

# Data analysis for DT40 RFD studies
## Laurent Lacroix (laurent.lacroix@inserm.fr)


### Scripts used for figures and analyses in Blin et al 2021 (NAR)

Scripts were used in R 4.0.3.  
Blin_et_al2021_data_preprocess.sh describes the preprocessing step used on the OKseq fastq data from the GEO dataset GSE164252.  
Blin_et_al2021.r describes all the analysis steps and produces bigwigs and plots used for the figures as well as the correlation coefficients reported in the article. 
Blin_et_al2021_Function.r describes the homemade functions used in the script.  
Original data are in the *data/* folder.  
This folder contains the following files:  
- *Track_coordinates.xlsx*: multisheet xlsx file containing the genomic coordinates of the oriented tracks for the 4 **R**egions **O**f **I**nterest (wtCCSER1, betaCCSER1, wtDMD, tetDMD) from the single molecule experiments reported in Blin *et al.* NSMB 2019. [1 file]  
- *Ini/Ter_ROI.bed*: Initiation and Termination events for the 4 ROI in a bed format (chromosome, start, end). [8 files]  
- *O1/E14/E15_ROI_OKreads_exp.bed*: OK-seq reads mapping on the extended wtCCSER1 and wtDMD loci in a bed format (chromosome, start, end, strand) for the trhee replicates (O1, E14, E15). [6 files]  

Data were mapped on *galGal4* genome using *bwa mem -M*, converted into bam with *samtools view -b -F2304 -q10*, sorted with *samtools sort* and indexed with *samtools index* then imported into *R*  as GenomicAlignmentsPairs, convert into GenomicRanges. Duplicates reads were removed and resulting reads were saved as bed files with the export function from the *rtracklayer* package (see Blin_et_al2020_data_preprocess.sh).  
(cutadapt 2.10 ; samtools v1.4.1 ; bwa v0.7.15-r1140)  
Typical results are in the *results/* folder.  

*session_info()* : see *session_info.txt*
