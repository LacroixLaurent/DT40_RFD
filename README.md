

# Data analysis for DT40 Single Molecule RFD studies
## Laurent Lacroix (laurent.lacroix@inserm.fr)


### Scripts used for figures and analysis in Blin et al 2020

Scripts were used in R 4.0.3  
Blin_et_al2020.r describes all the analyses steps and produces bigwigs and plots used for the figures as well as the correlation coefficients reported in the article. 
Blin_et_al2020_Function.r describes the homemade functions used in the script.  
Original data are in the *data/* folder.  
This folder contains the following files:  
- *Track_coordinates.xlsx*: multisheet xlsx file containing the genomic coordinates of the oriented tracks for the 4 Regions Of Interest (wtCCSER1, betaCCSER1, wtDMD, tetDMD) from the single molecule experiments reported in Blin *et al.* NSMB 2019. [1 file]  
- *Ini/Ter_ROI.bed*: Initiation and Termination events for the 4 ROI in a bed format (chromosome, start, end). [8 files]  
- *ROI_OKreads_exp.bed*: OK-seq reads mapping on the extended ROI in a bed format (chromosome, start, end, strand). [2 files]  
- *galGal6_DMD_OKreads_exp.bed*: OK-seq reads for the DMD locus mapped on the galGal6 genome version. [1 file]  

Data were mapped either on *galGal4* or *galGal6* genome using *bwa mem -M*, converted into bam with *samtools view -b -F2304 -q10*, sorted with *samtools sort* and indexed with *samtools index* then imported into *R* and saved as bed files with the export function from the *rtracklayer* package.  
(samtools v1.4.1 ; bwa v0.7.15-r1140)  
Typical results are in the *results/* folder.  
