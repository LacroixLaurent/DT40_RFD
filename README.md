

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

*sessionInfo()*
R version 4.0.3 (2020-10-10)  
Platform: x86_64-apple-darwin17.0 (64-bit)  
Running under: macOS Catalina 10.15.7  

Matrix products: default  
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib  
LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib  

locale:  
[1] fr_FR.UTF-8/fr_FR.UTF-8/fr_FR.UTF-8/C/fr_FR.UTF-8/fr_FR.UTF-8  

attached base packages:  
[1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods  
[9] base     

other attached packages:  
 [1] BSgenome.Ggallus.UCSC.galGal6_1.4.2 openxlsx_4.2.2                     
 [3] GenomicAlignments_1.24.0            Rsamtools_2.4.0                    
 [5] SummarizedExperiment_1.18.2         DelayedArray_0.14.1                
 [7] matrixStats_0.57.0                  Biobase_2.48.0                     
 [9] BSgenome.Ggallus.UCSC.galGal4_1.4.0 BSgenome_1.56.0                    
[11] rtracklayer_1.48.0                  Biostrings_2.56.0                  
[13] XVector_0.28.0                      GenomicRanges_1.40.0               
[15] GenomeInfoDb_1.24.2                 IRanges_2.22.2                     
[17] S4Vectors_0.26.1                    BiocGenerics_0.34.0   

loaded via a namespace (and not attached):  
 [1] Rcpp_1.0.5             rstudioapi_0.11        zlibbioc_1.34.0       
 [4] BiocParallel_1.22.0    lattice_0.20-41        tools_4.0.3           
 [7] grid_4.0.3             crayon_1.3.4           zip_2.1.1             
[10] Matrix_1.2-18          GenomeInfoDbData_1.2.3 BiocManager_1.30.10   
[13] bitops_1.0-6           RCurl_1.98-1.2         stringi_1.5.3         
[16] compiler_4.0.3         XML_3.99-0.5   
