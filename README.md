

# Data analysis for DT40 Single Molecule RFD studies
## Laurent Lacroix (laurent.lacroix@inserm.fr)


### Scripts used for figures and analysis in Blin et al 2020

Scripts were used in R 3.6.2.  
Blin_et_al2020.r describes all the analysis steps and produces bigwig and plot used for the figures as well as the correlation coefficient reported in the article. 
Blin_et_al2020_Function.r describes the the homemade function used in the script.  
Original data are in the *data/* folder. 
This folder contains the following files:  
- *Track_coordinates.xlsx*: multisheet xlsx file containing the genomic coordinates of the oriented tracks for the 4 Region Of Interest (wtCCSER1, betaCCSER1, wtDMD,tetDMD) from the single molecule experiments reported in Blin et al NSMB 2019. [1 file]  
- *Ini/Ter_ROI.bed*: Initiation and Terminaison events for the 4 ROI in a bed format (chromosome, start, end). [8 files]  
- *ROI_OKreads_exp.bed*: OKseq reads mapping on the extended ROI in a bed format (chromosme, start, end, strand). [2 files]  
Typical results are in the *results/* folder.  
