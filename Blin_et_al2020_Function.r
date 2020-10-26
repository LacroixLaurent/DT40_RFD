## helper function for Blin et al 2020
## LL 20200917
## laurent.lacroix@inserm.fr

library(GenomicAlignments)

makeRFD <- function(gr,bs=1000,lr=1,na2zero=F,bin=T,datatype='OKseq',export=F,saveRData=F,retur=T,outname='myRFDdata',OKcheck=T)
{
### Import data from GenomicRanges with the convention that strand=- for right progressing forks and strand=+ for left progressing forks
### Output total, R anf L coverage as well as RFD. It also filters out low reads regions from the RFD and export the coordinates of those regions.
### Output also BigWig files for the 3 coverages and the low reads regions filtered RFD
require(GenomicRanges)
require(rtracklayer)
corpm <- 1
if (datatype=='OKseq')
{
# check for strand symetry
pl <- length(gr[strand(gr)=='+'])
mi <- length(gr[strand(gr)=='-'])
(corpm <- pl/mi)
}
# strand inversion if not classical OKseq data
if (datatype=='nanopore')
	{
# invert 'strand' as in the datafile, -=left progressing and +=right progressing
	stra <- strand(gr)
	str2 <- stra
	str2[stra=='+'] <- '-'
	str2[stra=='-'] <- '+'
	strand(gr) <- str2
	}
# optionnal binning of the data to compute coverages
if (bin)
	{
	bingen <- tileGenome(seqinfo(gr),tilewidth=bs, cut.last.tile.in.chrom=T)

	bingenp <- countOverlaps(bingen,gr[strand(gr)=='+'])
	bingenm <- countOverlaps(bingen,gr[strand(gr)=='-'])
	cv_L <- coverage(bingen,weight=bingenp)	# left progressing fork, plus for Xia's data, W in OKseq paper
	cv_R <- coverage(bingen,weight=bingenm)	# right progressing fork, minus for Xia's data, C in OKseq paper
	}
else
	{
	cv_L <- coverage(gr[strand(gr)=='+'])
	cv_R <- coverage(gr[strand(gr)=='-'])
	bs=1
	}
# in a whole genome OKseq experiment, plus and minus strand shold be balanced, if not a correction is implemented
if (OKcheck)
{
if (abs(corpm-1)>0.01) {print('data are imbalanced, correction required');print(corpm); cv_R <- cv_R*corpm}
}

# compute RFD and filter for regions where total coverage is below a defined (lr) threshold
cv <- cv_L+cv_R
RFD <- (cv_R-cv_L)/(cv_R+cv_L)
lr_index <- which(cv<=lr)
RFD2 <- RFD
RFD2[lr_index] <- NA
lowread <- do.call(c,lapply(1:length(seqlevels(gr)), function(i) { toto <-Views(cv[[i]],cv[[i]]<=lr);if (length(toto)!=0) {nono=GRanges(seqnames=names(cv)[i], IRanges(start(toto),end(toto)), strand='*', seqinfo=seqinfo(gr))}else{nono=GRanges()};return(nono)}))

# Option to set out NA to 0 in RFD for easier import of the files in IGV but beware that it could change correlation
naname <- '_wiNA'
if (na2zero)
	{
	RFD[is.na(RFD)] <- 0
	RFD2[is.na(RFD2)] <- 0
	naname <- '_noNA'
	}
# Option to export resuls into BigWig files
if (export)
	{
		rtracklayer::export(cv,con=paste0(outname,'_cov_tot_bs',bs/1000,'k_lr',lr,'.bw'))
		rtracklayer::export(cv_L,con=paste0(outname,'_cov_2left_bs',bs/1000,'k_lr',lr,'.bw'))
		rtracklayer::export(cv_R,con=paste0(outname,'_cov_2right_bs',bs/1000,'k_lr',lr,'.bw'))
		rtracklayer::export(RFD2,con=paste0(outname,'_RFD_bs',bs/1000,'k_lr',lr,naname,'.bw'))
	}
res <- list(cv,cv_L,cv_R,RFD,RFD2,lowread)
names(res) <- c('cv','cv_L','cv_R','RFD','RFD2','lowread')
# Option to save the results in a RData file
if (saveRData)
{
save(res,file=paste0(outname,'_bs',bs/1000,'k_lr',lr,naname,'.RData'))
}
# Option to output the result of the funtion in the current R session
if (retur) {return(res)}
}

newXLSX_track_import <- function(xlfile="Coordonnees_polarite.xlsx",feuille=1)
{
# function used to import tracks coordinate from the XLSX file
require("openxlsx")
require("GenomicRanges")
df1 <- read.xlsx(xlfile,sheet=feuille)
# set strand/fork orientation as for OKseq reads
df1$strand <- '*'
df1$strand[df1$Orientation=='G'] <- '+'
df1$strand[df1$Orientation=='D'] <- '-'

df1.GR <- sort(GRanges(seqnames=df1$Chromosome,ranges=IRanges(df1$start,df1$end),strand=df1$strand))
return(df1.GR)
}

### a function to check correlation genome-wide between RFD disregarding places with NA 
### for binned values
cor.rfd.test2 <- function(a,b,met='s',binned=T,bs0=bs,...)
{
if (binned==T)
{
a.b <- endoapply(a, function(x) {y <-x ;runLength(y) <- runLength(x)%/%bs0;return(y)})
b.b <- endoapply(b, function(x) {y <-x ;runLength(y) <- runLength(x)%/%bs0;return(y)})
a <- a.b
b<- b.b
}
cor.test(as.numeric(unlist(a)[!is.na(unlist(a)) & !is.na(unlist(b))]),as.numeric(unlist(b)[!is.na(unlist(a)) & !is.na(unlist(b))]),method=met)
}
### for a single gene/chromosome with binning
cor.rfd.test3 <- function(a,b,met='s',binned=T,bs0=bs,...)
{
if (binned==T)
{
a.b <- a
runLength(a.b) <- runLength(a)/bs0
b.b <- b
runLength(b.b) <- runLength(b)/bs0
a <- a.b
b<- b.b
}
cor.test(as.numeric(unlist(a)[!is.na(unlist(a)) & !is.na(unlist(b))]),as.numeric(unlist(b)[!is.na(unlist(a)) & !is.na(unlist(b))]),method=met)
}

