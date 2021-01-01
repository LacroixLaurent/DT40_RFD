#### Data treatment for Blin et al 2020
### LL 20210101
### laurent.lacroix@inserm.fr

library("BSgenome.Ggallus.UCSC.galGal4")
genome <- BSgenome.Ggallus.UCSC.galGal4
seqinf <- seqinfo(genome)

### Set your working path
mypath <- "TypeYourPathHere"
setwd(mypath)

### Import home made function
source(paste0(mypath,"Blin_et_al2020_Function.r"))


### set coordinates for the region of interest
DMD <- GRanges(seqnames="chr1",ranges=IRanges(start=113880001,end=114960000), strand="*",seqinfo=seqinf)
CCSER1 <- GRanges(seqnames="chr4",ranges=IRanges(start=34950001,end=35660000), strand="*",seqinfo=seqinf)

## Import OKseq reads
gra_dmd <- import("data/galGal4_DMD_OKreads_exp.bed")
seqlevels(gra_dmd,pruning.mode="coarse") <- seqlevels(seqinf)
seqinfo(gra_dmd) <- seqinf
gra_ccser1 <- import("data/galGal4_CCSER1_OKreads_exp.bed")
seqlevels(gra_ccser1,pruning.mode="coarse") <- seqlevels(seqinf)
seqinfo(gra_ccser1) <- seqinf

gr_dmd <- gra_dmd[gra_dmd %within% DMD]
gr_ccser1 <- gra_ccser1[gra_ccser1 %within% CCSER1]

### work on DMD

# compute RFD for OKseq and SM  without binning and with 10 kb bins

smwtdmd <- newXLSX_track_import(xlfile="data/Track_coordinates.xlsx",feuille=3)
seqlevels(smwtdmd, pruning.mode="coarse") <- seqlevels(seqinf)
seqinfo(smwtdmd) <- seqinf
smwtdmd.g4 <- smwtdmd[smwtdmd %within% DMD]

sm.wtdmd.nt <- makeRFD(smwtdmd.g4,lr=5,na2zero=F,bin=F,datatype="OKseq",export=T,saveRData=F,retur=T,outname="results/sm_wtDMD_galGal4",OKcheck=F)
sm.wtdmd.10k <- makeRFD(smwtdmd.g4,lr=5,na2zero=F,bin=T,bs=10000,datatype="OKseq",export=T,saveRData=F,retur=T,outname="results/sm_wtDMD_galGal4",OKcheck=F)


smtetdmd <- newXLSX_track_import(xlfile="data/Track_coordinates.xlsx",feuille=4)
seqlevels(smtetdmd, pruning.mode="coarse") <- seqlevels(seqinf)
seqinfo(smtetdmd) <- seqinf
smtetdmd.g4 <- smtetdmd[smtetdmd %within% DMD]

sm.tetdmd.nt <- makeRFD(smtetdmd.g4,lr=5,na2zero=F,bin=F,datatype="OKseq",export=T,saveRData=F,retur=T,outname="results/sm_tetDMD_galGal4",OKcheck=F)
sm.tetdmd.10k <- makeRFD(smtetdmd.g4,lr=5,na2zero=F,bin=T,bs=10000,datatype="OKseq",export=T,saveRData=F,retur=T,outname="results/sm_tetDMD_galGal4",OKcheck=F)

# and OKseq
## larger DMD ROI
res_10k_dmdlarge <- makeRFD(gra_dmd,bs=10000,lr=20,na2zero=F,bin=T,datatype="OKseq",export=T,saveRData=F,retur=T,outname="results/DT40_OK_DMDlarge_galGal4",OKcheck=F)
## and for just DMD
res_10k_dmd <- makeRFD(gr_dmd,bs=10000,lr=20,na2zero=F,bin=T,datatype="OKseq",export=T,saveRData=F,retur=T,outname="results/DT40_OK_DMD_galGal4",OKcheck=F)

### compute correlation
cor.rfd.test2(res_10k_dmd[[4]], sm.wtdmd.10k[[4]],bs0=10000,binned=T)
#	Spearman's rank correlation rho
#data:  as.numeric(unlist(a)[!is.na(unlist(a)) & !is.na(unlist(b))]) and as.numeric(unlist(b)[!is.na(unlist(a)) & !is.na(unlist(b))])
#S = 35279, p-value < 2.2e-16
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#     rho
#0.831951


### try to evaluate the subsampling required to get the same coverage as the combing RFD

sm.dmd.cv.10k <- sm.wtdmd.10k[[1]][DMD][[1]]
summary(sm.dmd.cv.10k)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   14.0    21.0    26.0    24.9    29.0    36.0


### 100 sampling taking 2500 reads mapping on the ROI
n=2500
bs0=10000
set.seed(42)
res_list10k <- mclapply(1:100, function(i)
{gr_dmd_sub <- sample(gr_dmd,n)
resall_dmd_sub <- makeRFD(gr_dmd_sub,bs=bs0,lr=0,na2zero=F,bin=T,datatype="OKseq",export=F,saveRData=F,retur=T,OKcheck=F)},mc.cores=4L)
save(res_list10k,file="res_list10k_DMD_galGal4.RData")

dmd.rfds10k <- lapply(res_list10k, function(x) {
	y <- x[[4]][DMD][[1]]
	z <- y
	runLength(z) <- runLength(y)/bs0
	return(z)
})

dmd.rfds10k.cov <- lapply(res_list10k, function(x) {
	y <- x[[1]][DMD][[1]]
	z <- y
	runLength(z) <- runLength(y)/bs0
	return(z)
})
smdata <- sm.wtdmd.10k[[4]][DMD][[1]]
smdata2 <- smdata
runLength(smdata2) <- runLength(smdata)/bs0

smdata.cov <- sm.wtdmd.10k[[1]][DMD][[1]]
smdata2.cov <- smdata.cov
runLength(smdata2.cov) <- runLength(smdata.cov)/bs0

xmin <- 113880000/bs0
xmax <-114960000/bs0-1

pdf(file="results/subsample_DMD_10k_galGal4.pdf",height=4,width=6)
plot(smdata2,x=(xmin:xmax),col="red",type="l",lwd=2,ylim=c(-1,1),ylab="RFD",xlab=paste0("coordinates (chr1, *",bs0/1000,"kb)"))
dummy <- lapply(dmd.rfds10k, function(y) lines(y,x=(xmin:xmax),col=rgb(0,0,1,0.2)))
lines(smdata2,x=(xmin:xmax),col="red",type="l",lwd=2)
abline(h=0)
legend("topright",legend=c(paste0("subsampled(n=",n,")"),"single molecule RFD"), text.col=c("blue","red"),bty="n")

plot(smdata2.cov,x=(xmin:xmax),col="red",type="l",lwd=2,ylim=c(0,70),ylab="Coverage",xlab=paste0("coordinates (chr1, *",bs0/1000,"kb)"))
dummy <- lapply(dmd.rfds10k.cov, function(y) lines(y,x=(xmin:xmax),col=rgb(0,0,1,0.2)))
lines(smdata2.cov,x=(xmin:xmax),col="red",type="l",lwd=2)
legend("topright",legend=c(paste0("subsampled(n=",n,")"),"single molecule"), text.col=c("blue","red"),bty="n")
dev.off()


## REM (I-T profiles, Replication Efficiency Metric)
# WT DMD
ori_dmd <- import("data/Ini_wtDMD.bed")
ter_dmd <- import("data/Ter_wtDMD.bed")
seqinfo(ori_dmd) <- seqinf
seqinfo(ter_dmd) <- seqinf
ori_dmd <- resize(ori_dmd,fix="center",width=1)
ter_dmd <- resize(ter_dmd,fix="center",width=1)

bs=10000
bin10k <- tileGenome(seqinf,tilewidth=bs, cut.last.tile.in.chrom=T)
dmd10k <- bin10k[overlapsAny(bin10k,DMD)]
dmd10ki <- binnedAverage(dmd10k,coverage(ori_dmd),"ni",na.rm=T)
dmd10kt <- binnedAverage(dmd10k,coverage(ter_dmd),"nt",na.rm=T)
dmd10k$ni <- dmd10ki$ni*bs
dmd10k$nt <- dmd10kt$nt*bs
dmd10k$rem <- dmd10k$ni-dmd10k$nt
cvdmd10krem_wt <- coverage(dmd10k,weight=dmd10k$rem)
export(cvdmd10krem_wt,con="results/REM_WTDMD_10k_galGal4.bw")

# TET DMD
ori_dmd <- import("data/Ini_tetDMD.bed")
ter_dmd <- import("data/Ter_tetDMD.bed")
seqinfo(ori_dmd) <- seqinf
seqinfo(ter_dmd) <- seqinf
ori_dmd <- resize(ori_dmd,fix="center",width=1)
ter_dmd <- resize(ter_dmd,fix="center",width=1)

bs=10000
bin10k <- tileGenome(seqinf,tilewidth=bs, cut.last.tile.in.chrom=T)
dmd10k <- bin10k[overlapsAny(bin10k,DMD)]
dmd10ki <- binnedAverage(dmd10k,coverage(ori_dmd),"ni",na.rm=T)
dmd10kt <- binnedAverage(dmd10k,coverage(ter_dmd),"nt",na.rm=T)
dmd10k$ni <- dmd10ki$ni*bs
dmd10k$nt <- dmd10kt$nt*bs
dmd10k$rem <- dmd10k$ni-dmd10k$nt
cvdmd10krem_tet <- coverage(dmd10k,weight=dmd10k$rem)
export(cvdmd10krem_tet,con="results/REM_TETDMD_10k_galGal4.bw")

## OEM (Origin Efficiency Metrics)
## at 10k on 10k RFD
## for SM
win=10000
# WT DMD
cvL <- sm.wtdmd.10k[[2]]
cvT <- sm.wtdmd.10k[[1]]
cvLn <- cvL/cvT
oem0 <- endoapply(cvLn, function(cv) c(Rle(rep(NA,(win))),cv)-c(cv,Rle(rep(NA,(win)))))
oem1 <- endoapply(oem0, function(x) x[(1+win/2):(length(x)-(win/2))])
oem10k_wtDMD10k <- oem1
export(oem10k_wtDMD10k,con="results/OEM_WTDMD_10k_galGal4.bw")

oem.t <- oem10k_wtDMD10k[DMD-win][[1]]
rem.t <- cvdmd10krem_wt[DMD-win][[1]]
cor.rfd.test3(oem.t,rem.t,binned=T,bs0=5000)
#	Spearman's rank correlation rho
#data:  as.numeric(unlist(a)[!is.na(unlist(a)) & !is.na(unlist(b))]) and as.numeric(unlist(b)[!is.na(unlist(a)) & !is.na(unlist(b))])
#S = 906425, p-value = 6.527e-11
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho
#0.4291981

# TET DMD
cvL <- sm.tetdmd.10k[[2]]
cvT <- sm.tetdmd.10k[[1]]
cvLn <- cvL/cvT
oem0 <- endoapply(cvLn, function(cv) c(Rle(rep(NA,(win))),cv)-c(cv,Rle(rep(NA,(win)))))
oem1 <- endoapply(oem0, function(x) x[(1+win/2):(length(x)-(win/2))])
oem10k_tetDMD10k <- oem1
export(oem10k_tetDMD10k,con="results/OEM_TETDMD_10k_galGal4.bw")

oem.t <- oem10k_tetDMD10k[DMD-win][[1]]
rem.t <- cvdmd10krem_tet[DMD-win][[1]]
cor.rfd.test3(oem.t,rem.t,binned=T,bs0=5000)
#	Spearman's rank correlation rho
#data:  as.numeric(unlist(a)[!is.na(unlist(a)) & !is.na(unlist(b))]) and as.numeric(unlist(b)[!is.na(unlist(a)) & !is.na(unlist(b))])
#S = 664755, p-value < 2.2e-16
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho
#0.5813848


## for CCSER1

### compute RFD for SM  without binning and with 10 kb bins

smwtccser1 <- newXLSX_track_import(xlfile="data/Track_coordinates.xlsx",feuille=1)
seqlevels(smwtccser1, pruning.mode="coarse") <- seqlevels(seqinf)
seqinfo(smwtccser1) <- seqinf
smwtccser1.g4 <- smwtccser1[smwtccser1 %within% CCSER1]

sm.wtccser1.nt <- makeRFD(smwtccser1.g4,lr=5,na2zero=F,bin=F,datatype="OKseq",export=T,saveRData=F,retur=T,outname="results/sm_wtCCSER1_galGal4",OKcheck=F)
sm.wtccser1.10k <- makeRFD(smwtccser1.g4,lr=5,na2zero=F,bin=T,bs=10000,datatype="OKseq",export=T,saveRData=F,retur=T,outname="results/sm_wtCCSER1_galGal4",OKcheck=F)

smbetaccser1 <- newXLSX_track_import(xlfile="data/Track_coordinates.xlsx",feuille=2)
seqlevels(smbetaccser1, pruning.mode="coarse") <- seqlevels(seqinf)
seqinfo(smbetaccser1) <- seqinf
smbetaccser1.g4 <- smbetaccser1[smbetaccser1 %within% CCSER1]

sm.betaccser1.nt <- makeRFD(smbetaccser1.g4,lr=5,na2zero=F,bin=F,datatype="OKseq",export=T,saveRData=F,retur=T,outname="results/sm_betaCCSER1_galGal4",OKcheck=F)
sm.betaccser1.10k <- makeRFD(smbetaccser1.g4,lr=5,na2zero=F,bin=T,bs=10000,datatype="OKseq",export=T,saveRData=F,retur=T,outname="results/sm_betaCCSER1_galGal4",OKcheck=F)

# and OKseq with 10kb bins
## larger CCSER1 ROI
res_10k_ccser1large <- makeRFD(gra_ccser1,bs=10000,lr=20,na2zero=F,bin=T,datatype="OKseq",export=T,saveRData=F,retur=T,outname="results/DT40_OK_CCSER1large_galGal4",OKcheck=F)
## just CCSER1
res_10k_ccser1 <- makeRFD(gr_ccser1,bs=10000,lr=20,na2zero=F,bin=T,datatype="OKseq",export=T,saveRData=F,retur=T,outname="results/DT40_OK_CCSER1_galGal4",OKcheck=F)

## compute correlation

cor.rfd.test2(res_10k_ccser1[[4]], sm.wtccser1.10k[[4]],bs0=10000,binned=T)
#	Spearman's rank correlation rho
#data:  as.numeric(unlist(a)[!is.na(unlist(a)) & !is.na(unlist(b))]) and as.numeric(unlist(b)[!is.na(unlist(a)) & !is.na(unlist(b))])
#S = 20433, p-value = 4.741e-10
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho
#0.6573985


### Subsampling
n=2500
bs0=10000
set.seed(42)
res_list10kc <- mclapply(1:100, function(i)
{gr_ccser1_sub <- sample(gr_ccser1,n)
resall_ccser1_sub <- makeRFD(gr_ccser1_sub,bs=bs0,lr=0,na2zero=F,bin=T,datatype="OKseq",export=F,saveRData=F,retur=T,OKcheck=F)},mc.cores=4L)
save(res_list10kc,file="res_list10k_CCSER1_galGal4.RData")

ccsser1.rfds10k <- lapply(res_list10kc, function(x) {
	y <- x[[4]][CCSER1][[1]]
	z <- y
	runLength(z) <- runLength(y)/bs0
	return(z)
})

ccsser1.rfds10k.cov <- lapply(res_list10kc, function(x) {
	y <- x[[1]][CCSER1][[1]]
	z <- y
	runLength(z) <- runLength(y)/bs0
	return(z)
})
smdata <- sm.wtccser1.10k[[4]][CCSER1][[1]]
smdata2 <- smdata
runLength(smdata2) <- runLength(smdata)/bs0

smdata.cov <- sm.wtccser1.10k[[1]][CCSER1][[1]]
smdata2.cov <- smdata.cov
runLength(smdata2.cov) <- runLength(smdata.cov)/bs0

xmin <- 34950000/bs0
xmax <- 35660000/bs0-1

pdf(file="results/subsample_CCSER1_10k_galGal4.pdf",height=4,width=6)
plot(smdata2,x=(xmin:xmax),col="red",type="l",lwd=2,ylim=c(-1,1),ylab="RFD",xlab=paste0("coordinates (chr4, *",bs0/1000,"kb)"))
dummy <- lapply(ccsser1.rfds10k, function(y) lines(y,x=(xmin:xmax),col=rgb(0,0,1,0.2)))
lines(smdata2,x=(xmin:xmax),col="red",type="l",lwd=2)
legend('topright',legend=c('subsampled(n=2500)','single molecule RFD'), text.col=c('blue','red'),bty='n')

plot(smdata2.cov,x=(xmin:xmax),col="red",type="l",lwd=2,ylim=c(0,70),ylab="Coverage",xlab=paste0("coordinates (chr4, *",bs0/1000,"kb)"))
dummy <- lapply(ccsser1.rfds10k.cov, function(y) lines(y,x=(xmin:xmax),col=rgb(0,0,1,0.2)))
lines(smdata2.cov,x=(xmin:xmax),col="red",type="l",lwd=2)
legend('topright',legend=c('subsampled(n=2500)','single molecule'), text.col=c('blue','red'),bty='n')
dev.off()

## REM (I-T profiles, Replication Efficiency Metric)
ori_ccser1 <- import("data/Ini_wtCCSER1.bed")
ter_ccser1 <- import("data/Ter_wtCCSER1.bed")
seqlevels(ori_ccser1,pruning.mode="coarse") <- seqlevels(seqinf)
seqinfo(ori_ccser1) <- seqinf
seqlevels(ter_ccser1,pruning.mode="coarse") <- seqlevels(seqinf)
seqinfo(ter_ccser1) <- seqinf
ori_ccser1 <- resize(ori_ccser1,fix="center",width=1)
ter_ccser1 <- resize(ter_ccser1,fix="center",width=1)


bs=10000
bin10k <- tileGenome(seqinf,tilewidth=bs, cut.last.tile.in.chrom=T)
ccser110k <- bin10k[overlapsAny(bin10k,CCSER1)]
ccser110ki <- binnedAverage(ccser110k,coverage(ori_ccser1),"ni",na.rm=T)
ccser110kt <- binnedAverage(ccser110k,coverage(ter_ccser1),"nt",na.rm=T)
ccser110k$ni <- ccser110ki$ni*bs
ccser110k$nt <- ccser110kt$nt*bs
ccser110k$rem <- ccser110k$ni-ccser110k$nt
cvccser110krem_wt <- coverage(ccser110k,weight=ccser110k$rem)
export(cvccser110krem_wt,con="results/REM_WTCCSER1_10k_galGal4.bw")


ori_ccser1 <- import("data/Ini_betaCCSER1.bed")
ter_ccser1 <- import("data/Ter_betaCCSER1.bed")
seqlevels(ori_ccser1,pruning.mode="coarse") <- seqlevels(seqinf)
seqinfo(ori_ccser1) <- seqinf
seqlevels(ter_ccser1,pruning.mode="coarse") <- seqlevels(seqinf)
seqinfo(ter_ccser1) <- seqinf
ori_ccser1 <- resize(ori_ccser1,fix="center",width=1)
ter_ccser1 <- resize(ter_ccser1,fix="center",width=1)

bs=10000
bin10k <- tileGenome(seqinf,tilewidth=bs, cut.last.tile.in.chrom=T)
ccser110k <- bin10k[overlapsAny(bin10k,CCSER1)]
ccser110ki <- binnedAverage(ccser110k,coverage(ori_ccser1),"ni",na.rm=T)
ccser110kt <- binnedAverage(ccser110k,coverage(ter_ccser1),"nt",na.rm=T)
ccser110k$ni <- ccser110ki$ni*bs
ccser110k$nt <- ccser110kt$nt*bs
ccser110k$rem <- ccser110k$ni-ccser110k$nt
cvccser110krem_beta <- coverage(ccser110k,weight=ccser110k$rem)
export(cvccser110krem_beta,con="results/REM_BETACCSER1_10k_galGal4.bw")

## OEM (Origin Efficiency Metric)
## at 10k on 10k RFD
## for SM
win=10000
cvL <- sm.wtccser1.10k[[2]]
cvT <- sm.wtccser1.10k[[1]]
cvLn <- cvL/cvT
oem0 <- endoapply(cvLn, function(cv) c(Rle(rep(NA,(win))),cv)-c(cv,Rle(rep(NA,(win)))))
oem1 <- endoapply(oem0, function(x) x[(1+win/2):(length(x)-(win/2))])
oem10k_wtCCSER1_10k <- oem1
export(oem10k_wtCCSER1_10k,con="results/OEM_WTCCSER1_10k_galGal4.bw")

oem.t <- oem10k_wtCCSER1_10k[CCSER1-win][[1]]
rem.t <- cvccser110krem_wt[CCSER1-win][[1]]
cor.rfd.test3(oem.t,rem.t,binned=T,bs0=5000)
#	Spearman's rank correlation rho
#data:  as.numeric(unlist(a)[!is.na(unlist(a)) & !is.na(unlist(b))]) and as.numeric(unlist(b)[!is.na(unlist(a)) & !is.na(unlist(b))])
#S = 198105, p-value = 3.604e-12
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho
#0.5476937

cvL <- sm.betaccser1.10k[[2]]
cvT <- sm.betaccser1.10k[[1]]
cvLn <- cvL/cvT
oem0 <- endoapply(cvLn, function(cv) c(Rle(rep(NA,(win))),cv)-c(cv,Rle(rep(NA,(win)))))
oem1 <- endoapply(oem0, function(x) x[(1+win/2):(length(x)-(win/2))])
oem10k_betaCCSER1_10k <- oem1
export(oem10k_betaCCSER1_10k,con="results/OEM_BETACCSER1_10k_galGal4.bw")

oem.t <- oem10k_betaCCSER1_10k[CCSER1-win][[1]]
rem.t <- cvccser110krem_beta[CCSER1-win][[1]]
cor.rfd.test3(oem.t,rem.t,binned=T,bs0=5000)
#	Spearman's rank correlation rho
#data:  as.numeric(unlist(a)[!is.na(unlist(a)) & !is.na(unlist(b))]) and as.numeric(unlist(b)[!is.na(unlist(a)) & !is.na(unlist(b))])
#S = 179042, p-value = 2.27e-14
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho
#0.5912191

### Adding two OK seq replicates
gra14_dmd <- import("data/E14_galGal4_DMD_OKreads_exp.bed")
seqlevels(gra14_dmd,pruning.mode="coarse") <- seqlevels(seqinf)
seqinfo(gra14_dmd) <- seqinf
gra14_ccser1 <- import("data/E14_galGal4_CCSER1_OKreads_exp.bed")
seqlevels(gra14_ccser1,pruning.mode="coarse") <- seqlevels(seqinf)
seqinfo(gra14_ccser1) <- seqinf
gra15_dmd <- import("data/E15_galGal4_DMD_OKreads_exp.bed")
seqlevels(gra15_dmd,pruning.mode="coarse") <- seqlevels(seqinf)
seqinfo(gra15_dmd) <- seqinf
gra15_ccser1 <- import("data/E15_galGal4_CCSER1_OKreads_exp.bed")
seqlevels(gra15_ccser1,pruning.mode="coarse") <- seqlevels(seqinf)
seqinfo(gra15_ccser1) <- seqinf

grE14_dmd <- graE14_dmd[graE14_dmd %within% DMD]
grE14_ccser1 <- graE14_ccser1[graE14_ccser1 %within% CCSER1]
grE15_dmd <- graE15_dmd[graE15_dmd %within% DMD]
grE15_ccser1 <- graE15_ccser1[graE15_ccser1 %within% CCSER1]

res_10k_dmdE14 <- makeRFD(grE14_dmd,bs=10000,lr=20,na2zero=F,bin=T,datatype='OKseq',export=F,saveRData=F,retur=T,outname='DT40_OK_DMD_galGal4',OKcheck=F)

res_10k_dmdE15 <- makeRFD(grE15_dmd,bs=10000,lr=20,na2zero=F,bin=T,datatype='OKseq',export=F,saveRData=F,retur=T,outname='DT40_OK_DMD_galGal4',OKcheck=F)


cor.rfd.test2(sm.wtdmd.10k[[4]], res_10k_dmdE14[[4]],bs0=10000,binned=T)
#      rho
#0.8118004
cor.rfd.test2(sm.wtdmd.10k[[4]], res_10k_dmdE15[[4]],bs0=10000,binned=T)
#      rho
#0.8223484

cor.rfd.test2(res_10k_dmd[[4]], res_10k_dmdE14[[4]],bs0=10000,binned=T)
# 0.8818358
cor.rfd.test2(res_10k_dmd[[4]], res_10k_dmdE15[[4]],bs0=10000,binned=T)
# 0.8572687
cor.rfd.test2(res_10k_dmdE14[[4]], res_10k_dmdE15[[4]],bs0=10000,binned=T)
# 0.87142

res_10k_ccser1E14 <- makeRFD(grE14_ccser1,bs=10000,lr=20,na2zero=F,bin=T,datatype='OKseq',export=F,saveRData=F,retur=T,outname='DT40_OK_CCSER1_galGal4',OKcheck=F)

res_10k_ccser1E15 <- makeRFD(grE15_ccser1,bs=10000,lr=20,na2zero=F,bin=T,datatype='OKseq',export=F,saveRData=F,retur=T,outname='DT40_OK_CCSER1_galGal4',OKcheck=F)


cor.rfd.test2(sm.wtccser1.10k[[4]], res_10k_ccser1E14[[4]],bs0=10000,binned=T)
#      rho
#0.6445024
cor.rfd.test2(sm.wtccser1.10k[[4]], res_10k_ccser1E15[[4]],bs0=10000,binned=T)
#      rho
#0.6525184

cor.rfd.test2(res_10k_ccser1[[4]], res_10k_ccser1E14[[4]],bs0=10000,binned=T)
# .9661972
cor.rfd.test2(res_10k_ccser1[[4]], res_10k_ccser1E15[[4]],bs0=10000,binned=T)
# 0.9765258
cor.rfd.test2(res_10k_ccser1E14[[4]], res_10k_ccser1E15[[4]],bs0=10000,binned=T)
# 0.9672703

### exporting session info
library("devtools")
library(magrittr)
session_info() %>% capture.output(file="session_info.txt")
