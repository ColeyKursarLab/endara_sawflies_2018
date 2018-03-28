####### Load required packages #######

source("http://bioconductor.org/biocLite.R")
biocLite("xcms")
biocLite("multtest")
library(multtest)
library(xcms)
library(snow)
library(muma)
library(pvclust)
library(data.table)
library(vegan)
library(splitstackshape)
library(CAMERA)

### ARE THESE LINES ACTUALLY NECESSARY?????? ####
imports = parent.env(getNamespace("CAMERA"))
unlockBinding("groups", imports)
imports[["groups"]] = xcms::groups
lockBinding("groups", imports)

####### Load CAMERA file with list of possible adducts #######
rule_mod <- read.csv("./data/negative_mode_adducts_for_CAMERA.csv", header= TRUE)

###### Prepare data frame to fill with peak data #######
all_features_long <- data.frame("sample"=NULL,"Species_code_sample"=NULL,"RT"=NULL,"MZ"=NULL,"PC_ID"=NULL,"TIC"=NULL)

species <- list.files("./data/mzxml_files/")

####### Loop through directories containing .mzXML files. Samples and their associated blanks should be in separate folders named "Sample" and "Blank" #######
for(j in 1:length(species)) {
  print(species[j])

  files <- list.files(paste("./data/mzxml_files/", species[j], sep=""), recursive = TRUE, full.names = TRUE)
  sample_groups <- sapply(1:length(files), function(x) unlist(strsplit(files[x], split = "/"))[5])
  sample_names <- sapply(1:length(files), function(x) sub(unlist(strsplit(files[x], split = "/"))[6], pattern = ".mzXML", replacement = "", fixed = TRUE))
  
  # set xcms peak picking and grouping parameters
  cwparam <- CentWaveParam(ppm=15, peakwidth=c(4,12), snthresh=5, prefilter=c(10,500))
  dparam1 <- PeakDensityParam(sampleGroups = sample_groups, bw=10, binSize=0.05, minSamples=1, minFraction = 0.01)
  oparam <- ObiwarpParam(binSize=1)
  dparam2 <- PeakDensityParam(sampleGroups = sample_groups, bw = 3, binSize = 0.025, minSamples = 1, minFraction = 0.01)
  fpparam <- FillChromPeaksParam(expandMz = 0, expandRt = 0, ppm = 0)
  
  # read in data from .mzXML files
  raw_data <- readMSData(files, msLevel. = 1, mode="onDisk")
  row.names(raw_data@phenoData) <- sample_names
  # remove peaks before or after retention time cutoff (in seconds)
  minrt <- 60 # remove peaks before 1 minute
  maxrt <- 32 * 60 # remove peaks after 32 minutes
  fData(raw_data) <- fData(raw_data)[fData(raw_data)$retentionTime > minrt & fData(raw_data)$retentionTime < maxrt, ]
  
  
  xcmsexp <- findChromPeaks(object = raw_data, param = cwparam)   # peak identification
  xcmsexp <- groupChromPeaks(object = xcmsexp, param = dparam1)   # peak grouping across samples 1
  xcmsadjustrtime <- adjustRtime(object = xcmsexp, param = oparam)   # retention time correction
  xcmsexp <- groupChromPeaks(object = xcmsexp, param = dparam2)   # peak grouping across samples 2
  xcmsexp <- fillChromPeaks(xcmsexp, fpparam)   # fill peaks in all samples
  xset <- as(xcmsexp, "xcmsSet")   # convert xcmsExp object to xcmsSet object
  sampclass(xset) <- sample_groups
  
  ####### CAMERA annotation and grouping ########
  xset1 <- xsAnnotate(xs=xset,polarity="negative")   # create xsAnnotate object
  xset2 <- groupFWHM(xset1, perfwhm=0.7)   # group peaks based on overlap at FWHM
  xset3 <- findIsotopes(xset2, ppm=20, mzabs=0.015,intval="intb")   # annotate isotope peaks
  xset4 <- groupCorr(xset3,cor_eic_th=0.1, pval=1.0, graphMethod="lpc", calcIso = TRUE, calcCiS = TRUE, calcCaS = ifelse(
    sum(sampclass(xset) == "sample") > 3, TRUE, FALSE))   # regroup peaks based on correlation of peak intensities
  xsetFA <- findAdducts(xset4, polarity="negative", rules = rule_mod)   # annotate adducts
  xset5<- getPeaklist(xsetFA)   # return peak list
  
  ####### Delete peaks found in blank and rows where ratio of blank to sample is greater than 0.33 ####### 
  xset5[is.na(xset5)] <- 0
  if(length(which(sampclass(xset) == "Blank"))>1) {
    xset5$Blank_Average <- rowMeans(xset5[,c(7+length(unique(sampclass(xset)))+which(sampclass(xset) == "Blank"))],na.rm=T)
  } else {
    xset5$Blank_Average <- xset5[,c(7+length(unique(sampclass(xset)))+which(sampclass(xset) == "Blank"))]
  }
  xset5$TIC_Average <- rowMeans(xset5[,c(7+length(unique(sampclass(xset)))+which(sampclass(xset) == "Sample")), drop=FALSE],na.rm=T)
  xset5$TIC_ratio <- xset5$Blank_Average/xset5$TIC_Average
  xset6 <- xset5[xset5$Blank==0 & xset5$TIC_ratio < 0.33, ]
  
  ####### remove biochanin internal standard peak and formate dimer peaks #######
  if(nrow(xset6[which(abs(xset6$rt - 1350) <= 50 & abs(xset6$mz-268.036) <=0.1),]) >0) {xset6 <- xset6[-which(abs(xset6$rt - 1350) <= 50 & abs(xset6$mz-268.036) <=0.1),]} 
  if(nrow(xset6[which(abs(xset6$rt - 1350) <= 50 & abs(xset6$mz-269.042) <=0.1),]) >0) xset6 <- xset6[-which(abs(xset6$rt - 1350) <= 50 & abs(xset6$mz-269.042) <=0.1),]  
  if(nrow(xset6[which(abs(xset6$rt - 1350) <= 50 & abs(xset6$mz-283.062) <=0.1),]) >0) xset6 <- xset6[-which(abs(xset6$rt - 1350) <= 50 & abs(xset6$mz-283.062) <=0.1),]  
  if(nrow(xset6[which(abs(xset6$rt - 1350) <= 50 & abs(xset6$mz-284.062) <=0.1),]) >0) xset6 <-  xset6[-which(abs(xset6$rt - 1350) <= 50 & abs(xset6$mz-284.062) <=0.1),]  
  if(nrow(xset6[which(abs(xset6$rt - 1350) <= 50 & abs(xset6$mz-285.062) <=0.1),]) >0) xset6 <- xset6[-which(abs(xset6$rt - 1350) <= 50 & abs(xset6$mz-285.062) <=0.1),] 
  if(nrow(xset6[which(abs(xset6$mz-91.008) <=0.01),]) > 0) xset6 <- xset6[-which(abs(xset6$mz-91.008) <=0.01),] 
  
  ####### Add columns and change format of from wide to long #######
  xset6$rt_in_min<- (xset6$rt)/60  
  xset6$mz_round <- round((xset6$mz),4)
  xset6$rt_round <- round((xset6$rt_in_min),4)
  
  xset6_long <- melt(data=xset6, id.vars=c("mz_round","rt_round","pcgroup"),measure.vars=c(7+length(unique(sampclass(xset)))+which(sampclass(xset) == "Sample")))
  xset6_long$PC_ID <- sapply(1:nrow(xset6_long), function(x) paste(unlist(strsplit(as.character(xset6_long$variable)[x],"_"))[1], as.character(xset6_long$pcgroup)[x],sep="_"))
  xset6_long$sample <- sapply(1:nrow(xset6_long), function(x) paste("DN",unlist(strsplit(as.character(xset6_long$variable)[x],"_"))[2],sep=""))
  feature_table_long <- xset6_long[,c(7,4,2,1,6,5)]
  names(feature_table_long) <- c("sample","Species_code_sample","RT","MZ","PC_ID","TIC")
  
  ####### Append to peak table from previously processed samples #######
  all_features_long <- rbind(all_features_long, feature_table_long)
}

####### Write peak table to CSV #######
write.csv(all_features_long, "./results/all_peaks_long.csv", row.names=FALSE)


