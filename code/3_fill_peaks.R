library(xcms)

####### load feature table built in 'run_compound_table_code.py' #######
ms_features <- read.csv("./results/features_matched_across_samples.csv")
# make list of unique fatures, and average mz and rt for all times that feature was found
features <- data.frame("feature_number" = unique(ms_features$feature_number), mz = numeric(length(unique(ms_features$feature_number))), rt = numeric(length(unique(ms_features$feature_number))))
for(i in 1:nrow(features)) {
  features$mz[i] <- mean(ms_features[ms_features$feature_number == features$feature_number[i], "mz"])
  features$rt[i] <- mean(ms_features[ms_features$feature_number == features$feature_number[i], "rt"])
}

####### load file associating sample names blank run on same day #######
sample_blank <- read.csv("./data/sample_blank_association.csv")

####### set mz and retention time tolerance to use for finding features #######
mz_tolerance <- 25 # mz tolerance (ppm)
rt_tolerance <- 30 # rt tolerance (seconds)
features$mzmin <- features$mz*(1-mz_tolerance/1000000)
features$mzmax <- features$mz*(1+mz_tolerance/1000000)
features$rtmin <- features$rt*60 - rt_tolerance
features$rtmax <- features$rt*60 + rt_tolerance
ms_features_matrix <- as.matrix(features[, c("mzmin","mzmax","rtmin","rtmax")])


####### Load .mzXML file for each file and search for all features using fillPeaks function from package xcms #######
species <- list.files("./data/mzxml_files/")
for(k in 1:length(species)) {
  samples <- list.files(paste("./data/mzxml_files/",species[k],"/Sample/",sep=""), full.names = TRUE)
    for(l in 1:length(samples)) {
      print(samples[l])
      
      # fill peaks for sample current sample
      sample_xcms <- xcmsRaw(samples[l], profstep = 1, profmethod = , includeMSn = FALSE, mslevel = 1)
      sample_peaks <- getPeaks(object = sample_xcms, peakrange = ms_features_matrix, step = 0.05)
      sample_peaks_1 <- as.data.frame(sample_peaks)
      sample_peaks_2 <- sample_peaks_1[sample_peaks_1$into >= 1000,]  # only keep if TIC above 1000
      sample_peaks_2$mass_matches <- character(nrow(sample_peaks_2))
      sample_peaks_2$feature_number <- numeric(nrow(sample_peaks_2))
      sample_peaks_2$feature_mz <- numeric(nrow(sample_peaks_2))
      sample_peaks_2$feature_rt <- numeric(nrow(sample_peaks_2))
      sample_peaks_2$actual_mz <- numeric(nrow(sample_peaks_2))
      sample_peaks_2$actual_rt <- numeric(nrow(sample_peaks_2))
      sample_peaks_2$TIC <- numeric(nrow(sample_peaks_2))
      
      # confirm that mass matches within defined tolerance
      for(i in 1:nrow(sample_peaks_2)) {
        # extract raw scans identified as containing current feature
        feature_number <- features[features$mzmin == sample_peaks_2$mzmin[i] & features$rtmin == sample_peaks_2$rtmin[i], "feature_number"]
        sample_peaks_2$feature_number[i] <- feature_number
        sample_peaks_2$feature_mz[i] <- features[features$feature_number == feature_number, "mz"]
        sample_peaks_2$feature_rt[i] <- features[features$feature_number == feature_number, "rt"]
        peakscan <- which.min(abs(sample_xcms@scantime - sample_peaks_2$rt[i]))
        test2 <- getScan(sample_xcms, peakscan, mzrange = c(sample_peaks_2$mzmin[i],sample_peaks_2$mzmax[i]))
        
        # if feature is not present in scan, there was a mismatch. Remove feature from list and continue to next feature.
        if(length(test2) == 0) {
          sample_peaks_2$mass_matches[i] <- "mismatch"
          next }
        
        # if mz is within defined tolerance, make sure feature is roughly peak-shaped and not just noise
        if(abs(sample_peaks_2$feature_mz[i] - test2[1,1]) * 1000000 / sample_peaks_2$feature_mz[i]  < mz_tolerance) {
          # extract vector of intensities of feature for 20 scans surrounding peak
          intensities <- sapply((peakscan-20):(peakscan+20), function(x) {
            tempscan = getScan(sample_xcms, x, mzrange = c(sample_peaks_2$mzmin[i],sample_peaks_2$mzmax[i]))
            ifelse(nrow(tempscan) == 0, 0, tempscan[1,2]) } )
          midpoint <- length(intensities)/2 + 0.5
          # midpoint of peak must be at least 1.2x as intense as edges of peak, and the correct mass must be present in the middle seven scans of the peak.
          if(intensities[midpoint] / intensities[1] > 1.2 & 
             intensities[midpoint] / intensities[length(intensities)] > 1.2 &
             sum((midpoint-3):(midpoint+3) %in% which(intensities == 0)) <= 1) {
          # if these criteria have been met, figure out where peak begins and ends and sum intensities
            diffs <- which(intensities == 0) - midpoint
          if(sum(diffs<0) == 0) first_scan = -20
          else first_scan <- diffs[diffs<0][sum(diffs<0)]
          if(sum(diffs>0) == 0) last_scan = 20
          else last_scan <- diffs[diffs>0][1]
            sample_peaks_2$mass_matches[i] <- "match"
          # get actual mz, rt, and TIC values for current peak
            sample_peaks_2$actual_mz[i] <- test2[1,1]
            sample_peaks_2$actual_rt[i] <- sample_peaks_2$rt[i]/60
            sample_peaks_2$TIC[i] <- sum(intensities[(midpoint+first_scan):(midpoint+last_scan)])}
          # if critera aren't met, peak is considered noise and is removed
          else sample_peaks_2$mass_matches[i] <- "noise" }
        else sample_peaks_2$mass_matches[i] <- "mismatch"
      }
sample_peaks_3 <- sample_peaks_2[sample_peaks_2$mass_matches == "match", c("feature_number","feature_mz","feature_rt", "into", "actual_mz", "actual_rt", "TIC")]

sample_name <- unlist(strsplit(samples[l], split = "/"))[6]
blank_name <- as.character(sample_blank[sample_blank$file_name == sample_name, "associated_blank"])

# use same process to check associated blank for all features
blank_xcms <- xcmsRaw(paste("./data/mzxml_files/", species[k], "/Blank/", blank_name, ".mzXML", sep = ""), profstep = 1, profmethod = , includeMSn = FALSE, mslevel = 1)
blank_findpeaks <- getPeaks(object = blank_xcms, peakrange = ms_features_matrix, step = 0.05)
blank_findpeaks_1 <- as.data.frame(blank_findpeaks)
blank_findpeaks_2 <- blank_findpeaks_1[blank_findpeaks_1$into > 100,]
blank_findpeaks_2$mass_matches <- character(nrow(blank_findpeaks_2))
blank_findpeaks_2$feature_number <- numeric(nrow(blank_findpeaks_2))
blank_findpeaks_2$feature_mz <- numeric(nrow(blank_findpeaks_2))
blank_findpeaks_2$feature_rt <- numeric(nrow(blank_findpeaks_2))
blank_findpeaks_2$TIC <- numeric(nrow(blank_findpeaks_2))

for(i in 1:nrow(blank_findpeaks_2)) {
  feature_number <- features[features$mzmin == blank_findpeaks_2$mzmin[i] & features$rtmin == blank_findpeaks_2$rtmin[i], "feature_number"]
  blank_findpeaks_2$feature_number[i] <- feature_number
  blank_findpeaks_2$feature_mz[i] <- features[features$feature_number == feature_number, "mz"]
  blank_findpeaks_2$feature_rt[i] <- features[features$feature_number == feature_number, "rt"]
  peakscan <- as.numeric(blank_xcms@scanindex[which.min(abs(blank_xcms@scantime - blank_findpeaks_2$rt[i]))])
  peakscan <- which.min(abs(blank_xcms@scantime - blank_findpeaks_2$rt[i]))
  test2 <- getScan(blank_xcms, peakscan, mzrange = c(blank_findpeaks_2$mzmin[i],blank_findpeaks_2$mzmax[i]))
  if(length(test2) == 0) {
    blank_findpeaks_2$mass_matches[i] <- "mismatch"
    next }
  if(abs(blank_findpeaks_2$feature_mz[i] - test2[1,1]) * 1000000 / blank_findpeaks_2$feature_mz[i]  < mz_tolerance) {
    intensities <- sapply((peakscan-20):(peakscan+20), function(x) {
      tempscan = getScan(blank_xcms, x, mzrange = c(blank_findpeaks_2$mzmin[i],blank_findpeaks_2$mzmax[i]))
      ifelse(nrow(tempscan) == 0, 0, tempscan[1,2]) } )
    midpoint <- length(intensities)/2 + 0.5
      diffs <- which(intensities == 0) - midpoint
      if(sum(diffs<0) == 0) first_scan = -20
      else first_scan <- diffs[diffs<0][sum(diffs<0)]
      if(sum(diffs>0) == 0) last_scan = 20
      else last_scan <- diffs[diffs>0][1]
      blank_findpeaks_2$mass_matches[i] <- "match"
      blank_findpeaks_2$TIC[i] <- sum(intensities[(midpoint+first_scan):(midpoint+last_scan)])}
  else blank_findpeaks_2$mass_matches[i] <- "mismatch"
}
blank_findpeaks_3 <- blank_findpeaks_2[blank_findpeaks_2$mass_matches == "match",c("feature_number","feature_mz","feature_rt", "TIC")]

# delete all peaks that are not at least 5x as abundant in sample as in blank.
sample_blank_peaks <- merge(sample_peaks_3, blank_findpeaks_3[,c("feature_number", "TIC")], by = "feature_number", all.x = TRUE, all.y = FALSE)
sample_blank_peaks$TIC.y[is.na(sample_blank_peaks$TIC.y)] <- 0
sample_peaks <- sample_blank_peaks[sample_blank_peaks$TIC.x / sample_blank_peaks$TIC.y > 5, ]

sample_peaks_1 <- data.frame("feature_number" = sample_peaks$feature_number, "TIC" = sample_peaks$TIC.x, "actual_mz" = sample_peaks$actual_mz, "actual_rt" = sample_peaks$actual_rt, "sample_name" = unlist(strsplit(sample_name, split = "[.]"))[1])

####### write results to CSV #######
if(k == 1 & l == 1) {
  write.table(sample_peaks_1, "./results/filled_features.csv", sep = ",", append = FALSE, row.names = FALSE, col.names = TRUE)
}
else {
write.table(sample_peaks_1, "./results/filled_features.csv", sep = ",", append = TRUE, row.names = FALSE, col.names = FALSE) }
    }}
