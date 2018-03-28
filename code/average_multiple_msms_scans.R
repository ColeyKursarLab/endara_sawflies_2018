library(mzR)
library(MSnbase)

# function that extracts mz/rt/TIC data for peaks from mzXML files. Requires vector containing full or relative filepaths for .mzXML files of interest and MSlevel to extract data for
extract.peakdata <- function(files, MSlevel) {
  curr_header <- header(openMSfile(files[1]))
  msms_scans <- curr_header[curr_header$peaksCount >= 5 & curr_header$totIonCurrent >= 2000 & curr_header$basePeakIntensity >=1000 & curr_header$msLevel %in% MSlevel,]
  if(nrow(msms_scans) > 0) msms_scans$file_idx <- 1
  else msms_scans$file_idx <- numeric(0)
  if(length(files) > 1) {
  for(j in 2:length(files)) {
    curr_header <- header(openMSfile(files[j]))
    curr_header$file_idx <- j
    msms_scans <- rbind(msms_scans, curr_header[curr_header$peaksCount >= 5 & curr_header$totIonCurrent >= 2000 & curr_header$basePeakIntensity >=1000 & curr_header$msLevel %in% MSlevel,])
  }
  }
  return(msms_scans)
}

# Function finds weighted average of multiple msms scans. Requires msms_scans dataframe produced by extract.peakdata, vector containing file paths for MSMS .mzXML files (IN SAME ORDER AS USED FOR extract.peakdata function) as well as mz and rt of peak of interest.
avg.msms.spec <- function(files, peak.data, rt, mz) {
  scans_to_merge <- peak.data[abs(peak.data$retentionTime - rt*60) <= 20 & abs(peak.data$precursorMZ - mz) <= 0.01, c("acquisitionNum", "file_idx"), drop = FALSE]
  if(nrow(scans_to_merge) < 1) return(NULL)
  spec1 <- lapply(1:nrow(scans_to_merge), function(j) mzR::peaks(openMSfile(files[scans_to_merge$file_idx[j]]),scans_to_merge$acquisitionNum[j]))
  spec2 <- do.call(rbind, spec1)
  spec3 <- spec2[order(spec2[,2], decreasing = T),]
  row.names(spec3) <- 1:nrow(spec3)
  
  rows_done = numeric()
  spec4 <- array(numeric(), dim=c(0,2))
  while(length(setdiff(1:nrow(spec3), rows_done)) > 0) {
    i = min(setdiff(1:nrow(spec3), rows_done))
    to_merge = setdiff(which(abs(spec3[,1] - spec3[i,1]) <= 0.25), rows_done)
    
    mz = weighted.mean(spec3[to_merge, 1], spec3[to_merge, 2])
    TIC = sum(spec3[to_merge, 2])
    
    spec4 <- rbind(spec4, c(mz, TIC))
    
    rows_done <- c(rows_done, to_merge)
  }
  spec5 <- spec4[order(spec4[,1]),]
  TICsum <- sum(spec5[,2])
  spec5 <- spec5[spec5[,2] >= TICsum/1000,,drop=FALSE]
  return(spec5)
}