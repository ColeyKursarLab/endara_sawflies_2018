library(splitstackshape)
source("./code/average_multiple_msms_scans.R")

####### Load filled compound table created in 'fill_compounds.py' #######
compound_table <- read.csv("./results/filled_compounds.csv", stringsAsFactors=FALSE)
# make row for species
compound_table$species <-  sapply(1:nrow(compound_table), function(x) unlist(strsplit(compound_table$species_sample[x], split="_"))[1])
# aggregate compound table by species
compound_table_2 <- aggregate(compound_table$TIC, by=list(compound_table$species, compound_table$compound_number), FUN=mean)
names(compound_table_2) <- c("species","compound_number","TIC")

####### Load table listing all features associated with each compound created in 'run_compound_table_code.py' #######
compound_feature <- read.csv("./results/compound_feature_table.csv")
# Select most abundant feature associated with each compound
compound_max_feature <- aggregate(compound_feature$TIC, by=list(Category=compound_feature$compound_number),FUN=max)
names(compound_max_feature) <- c("compound_number","TIC_max")
compound_max_feature_2 <- merge(compound_max_feature, compound_feature, by.x=c("compound_number","TIC_max"), by.y=c("compound_number","TIC"))
compound_table_3 <- merge(compound_table_2, compound_max_feature_2, by="compound_number", all.x = FALSE)
head(compound_table_3)


####### Search DDA and MSMS files for msms spectra of each compound #######
species <- unique(compound_table_3$species)
msms_spec <- vector("list", max(compound_table_3$compound))

for(i in 1:length(species)) {
  print(as.character(species[i]))

  file.paths <- list.files(paste("./data/DDA_and_MSMS_files/", species[i], sep=""), full.names=TRUE) # list of MSMS and DDA files from current species
  peakdata <- extract.peakdata(file.paths, MSlevel = 2) # use function from 'average_multiple_msms_scans.R' to get precursor ion information from all relevant files
  spec_i <- compound_table_3[compound_table_3$species == species[i], ] # compounds to search for in these files
  # for each compound, search for MSMS spectra
  for(k in 1:nrow(spec_i)) {
    rt_k <- spec_i$rt[k]
    mz_k <- spec_i$mz[k]
    current_spec <- avg.msms.spec(file.paths, peakdata, rt = rt_k, mz = mz_k) # average all MSMS spectra with this mz/rt combo
    if(is.null(current_spec)) next # if current_spec is NULL, a spectrum was not collected for this compound
    if(nrow(current_spec) < 5) next # filter out if spectrum contains fewer than 5 peaks
    if(max(msms_spec[[spec_i$compound[k]]][,2]) > max(current_spec[,2])) next # if spectrum of same compound was already collected from different species, keep the one with the highest signal
    msms_spec[[spec_i$compound[k]]] <- current_spec
  }
}

####### For compounds that don't have associated spectra with their most abundant feature, look for spectra associated with other features #######
no_spec <- which(sapply(1:length(msms_spec), function(x) is.null(msms_spec[[x]])))

for(i in no_spec) {
  print(i)
  compound_table_i <- compound_table_3[compound_table_3$compound == i,] # subset compound table to only include current compound
  if(nrow(compound_table_i) == 0) next
  features_i <- compound_feature[compound_feature$compound_number == i,] # get features associated with current compound
  features_i$relTIC <- features_i$TIC / max(features_i$TIC) # calculate abundance of each feature relative to most abundant feature
  features_i <- features_i[features_i$relTIC > 0.33 & features_i$relTIC != 1,] # keep features at least 1/3 as abundant as most abundant feature
  if(nrow(features_i) == 0) next 
  features_i <- features_i[order(features_i$relTIC, decreasing=TRUE),]
  compound_matched = FALSE
  # search DDA and MSMS files associated with these features, starting with the second most abundant feature, and looping through each species in which the compound is found
  for(curr_species in compound_table_i$species) {
    file.paths <- list.files(paste("./data/DDA_and_MSMS_files/", curr_species, sep=""), full.names=TRUE)
    peakdata <- extract.peakdata(file.paths, MSlevel = 2)
    for(k in 1:nrow(features_i)) {
      rt_k <- features_i$rt[k]
      mz_k <- features_i$mz[k]
      current_spec <- avg.msms.spec(file.paths, peakdata, rt = rt_k, mz = mz_k)
      if(is.null(current_spec)) next
      if(nrow(current_spec) < 5) next
      msms_spec[[features_i$compound_number[k]]] <- current_spec
      compound_matched = TRUE
      break
    }
    if(compound_matched == TRUE) break # if msms spectrum is found, continue on to next compound. If not, search through remaining species
  }
}

####### Save list of MSMS spectra as .rds object #######
saveRDS(msms_spec, "./results/msms_spectra_list.rds")

####### Create vector with compound numbers that have MSMS spectra #######
comps_with_spec <- which(sapply(1:length(msms_spec), function(x) !is.null(msms_spec[[x]])))
# Get mz and rt data for these compounds 
compound_max_feature_3 <- compound_max_feature_2[compound_max_feature_2$compound %in% comps_with_spec, ]

####### Write .mgf containing MSMS spectrum for each compound #######
write_txt_file <- file("./results/all_msms_spectra.mgf")
writeLines(
  unlist(lapply(1:nrow(compound_max_feature_3), function(i)
    c("BEGIN IONS",paste("PEPMASS=",compound_max_feature_3$mz[i],sep=""),"CHARGE=1-",paste("SCANS=",compound_max_feature_3$compound[i],sep=""),sapply(1:nrow(msms_spec[[compound_max_feature_3$compound[i]]]), function(x) paste(msms_spec[[compound_max_feature_3$compound[i]]][x,1],msms_spec[[compound_max_feature_3$compound[i]]][x,2],sep="\t")),"END IONS")))
  , write_txt_file)
close(write_txt_file)

####### Upload .mgf to GNPS (gnps.ucsd.edu) for molecular networking #######