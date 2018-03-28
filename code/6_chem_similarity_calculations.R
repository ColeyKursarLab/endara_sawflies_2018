library(foreach)
library(doParallel)
source("./code/create_pairwiseComps_sampsByComps.R")
source("./code/chem_similarity_function.R")

####### Download 'veiw raw spectra' from GNPS results page and use make_pairwisecomps function from 'create_pairwiseComps_sampsByComps.R' to creating pairwise compound similarity matrix #######
pairwise.comps <- make_pairwisecomps("./results/GNPS_view_raw_spectra/")

####### Load filled compound table created in 'fill_compounds.py' and alter format using make_sampsByCompound function from 'create_pairwiseComps_sampsByComps.R' #######
sampsByCompounds <- make_sampsByCompounds("./results/filled_compounds.csv", by_species = FALSE)
sampsByCompoundsLog <- log(sampsByCompounds) # take log of TIC values for each compound
sampsByCompoundsLog[sampsByCompoundsLog <= 0] <- 0 # change any negative values to 0
write.csv(t(sampsByCompounds), "./data/filled_samps_by_comps_2017_11_20.csv")

####### Create list of phenolics and saponins based on mz and rt #######
# Load table listing all features associated with each compound created in 'run_compound_table_code.py'
compound_feature <- read.csv("./results/compound_feature_table.csv")
# Select most abundant feature associated with each compound
compound_max_feature <- aggregate(compound_feature$TIC, by=list(Category=compound_feature$compound_number),FUN=max)
names(compound_max_feature) <- c("compound_number","TIC_max")
compound_max_feature_2 <- merge(compound_max_feature, compound_feature, by.x=c("compound_number","TIC_max"), by.y=c("compound_number","TIC"))
# saponins approximated as greater than 580 Da and eluting after 18 minutes
sap_compounds <- compound_max_feature_2[compound_max_feature_2$mz >=580 & compound_max_feature_2$rt >= 18 & compound_max_feature_2$compound_number %in% names(pairwise.comps), "compound_number"]
# everything else with phenolics
phen_compounds <- compound_max_feature_2[!compound_max_feature_2$compound %in% sap_compounds & compound_max_feature_2$compound_number %in% names(pairwise.comps), "compound_number"]

###### Create sample by compound matrix for each compound class #######
# with raw TIC values
sampsByCompoundsSap <- sampsByCompounds[,names(sampsByCompounds) %in% as.character(sap_compounds), drop=FALSE]
sampsByCompoundsPhen <- sampsByCompounds[,names(sampsByCompounds) %in% as.character(phen_compounds), drop=FALSE]
# and with log of TIC values
sampsByCompoundsSapLOG <- sampsByCompoundsLog[,names(sampsByCompoundsLog) %in% as.character(sap_compounds), drop=FALSE]
sampsByCompoundsPhenLOG <- sampsByCompoundsLog[,names(sampsByCompoundsLog) %in% as.character(phen_compounds), drop=FALSE]
# standardize by sample so that total compound abundance in each sample sums to 1.0
sampsCompsStandSap <- standardizeByRow(sampsByCompoundsSapLOG)
sampsCompsStandPhen <- standardizeByRow(sampsByCompoundsPhenLOG)

####### Detect available cores and set up environment for parallelization #######
cores = detectCores()
cl <- makeCluster(cores[1]) #not to overload your computer
registerDoParallel(cl)

####### Run similarity calculation for saponins using chemical_similarity_single function from 'chemical_similarity_function.R' code #######
similarity_matrix_sap <- foreach(i = 1:nrow(sampsCompsStandSap), .combine = rbind) %:% foreach(j = 1:nrow(sampsCompsStandSap)) %dopar% {
  if(i <= j) {
  chemical_similarity_single(row.names(sampsCompsStandSap)[i], row.names(sampsCompsStandSap)[j], sampsCompsStandSap, pairwise.comps)
  }
  else NA
}
# change format and fill in rest of pairwise sample similarity matrix
similarity_matrix_sap <- as.data.frame(similarity_matrix_sap)
names(similarity_matrix_sap) <- row.names(sampsCompsStandSap)
row.names(similarity_matrix_sap) <- row.names(sampsCompsStandSap)
for(i in 1:ncol(similarity_matrix_sap)) similarity_matrix_sap[,i] <- unlist(similarity_matrix_sap[,i])
for(i in 1:nrow(similarity_matrix_sap)) {
  for(j in i:nrow(similarity_matrix_sap)) {
    similarity_matrix_sap[j,i] <- similarity_matrix_sap[i,j]
  }
}
# write saponin similarity to csv
write.csv(similarity_matrix_sap, "./results/saponin_similarity.csv")


####### Repeat similarity calculation for phenolics #######
similarity_matrix_phen <- foreach(i = 1:nrow(sampsCompsStandPhen), .combine = rbind) %:% foreach(j = 1:nrow(sampsCompsStandPhen)) %dopar% {
  if(i <= j) {
    chemical_similarity_single(row.names(sampsCompsStandPhen)[i], row.names(sampsCompsStandPhen)[j], sampsCompsStandPhen, pairwise.comps)
  }
  else NA
}
similarity_matrix_phen <- as.data.frame(similarity_matrix_phen)
names(similarity_matrix_phen) <- row.names(sampsCompsStandPhen)
row.names(similarity_matrix_phen) <- row.names(sampsCompsStandPhen)
for(i in 1:ncol(similarity_matrix_phen)) similarity_matrix_phen[,i] <- unlist(similarity_matrix_phen[,i])
for(i in 1:nrow(similarity_matrix_phen)) {
  for(j in i:nrow(similarity_matrix_phen)) {
    similarity_matrix_phen[j,i] <- similarity_matrix_phen[i,j]
  }
}
write.csv(similarity_matrix_phen, "./results/phenolic_similarity.csv")

####### Load data on tyrosine content and extraction weights #######
tyr.pct <- read.csv("./data/tyrosine_content.csv")
extr.pct <- read.csv("./data/extraction_weights.csv")

####### Combine phenolic and saponin similarity with tyrosine data #######
# calculate total investment in phenolics and saponins based on sum of TIC for all compounds in each class 
comp.class.pcts <- data.frame("sample" = row.names(sampsByCompounds), 
                              "phenSumTIC" = sapply(1:nrow(sampsByCompoundsPhen), function(x) sum(sampsByCompoundsPhen[x,])), 
                              "sapSumTIC" = sapply(1:nrow(sampsByCompoundsSap), function(x) sum(sampsByCompoundsSap[x,])), 
                              "species_code" = sapply(1:nrow(sampsByCompoundsPhen), 
                                                      function(x) unlist(strsplit(row.names(sampsByCompoundsPhen)[x], split = "_"))[1]),
                              stringsAsFactors=FALSE)
# take log of TIC for each class
comp.class.pcts$logPhen <- log(comp.class.pcts$phenSumTIC)
comp.class.pcts$logSap <- log(comp.class.pcts$sapSumTIC)
# calculate percent investment in saponins/phenolics
comp.class.pcts$phenTICpct <- comp.class.pcts$phenSumTIC / (comp.class.pcts$phenSumTIC + comp.class.pcts$sapSumTIC)
comp.class.pcts$sapTICpct <- comp.class.pcts$sapSumTIC / (comp.class.pcts$phenSumTIC + comp.class.pcts$sapSumTIC)
# get intermediate polarity (phenolics+saponins) extraction weight and tyrosine % dry weight
comp.class.pcts <- join(comp.class.pcts, extr.pct, by="species_code", type="left", match="all")
comp.class.pcts <- join(comp.class.pcts, tyr.pct, by="species_code", type="left", match="all")
comp.class.pcts$percent_tyrosine[is.na(comp.class.pcts$percent_tyrosine)] <- 0
# calculate tyrosine as % of total secondary metabolite extractions
comp.class.pcts$tyr.final.pct <- comp.class.pcts$percent_tyrosine / (comp.class.pcts$percent_extracted * 0.01 + comp.class.pcts$percent_tyrosine)
# calculate phenolics+saponins as % of total secondary metabolite extractions
comp.class.pcts$phensap.final.pct <- comp.class.pcts$percent_extracted * 0.01 / (comp.class.pcts$percent_extracted * 0.01 + comp.class.pcts$percent_tyrosine)
comp.class.pcts$phen.final.pct <- comp.class.pcts$phenTICpct * comp.class.pcts$phensap.final.pct
comp.class.pcts$sap.final.pct <- comp.class.pcts$sapTICpct * comp.class.pcts$phensap.final.pct

####### Create matrices of pairwise minimum phenolic, saponin, and tyrosine investment #######
pairwise.phen.min <- outer(comp.class.pcts$phen.final.pct, comp.class.pcts$phen.final.pct, FUN = function(X,Y) sapply(1:length(X), function(i) min(X[i],Y[i])))
pairwise.sap.min <- outer(comp.class.pcts$sap.final.pct, comp.class.pcts$sap.final.pct, FUN = function(X,Y) sapply(1:length(X), function(i) min(X[i],Y[i])))
pairwise.tyr.min <- outer(comp.class.pcts$tyr.final.pct, comp.class.pcts$tyr.final.pct, FUN = function(X,Y) sapply(1:length(X), function(i) min(X[i],Y[i])))

####### Combine matrices, using minimum investment in each compound class as weights
pairwise.spp <- similarity_matrix_phen*pairwise.phen.min + similarity_matrix_sap*pairwise.sap.min + pairwise.tyr.min
write.csv(pairwise.spp, "./results/chemical_similarity_matrix.csv")
