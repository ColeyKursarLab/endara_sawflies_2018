library(reshape2)
library(igraph)
library(splitstackshape)

# function takes network from GNPS and makes compound x compound similarity matrix using cosine scores and calculating through-network similarity (default), can prevent this by using through_network = FALSE
# requires path to folder that come from 'View Raw Spectra' in GNPS
make_pairwisecomps <- function(network_filepath, through_network = TRUE) {
# load .pairsinfo file from networkedges_selfloop folder and raw spectra file from main folder
  if(is.character(network_filepath)) {
  network_edges_path <- list.files(paste(network_filepath, "/networkedges_selfloop", sep = ""), full.names = TRUE)
  raw_spectra_path <- list.files(network_filepath, full.names = TRUE)[grep(".tsv", list.files(network_filepath))]
  network_clean <- read.delim(file = network_edges_path, header=TRUE, sep="\t")
  compound_cluster_match <- read.delim(file = raw_spectra_path,header=TRUE, sep="\t")[,c("cluster.index","ScanNumber")]
  names(network_clean) <- c("CLUSTERID1","CLUSTERID2","DeltaMZ","V4","Cosine","OtherScore","V7")

  # convert network from clusterid to scan number (same as compound number)
  names(compound_cluster_match) <- c("CLUSTERID1","NODE1")
  compound_network <- merge(network_clean, compound_cluster_match, by="CLUSTERID1",all.x=TRUE)
  names(compound_cluster_match) <- c("CLUSTERID2","NODE2")
  compound_network <- merge(compound_network, compound_cluster_match, by="CLUSTERID2",all.x=TRUE)
  compound_network <- compound_network[,c("NODE1","NODE2","DeltaMZ","Cosine")]

  compound_network$Cosine <- as.numeric(as.character(compound_network$Cosine))
  } else {
   compound_network <- network_filepath 
  }
  
  # make changes to structure of network, calculate inverse cosine score for through-network calculations
  compound_network$inv_cos_score <- 1/compound_network$Cosine
  nodes <- unique(c(compound_network$NODE1, compound_network$NODE2))
  
  # compute shortest paths using inverse of cosine, then taking the inverse again
  g2 <- graph_from_data_frame(d=compound_network, vertices=data.frame(nodes), directed=FALSE)
  if(through_network == TRUE) {
  shortest_paths_list_inverse <- lapply(1:length(nodes), function(x) shortest.paths(g2,v=as.character(nodes[x]),weights=E(g2)$inv_cos_score))
  shortest_paths_df_inverse <- do.call(rbind.data.frame, shortest_paths_list_inverse)
  short_paths_inv_inv <- 1/shortest_paths_df_inverse
  short_paths_inv_inv[short_paths_inv_inv == Inf] <- 1
  pairwise.comps <- short_paths_inv_inv
  }
  else{
    pairwise.comps <- as.data.frame(matrix(0, ncol = length(nodes), nrow = length(nodes)))
    names(pairwise.comps) <- as.character(sort(nodes))
    row.names(pairwise.comps) <- as.character(sort(nodes))
    for(i in 1:nrow(compound_network)) {
      pairwise.comps[row.names(pairwise.comps) == as.character(compound_network$NODE1[i]),
                     names(pairwise.comps) == as.character(compound_network$NODE2[i])] <- compound_network$Cosine[i]
      pairwise.comps[row.names(pairwise.comps) == as.character(compound_network$NODE2[i]),
                     names(pairwise.comps) == as.character(compound_network$NODE1[i])] <- compound_network$Cosine[i]
    }
    for(i in 1:nrow(pairwise.comps)) {
      pairwise.comps[i,i] <- 1
    }
    }
  
  return(pairwise.comps) }

# function to create sampsByCompounds table
# requires compound tic table, which was created in python
# argument by_species = TRUE specifies that samples should be grouped by species. If FALSE (default), each sample will remain separate.
# if you want to remove some samples or species, specify them in samps_to_remove

make_sampsByCompounds <- function(compound_tic_table_path, samps_to_remove = character(0), by_species = FALSE) {
  if(is.character(compound_tic_table_path)) compound_tic_table <- read.csv(compound_tic_table_path)
  else compound_tic_table <- compound_tic_table_path
  sampsByCompounds <- dcast(compound_tic_table, compound_number  ~ species_sample, value.var = "TIC")

  if(by_species == TRUE) {
  compound_tic_table_species <- data.frame(compound_tic_table,cSplit(compound_tic_table, 'species_sample', sep="_", type.convert=FALSE))
  compound_tic_table_species_1 <- aggregate(compound_tic_table_species$TIC, by=list(Category=compound_tic_table_species$compound_number,compound_tic_table_species$species_sample_1),FUN=mean)
  names(compound_tic_table_species_1) <- c( "compound_number","species","TIC")
  sampsBySpecies <- dcast(compound_tic_table_species_1, compound_number  ~ species, value.var = "TIC", fun.aggregate = sum)
  sampsByCompounds <-sampsBySpecies }

  sampsByCompounds <- as.data.frame(sampsByCompounds)
  sampsByCompounds <- t(sampsByCompounds)
  sampsByCompounds[is.na(sampsByCompounds)] <- 0
  sampsByCompounds <- as.data.frame(sampsByCompounds)
  row.names(sampsByCompounds)
  names(sampsByCompounds) <- sampsByCompounds["compound_number",]
  sampsByCompounds <- sampsByCompounds[!row.names(sampsByCompounds)%in%c("compound_number",samps_to_remove),]

  return(sampsByCompounds) }

standardizeByRow <- function(dataframe) {
  for(i in 1:nrow(dataframe)) {
    dataframe[i, ] <- dataframe[i, ] / sum(dataframe[i, ])
  }
  return(dataframe)
}
