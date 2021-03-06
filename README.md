# About
This repository contains R and Python code used for the chemical similarity components of the 2018 Endara et al. paper entitled 
'Co-diversification of leaf-feeding sawflies (Hymenoptera) and their Inga (Fabaceae) hosts in the Neotropics'.

# Dependencies
R code operates using R version 3.4.3 and requires packages XCMS, multtest, snow, muma, pvclust, data.table, vegan, splitstackshape, CAMERA, foreach, doParallel, gdata, reshape2, and igraph.
Python code operates using Python 2.7.11 and uses the module imp.

# Structure
The six numbered files in the code folder should be run in sequential order and will complete analysis from peak picking from raw chromatograms to calculating a pairwise chemical similarity score for all samples. 

# Data requirements
Sample .mzXML files will be uploaded to MetaboLights. In order to run this code from the beginning, single MS files should be grouped 
into folders based on species. Within each species folder, data files should be placed in a folder called 'Sample', and associated
blanks should be placed in a folder called 'Blank'. These folders should be saved in a folder entitled 'mzxml_files' within the 'data' 
folder. 
MSMS files should again be grouped into folders by species saved in a folder entitled 'DDA_and_MSMS_files' within the 'data' folder.

Other required data files include a list of extraction weights (as % of dry weight), associations between samples and a blank run on
the same day, tyrosine content (as % of dry weight), and a list of adducts for CAMERA. Sample files are included in the 'data' folder
for species M4 and T69. Some of these files might not be included if not relevant to the analysis. For example, if analyzing species
that do not overexpress tyrosine, the tyrosine content file could be omitted, and the related lines of code would be skipped over.

Intermediate data files are stored in the 'results' folder and are derived from completing the analysis using only
species M4 and T69.
