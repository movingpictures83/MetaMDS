


library(vegan)
library(ggplot2)
library(gridExtra)
library(grid)
library(lattice)
library(scales)
library(stringi)
library(MASS)
#library(stargazer)
library(reporttools)
library(epitools)
library(gdata)
library(car)
library(plyr)
library(dplyr)
library(data.table)
library(tibble)
library(psych)
library(tidyr)
library(janitor)
library(psych)
library(plotrix)
library(slopegraph)
library(Lock5Data)
library(tidyverse)
library(ggrepel)
library(RColorBrewer)
library(treemap)
library (treemapify)
library(ggraph)
library(igraph)



###The code, uses the following "INPUT" files:

# CARD_read_count_specific.tsv
# ResPipe_CARD-3.0.3.meta.tsv
# CARD_lateral_coverage_specific.tsv
# CARD_read_lengths_specific.tsv
# Metagenomics_Metadata.csv
# Infection_Data_For_Bayesian_Model.csv
# bracken_combined_reads.tsv

# to produce the corrected resistance gene counts, plus a matrix that links each resistance gene and antibiotic based on the "Confers_Resistance_to_Antibiotic" relationship ontology term.
    # in the matrix, 1 == the gene is associated with clear experimental evidence of elevated MIC for that antibiotic but the "Confers_Resistance_to_Antibiotic" relationship ontology term is missing
    # in the matrix, 2 == the gene is associated with demonstrably elevated MIC for that antibiotic and is known to confer or contribute to clinically relevant resistance to that antibiotic ("Confers_Resistance_to_Antibiotic" relationship ontology term is present)

#OUTPUT FILES:

# Corrected_Counts.csv
# AB_Matrix_1_or_2.csv

# AMR_DEF.csv
# AMR_ALL.csv

# And to produce the final dataset for the bayesian modelling (i.e. "OUTPUT" file)

# Dataset_For_Bayesian_Model.csv

# The code for the non-metric multidimensional scaling (NMDS) ordination method is also presented (see supplementary results for the validation of the pooling using 30-sample pools)


##################PRODUCE THE DATASET WITH CORRECTED GENE COUNTS#############################

input <- function(inputfile) {

#ARO_Counts <- readRDS("RDS/ARO_Counts.rds")
ARO_Counts <<- readRDS(inputfile)
}

run <- function() {}

output <- function(outputfile) {


#########################################################non-metric multidimensional scaling (NMDS) ############################

Intrg_3<-ARO_Counts


Intrg_4<-Intrg_3[, -grep("POOL", colnames(Intrg_3))]

patterns <- unique(substr(names(Intrg_4), 1, 2))
new <- sapply(patterns, function(xx) rowSums(Intrg_4[,grep(xx, names(Intrg_4)), drop=FALSE]))

new<-as.data.frame(new)
new$FC<-as.numeric(new$FC)
new$PK<-as.numeric(new$PK)
new$FCPK <- new$FC + new$PK
new$FC <- NULL
new$PK<-NULL
colnames(new)[colnames(new) == "FCPK"] <- "KENYA_Sum_Ind"
colnames(new)[colnames(new) == "NC"] <- "CAMBODIA_sum_Ind"
colnames(new)[colnames(new) == "UK"] <- "UK_Sum_Ind"

Intrg_5 <- Intrg_3[ , grepl( "POOL" , names( Intrg_3 ) ) ]
Intrg_6 <- merge(Intrg_5, new, by="row.names", all=TRUE)


#/transpose the matrix

Intrg_7 <- Intrg_6[,-1]##Intrg_7 <- Intrg_6[,-1]
rownames(Intrg_7) <- Intrg_6[,1]

counts_lat_cov<-t(Intrg_7)
#fix(counts_lat_cov)


counts_lat_cov <- counts_lat_cov[ order(row.names(counts_lat_cov)), ]

classification=c(rep("Cambodia",3),rep("Kenya",3),rep("UK",3))

counts_lat_cov<-as.data.frame(counts_lat_cov)

total_col<-apply(counts_lat_cov[,], 1, sum)

relative_abundance = lapply(counts_lat_cov[,], function(x) {
  (x / total_col)*100
})

relative_abundance<-as.data.frame(relative_abundance)


set.seed(10)
nmds1 <- metaMDS(counts_lat_cov, binary = F,k = 2,try=1000)


set.seed(10)
nmds2 <- metaMDS(relative_abundance, binary = F,k = 2,try=1000)


}


