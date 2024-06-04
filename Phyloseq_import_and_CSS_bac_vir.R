library(phyloseq)
library(ggplot2)      # graphics
library(readxl)       # necessary to import the data from Excel file
library(dplyr)        # filter and reformat data frames
library(stringr)
library(vegan)
library(ggpubr)
library(metagenomeSeq)
library(tidyr)
library(RColorBrewer)
library(reshape2)


#Set working directory to r-script directory location - will not be run when script is run through source command

#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

######################################Construct Phyloseq objects##########################################

###########Bacteria##############

##Load raw data

otu_mat <- read.delim('bacteria_otu_result.tsv', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)
#generate taxonomy table from otu-table
tax_mat <- cbind(rownames(otu_mat),otu_mat$taxonomy)
rownames(tax_mat) <- rownames(otu_mat)

#Remove "k__", "p__" etc from taxonomy column

tax_mat[,2] <- str_remove_all(tax_mat[,2], str_c(c("k__", "p__", "c__", "o__", "f__", "g__", "s__"), collapse="|"))

#Split taxonomy column by ";"

split <- colsplit(tax_mat[,2], ";", c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))

#Remove olt taxonomy column
tax_mat <- tax_mat[,-c(2)]

tax_mat <- cbind(tax_mat,split)

tax_mat <- subset( tax_mat, select = -c(tax_mat) )

tax_mat <- as.matrix.data.frame(tax_mat)

#Remove taxonomy from otu table
otu_mat <- within(otu_mat,rm("taxonomy"))
otu_mat <- as.matrix(otu_mat)
map_mat <- read.delim('16S_2801_updated_mapping.txt', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)

##Construct individual tables is phyloseq format

OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)
TAX <- phyloseq::tax_table(tax_mat)
samples <- sample_data(map_mat)

##Combine in phyloseq object

PSB <- phyloseq(OTU, TAX, samples)

#### Perform CSS normalisation

PSB.CSS <- metagMisc::phyloseq_transform_css(PSB)

########### Virome ###########

##Load raw data 

#OTU table
otu_mat <- read.delim('virome_otu_result.tsv', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)

#04/11-2021 taxonomy from Josue
#tax_mat <- read.delim('taxonomy_lineage_tax.txt', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)
#2023 taxonomy from Frej Larsen 
# tax_mat <- read.delim('taxonomy_merged.tsv', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE,
#                       row.names = 1,
#                       header = F)
tax_mat <- read.delim('h_qual_clean_contigs_virus_summary.tsv', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE,
                      row.names = 1,
                      header = T) %>%
  dplyr::select(taxonomy) 

colnames(tax_mat) <- NULL

#04/11-2021 host predictions form Josue
host_mat <- read.delim('taxonomy_lineage_host.txt', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE,row.names = 1)


##Processing for phyloseq

#Merge old and new taxonomy columns. Updated taxonomy file contains only hits so otherwise all ucharacterised ASV will be lost when assembling phyloseq object.
#tax_mat <- dplyr::select(otu_mat,"taxonomy")

#tax_marged <- left_join(
#  rownames_to_column(tax_old,"ID"),
#  rownames_to_column(tax_mat,"ID"),
#  by="ID")

#Split taxonomy column by ";" - NB!: If any tax assignments have species add species column below

split <- colsplit(tax_mat[,1], ";", c("Domain","Realm","Kingdom","Phylum","Class","Order","Family","Genus","Species")) #%>% dplyr::select(-Species2) #Remove redundant row
#Remove lefter over "," - only if no species level assignments
split[] <- lapply(split, gsub, pattern=';', replacement='')

#Replace "" with "Uncalssified
split <- split %>% 
  mutate_all(~if_else(is.na(.), "", as.character(.))) %>%
  mutate_if(is.character, ~if_else(. == "", "Unclassified", .))

#Add split columns to tax_mat dataframe

tax_mat <- cbind(tax_mat,split)

#Remove old taxonomy and extra otu-id column

tax_mat <- tax_mat[,-c(1)]

#Convert to dataframe

tax_mat <- as.matrix.data.frame(tax_mat)

#Remove taxonomy from otu table
otu_mat <- within(otu_mat,rm("taxonomy"))
#Round all values
rownames <- rownames(otu_mat)
otu_mat <- mutate_all(otu_mat, round) %>% mutate_all(as.numeric)
rownames(otu_mat) <- rownames

map_mat <- read.delim('virome_cleaned_2801_mapping.txt', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)


##Construct individual tables is phyloseq format

OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)
TAX <- phyloseq::tax_table(tax_mat)
samples <- sample_data(map_mat)

##Combine in phyloseq object

PSV <- phyloseq(OTU, TAX, samples)

### CSS normalisation

PSV.CSS <- metagMisc::phyloseq_transform_css(PSV)

#Make metadata files for easy overview

meta.bac <- as.data.frame(as.matrix(sample_data(PSB)))

meta.vir <- as.data.frame(as.matrix(sample_data(PSV)))

#############Make host prediction file

split <- colsplit(host_mat[,1], ";", c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))

#Replace NA's with "Unknown" for species column
#split <- tidyr::replace_na(split, list(Species = "Unknown"))
split$Species <- "Unknown"

#Remove lefter over "," - only if no species level assignments
split[] <- lapply(split, gsub, pattern=';', replacement='')

#Add split columns to tax_mat dataframe

host_mat <- cbind(host_mat,split)

#Remove old taxonomy and extra otu-id column

host_mat <- host_mat[,-c(1)]

#Convert to dataframe

host_mat <- as.matrix.data.frame(host_mat)

##Combine in phyloseq object

HOST <- phyloseq::tax_table(host_mat)

PSV.host <- phyloseq(OTU, HOST, samples)

### CSS normalisation

PSV.CSS.host <- metagMisc::phyloseq_transform_css(PSV.host)

#########################PSV table without Realm - for Ampvis2

TAX.no.realm <- subset(tax_mat, select = -c(Realm)) %>% phyloseq::tax_table()

##Combine in phyloseq object

PSV.no.realm <- phyloseq(OTU, TAX.no.realm, samples)


#############Remove left over objects

rm(tax_mat, map_mat, otu_mat, samples, split, OTU, TAX)



