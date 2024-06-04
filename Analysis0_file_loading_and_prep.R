###################Listeria phage data analysis###################

#install.packages("tidyverse")
library(tidyverse)
library(ggpubr)
library(ggsci)
library(rstatix)
library(ampvis2)
library(plyr)
library(cowplot)
library(RVAideMemoire)
library(data.table)
library(DESeq2)
library(viridis)
library(vegan)
library(directlabels)
library(usdm)
#library(caret)
library(mixOmics)
#library(Rarefy)
#library(EnhancedVolcano)

#Load phyloseq to amvis2 converter
devtools::source_gist("8d0ca4206a66be7ff6d76fc4ab8e66c6")

#Coulour pallettes
col_fil <- pal_jco("default")(10)

col_scale <- scale_color_jco()

#Set working directory to script directory

#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##############################Load sequencing data#######################################

##Import files

#Bacteria

source("Phyloseq_import_and_CSS_bac_vir.R")

############Set data categories - and order

##Days as factor

PSB@sam_data$Days <- factor(PSB@sam_data$Days)
PSB.CSS@sam_data$Days <- factor(PSB.CSS@sam_data$Days)
PSV@sam_data$Days <- factor(PSV@sam_data$Days)
PSV.no.realm@sam_data$Days <- factor(PSV.no.realm@sam_data$Days)
PSV.CSS@sam_data$Days <- as.factor(PSV.CSS@sam_data$Days)

############Functions

# Save/load function for 'heavy' objects
sl <- function(name, ..., overwrite = FALSE, dir_path = here::here("results", "RData", subdir_name)) {
  # Possibility to add name as name or literal character string
  name <- as.character(substitute(name))
  assign(name, 
         if(file.exists(glue::glue("{dir_path}/{name}.Rds")) && !overwrite) {
           readRDS(glue::glue("{dir_path}/{name}.Rds"))
         }
         else { 
           dir.create(dir_path, showWarnings = F, recursive = T)
           saveRDS(..., file=glue::glue("{dir_path}/{name}.Rds"))
           readRDS(glue::glue("{dir_path}/{name}.Rds"))
         }, envir=.GlobalEnv)
}

#######################################Remove low read samples###############################################################

PSB.CSS = prune_samples(sample_sums(PSB) >= 100, PSB.CSS)
PSB = prune_samples(sample_sums(PSB) >= 100, PSB)

PSV.CSS = prune_samples(sample_sums(PSV)>=500, PSV.CSS)
PSV = prune_samples(sample_sums(PSV)>=500, PSV) # Not used since all samples had very high reads


