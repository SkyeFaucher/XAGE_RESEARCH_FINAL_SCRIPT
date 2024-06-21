######################################################################################################################
#                                   NEW XAGE ANALYSIS WITH UPDATED DATA
######################################################################################################################

# set working directory
setwd("/Users/Skye/Desktop/New_Xage/Second_Analysis/")

# loading the required packages
library(dplyr)
library(readxl)
library(writexl)

######################################################################################################################
#                               READING DIFFERENTIALLY EXPRESSED GENES INTO R 
######################################################################################################################

# xage1a differentially expressed genes
# the gene ID column was not labelled, so label was created
# remove identical xages from their lists of differentially expressed genes
xage1a_72_upregulated <- read.csv("xage1a_72_upregulated.csv", header = T)
colnames(xage1a_72_upregulated)[1] <- "Ensembl_gene_ID"

xage1a_120_upregulated <- read.csv("xage1a_120_upregulated.csv", header = T)
colnames(xage1a_120_upregulated)[1] <- "Ensembl_gene_ID"

xage1a_72_downregulated <- read.csv("xage1a_72_downregulated.csv", header = T)
colnames(xage1a_72_downregulated)[1] <- "Ensembl_gene_ID"
# removing xage1a from the list of differentially expressed genes
xage1a_72_downregulated <- xage1a_72_downregulated %>% dplyr::filter(xage1a_72_downregulated$Ensembl_gene_ID != "ENSG00000204379")

xage1a_120_downregulated <- read.csv("xage1a_120_downregulated.csv", header = T)
colnames(xage1a_120_downregulated)[1] <- "Ensembl_gene_ID"
xage1a_120_downregulated <- xage1a_120_downregulated %>% dplyr::filter(xage1a_120_downregulated$Ensembl_gene_ID != "ENSG00000204379")

# uploading the tables to desktop files in excel format
write_xlsx(xage1a_72_upregulated, "/Users/Skye/Desktop/New_Xage/Second_Analysis/Results/xage1a_72_upregulated.xlsx")
write_xlsx(xage1a_120_upregulated, "/Users/Skye/Desktop/New_Xage/Second_Analysis/Results/xage1a_120_upregulated.xlsx")
write_xlsx(xage1a_72_downregulated, "/Users/Skye/Desktop/New_Xage/Second_Analysis/Results/xage1a_72_downregulated.xlsx")
write_xlsx(xage1a_120_downregulated, "/Users/Skye/Desktop/New_Xage/Second_Analysis/Results/xage1a_120_downregulated.xlsx")


# xage2 differentially expressed genes
xage2_72_upregulated <- read.csv("xage2_72_upregulated.csv", header = T)
colnames(xage2_72_upregulated)[1] <- "Ensembl_gene_ID"

xage2_120_upregulated <- read.csv("xage2_120_upregulated.csv", header = T)
colnames(xage2_120_upregulated)[1] <- "Ensembl_gene_ID"

xage2_72_downregulated <- read.csv("xage2_72_downregulated.csv", header = T)
colnames(xage2_72_downregulated)[1] <- "Ensembl_gene_ID"
xage2_72_downregulated <- xage2_72_downregulated %>% dplyr::filter(xage2_72_downregulated$Ensembl_gene_ID != "ENSG00000155622")

xage2_120_downregulated <- read.csv("xage2_120_downregulated.csv", header = T)
colnames(xage2_120_downregulated)[1] <- "Ensembl_gene_ID"
xage2_120_downregulated <- xage2_120_downregulated %>% dplyr::filter(xage2_120_downregulated$Ensembl_gene_ID != "ENSG00000155622")

# uploading the tables to desktop files in excel format
write_xlsx(xage2_72_upregulated, "/Users/Skye/Desktop/New_Xage/Second_Analysis/Results/xage2_72_upregulated.xlsx")
write_xlsx(xage2_120_upregulated, "/Users/Skye/Desktop/New_Xage/Second_Analysis/Results/xage2_120_upregulated.xlsx")
write_xlsx(xage2_72_downregulated, "/Users/Skye/Desktop/New_Xage/Second_Analysis/Results/xage2_72_downregulated.xlsx")
write_xlsx(xage2_120_downregulated, "/Users/Skye/Desktop/New_Xage/Second_Analysis/Results/xage2_120_downregulated.xlsx")


######################################################################################################################
#                     OBTAINING FULL LIST OF DIFFERENTIALLY EXPRESSED GENES                                               
######################################################################################################################

# installing and loading the required packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install('grimbough/biomaRt')

# loading the required packages
library(biomaRt)
library(annotables)
library(tidyverse)
library(EnsDb.Hsapiens.v86)
library(org.Hs.eg.db)

# using the gene ids to obtain the gene names and gene symbols
listEnsembl()
ensembl <- useEnsembl(biomart = "genes")
datasets <- listDatasets(ensembl)
ensembl.connected <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# observing what information about the genes we can obtain from the database
attributes <- listAttributes(ensembl.connected)

# identifying what format our input data is in
filters <- listFilters(ensembl.connected)

# xage1a differentially expressed genes
# compiling new data tables with gene ids, gene symbols and gene description
xage1a_72_upregulated_gene_list <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", 
                                                                   "uniprot_gn_symbol", "description"),  
                                                filters = "ensembl_gene_id", 
                                                values = xage1a_72_upregulated$Ensembl_gene_ID, 
                                                mart = ensembl.connected))

xage1a_120_upregulated_gene_list <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", 
                                                                   "uniprot_gn_symbol", "description"),  
                                                    filters = "ensembl_gene_id", 
                                                    values = xage1a_120_upregulated$Ensembl_gene_ID, 
                                                    mart = ensembl.connected))

xage1a_72_downregulated_gene_list <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", 
                                                                   "uniprot_gn_symbol", "description"),  
                                                    filters = "ensembl_gene_id", 
                                                    values = xage1a_72_downregulated$Ensembl_gene_ID, 
                                                    mart = ensembl.connected))

xage1a_120_downregulated_gene_list <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", 
                                                                   "uniprot_gn_symbol", "description"),  
                                                    filters = "ensembl_gene_id", 
                                                    values = xage1a_120_downregulated$Ensembl_gene_ID, 
                                                    mart = ensembl.connected))

# uploading the new data tables to research folder on desktop
write_xlsx(xage1a_72_upregulated_gene_list, "xage1a_72_upregulated_gene_list.xlsx")
write_xlsx(xage1a_120_upregulated_gene_list, "xage1a_120_upregulated_gene_list.xlsx")
write_xlsx(xage1a_72_downregulated_gene_list, "xage1a_72_downregulated_gene_list.xlsx")
write_xlsx(xage1a_120_downregulated_gene_list, "xage1a_120_downregulated_gene_list.xlsx")


# xage2 differentially expressed genes
# compiling new data tables with gene ids, gene symbols and gene description
xage2_72_upregulated_gene_list <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", 
                                                                   "uniprot_gn_symbol", "description"),  
                                                    filters = "ensembl_gene_id", 
                                                    values = xage2_72_upregulated$Ensembl_gene_ID, 
                                                    mart = ensembl.connected))

xage2_120_upregulated_gene_list <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", 
                                                                    "uniprot_gn_symbol", "description"),  
                                                     filters = "ensembl_gene_id", 
                                                     values = xage2_120_upregulated$Ensembl_gene_ID, 
                                                     mart = ensembl.connected))

xage2_72_downregulated_gene_list <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", 
                                                                     "uniprot_gn_symbol", "description"),  
                                                      filters = "ensembl_gene_id", 
                                                      values = xage2_72_downregulated$Ensembl_gene_ID, 
                                                      mart = ensembl.connected))

xage2_120_downregulated_gene_list <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", 
                                                                      "uniprot_gn_symbol", "description"),  
                                                       filters = "ensembl_gene_id", 
                                                       values = xage2_120_downregulated$Ensembl_gene_ID, 
                                                       mart = ensembl.connected))

# uploading the new data tables to research folder on desktop
write_xlsx(xage2_72_upregulated_gene_list, "xage2_72_upregulated_gene_list.xlsx")
write_xlsx(xage2_120_upregulated_gene_list, "xage2_120_upregulated_gene_list.xlsx")
write_xlsx(xage2_72_downregulated_gene_list, "xage2_72_downregulated_gene_list.xlsx")
write_xlsx(xage2_120_downregulated_gene_list, "xage2_120_downregulated_gene_list.xlsx")



######################################################################################################################
#                       OBTAINING TOP 10 LIST OF DIFFERENTIALLY EXPRESSED GENES 
######################################################################################################################

# sorting based on p-value for significance
# xage1a differentially expressed genes
xage1a_72_upregulated_top_10 <- xage1a_72_upregulated[order(xage1a_72_upregulated$PValue), ][1:10, 1:8]

xage1a_120_upregulated_top_10 <- xage1a_120_upregulated[order(xage1a_120_upregulated$PValue), ][1:10, 1:8]

xage1a_72_downregulated_top_10 <- xage1a_72_downregulated[order(xage1a_72_downregulated$PValue), ][1:10, 1:8]

xage1a_120_downregulated_top_10 <- xage1a_120_downregulated[order(xage1a_120_downregulated$PValue), ][1:10, 1:8]

# compiling new data tables with gene ids, gene symbols and gene description
xage1a_72_upregulated_top_10_gene_list <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", 
                                                                          "uniprot_gn_symbol", "description"),  
                                                           filters = "ensembl_gene_id", 
                                                           values = xage1a_72_upregulated_top_10$Ensembl_gene_ID, 
                                                           mart = ensembl.connected))

xage1a_120_upregulated_top_10_gene_list <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", 
                                                                    "uniprot_gn_symbol", "description"),  
                                                     filters = "ensembl_gene_id", 
                                                     values = xage1a_120_upregulated_top_10$Ensembl_gene_ID, 
                                                     mart = ensembl.connected))

xage1a_72_downregulated_top_10_gene_list <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", 
                                                                     "uniprot_gn_symbol", "description"),  
                                                      filters = "ensembl_gene_id", 
                                                      values = xage1a_72_downregulated_top_10$Ensembl_gene_ID, 
                                                      mart = ensembl.connected))

xage1a_120_downregulated_top_10_gene_list <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", 
                                                                      "uniprot_gn_symbol", "description"),  
                                                       filters = "ensembl_gene_id", 
                                                       values = xage1a_120_downregulated_top_10$Ensembl_gene_ID, 
                                                       mart = ensembl.connected))

# uploading the new data tables to research folder on desktop
write_xlsx(xage1a_72_upregulated_top_10_gene_list, "xage1a_72_upregulated_top_10_gene_list.xlsx")
write_xlsx(xage1a_120_upregulated_top_10_gene_list, "xage1a_120_upregulated_top_10_gene_list.xlsx")
write_xlsx(xage1a_72_downregulated_top_10_gene_list, "xage1a_72_downregulated_top_10_gene_list.xlsx")
write_xlsx(xage1a_120_downregulated_top_10_gene_list, "xage1a_120_downregulated_top_10_gene_list.xlsx")


# xage2 differentially expressed genes
xage2_72_upregulated_top_10 <- xage2_72_upregulated[order(xage2_72_upregulated$PValue), ][1:10, 1:8]

xage2_120_upregulated_top_10 <- xage2_120_upregulated[order(xage2_120_upregulated$PValue), ][1:10, 1:8]

xage2_72_downregulated_top_10 <- xage2_72_downregulated[order(xage2_72_downregulated$PValue), ][1:10, 1:8]

xage2_120_downregulated_top_10 <- xage2_120_downregulated[order(xage2_120_downregulated$PValue), ][1:10, 1:8]

# compiling new data tables with gene ids, gene symbols and gene description
xage2_72_upregulated_top_10_gene_list <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", 
                                                                          "uniprot_gn_symbol", "description"),  
                                                           filters = "ensembl_gene_id", 
                                                           values = xage2_72_upregulated_top_10$Ensembl_gene_ID, 
                                                           mart = ensembl.connected))

xage2_120_upregulated_top_10_gene_list <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", 
                                                                           "uniprot_gn_symbol", "description"),  
                                                            filters = "ensembl_gene_id", 
                                                            values = xage2_120_upregulated_top_10$Ensembl_gene_ID, 
                                                            mart = ensembl.connected))

xage2_72_downregulated_top_10_gene_list <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", 
                                                                            "uniprot_gn_symbol", "description"),  
                                                             filters = "ensembl_gene_id", 
                                                             values = xage2_72_downregulated_top_10$Ensembl_gene_ID, 
                                                             mart = ensembl.connected))

xage2_120_downregulated_top_10_gene_list <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", 
                                                                             "uniprot_gn_symbol", "description"),  
                                                              filters = "ensembl_gene_id", 
                                                              values = xage2_120_downregulated_top_10$Ensembl_gene_ID, 
                                                              mart = ensembl.connected))

# uploading the new data tables to research folder on desktop
write_xlsx(xage2_72_upregulated_top_10_gene_list, "xage2_72_upregulated_top_10_gene_list.xlsx")
write_xlsx(xage2_120_upregulated_top_10_gene_list, "xage2_120_upregulated_top_10_gene_list.xlsx")
write_xlsx(xage2_72_downregulated_top_10_gene_list, "xage2_72_downregulated_top_10_gene_list.xlsx")
write_xlsx(xage2_120_downregulated_top_10_gene_list, "xage2_120_downregulated_top_10_gene_list.xlsx")



######################################################################################################################
#                                         INTER GENE INTERSECTIONS
######################################################################################################################

######################################################################################################################
#           INTERSECTING COMMON DIFFERENTIALLY EXPRESSED GENES BETWEEN BOTH XAGES
######################################################################################################################

# isolating genes differentially expressed in both xage1a and xage2 after 72 and 120 hrs
xage_72_upregulated_intersection <- data.frame(intersect(xage1a_72_upregulated_gene_list$ensembl_gene_id,
                                                         xage2_72_upregulated_gene_list$ensembl_gene_id))
colnames(xage_72_upregulated_intersection)[1] <- "Ensembl_gene_ID"


xage_120_upregulated_intersection <- data.frame(intersect(xage1a_120_upregulated_gene_list$ensembl_gene_id,
                                                         xage2_120_upregulated_gene_list$ensembl_gene_id))
colnames(xage_120_upregulated_intersection)[1] <- "Ensembl_gene_ID"


xage_72_downregulated_intersection <- data.frame(intersect(xage1a_72_downregulated_gene_list$ensembl_gene_id,
                                                          xage2_72_downregulated_gene_list$ensembl_gene_id))
colnames(xage_72_downregulated_intersection)[1] <- "Ensembl_gene_ID"


xage_120_downregulated_intersection <- data.frame(intersect(xage1a_120_downregulated_gene_list$ensembl_gene_id,
                                                          xage2_120_downregulated_gene_list$ensembl_gene_id))
colnames(xage_120_downregulated_intersection)[1] <- "Ensembl_gene_ID"


# obtaining the gene names and descriptions for the intersecting gene lists
xage_72_upregulated_intersecting_genes <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", 
                                                                          "uniprot_gn_symbol", "description"),  
                                                           filters = "ensembl_gene_id", 
                                                           values = xage_72_upregulated_intersection$Ensembl_gene_ID, 
                                                           mart = ensembl.connected))

xage_120_upregulated_intersecting_genes <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", 
                                                                           "uniprot_gn_symbol", "description"),  
                                                            filters = "ensembl_gene_id", 
                                                            values = xage_120_upregulated_intersection$Ensembl_gene_ID, 
                                                            mart = ensembl.connected))

xage_72_downregulated_intersecting_genes <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", 
                                                                            "uniprot_gn_symbol", "description"),  
                                                             filters = "ensembl_gene_id", 
                                                             values = xage_72_downregulated_intersection$Ensembl_gene_ID, 
                                                             mart = ensembl.connected))

xage_120_downregulated_intersecting_genes <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", 
                                                                             "uniprot_gn_symbol", "description"),  
                                                              filters = "ensembl_gene_id", 
                                                              values = xage_120_downregulated_intersection$Ensembl_gene_ID, 
                                                              mart = ensembl.connected))


# uploading the new data tables to research folder on desktop
write_xlsx(xage_72_upregulated_intersecting_genes, "xage_72_upregulated_intersecting_genes.xlsx")

write_xlsx(xage_120_upregulated_intersecting_genes, "xage_120_upregulated_intersecting_genes.xlsx")

write_xlsx(xage_72_downregulated_intersecting_genes, "xage_72_downregulated_intersecting_genes.xlsx")

write_xlsx(xage_120_downregulated_intersecting_genes, "xage_120_downregulated_intersecting_genes.xlsx")



######################################################################################################################
#      INTERSECTING COMMON DIFFERENTIALLY EXPRESSED GENES BETWEEN XAGES WITH CONTRASTING REGULATORY EFFECTS
######################################################################################################################

# isolating differentially expressed genes with different regulatory effects on xage1a and xage2
xage1a_up_xage2_down_72_intersection <- data.frame(intersect(xage1a_72_upregulated_gene_list$ensembl_gene_id, 
                                                             xage2_72_downregulated_gene_list$ensembl_gene_id))
colnames(xage1a_up_xage2_down_72_intersection)[1] <- "Ensembl_gene_ID"


xage1a_up_xage2_down_120_intersection <- data.frame(intersect(xage1a_120_upregulated_gene_list$ensembl_gene_id, 
                                                              xage2_120_downregulated_gene_list$ensembl_gene_id))
colnames(xage1a_up_xage2_down_120_intersection)[1] <- "Ensembl_gene_ID"


xage1a_down_xage2_up_72_intersection <- data.frame(intersect(xage1a_72_downregulated_gene_list$ensembl_gene_id, 
                                                             xage2_72_upregulated_gene_list$ensembl_gene_id))
colnames(xage1a_down_xage2_up_72_intersection)[1] <- "Ensembl_gene_ID"


xage1a_down_xage2_up_120_intersection <- data.frame(intersect(xage1a_120_downregulated_gene_list$ensembl_gene_id, 
                                                              xage2_120_upregulated_gene_list$ensembl_gene_id))
colnames(xage1a_down_xage2_up_120_intersection)[1] <- "Ensembl_gene_ID"


# obtaining the gene names and descriptions for the intersecting gene lists
xage1a_up_xage2_down_120_intersecting_genes <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", 
                                                                               "uniprot_gn_symbol", "description"),  
                                                                filters = "ensembl_gene_id", 
                                                                values = xage1a_up_xage2_down_120_intersection$Ensembl_gene_ID, 
                                                                mart = ensembl.connected))

# uploading the new data tables to research folder on desktop
write_xlsx(xage1a_up_xage2_down_120_intersecting_genes, "xage1a_up_xage2_down_120_intersecting_genes.xlsx")



######################################################################################################################
# INTERSECTING COMMON DIFFERENTIALLY EXPRESSED GENES BETWEEN XAGES WITH CONTRASTING REGULATORY EFFECTS & TIME POINTS
######################################################################################################################

# isolating differentially expressed genes with contrasting regulatory effects on xage1a and xage2 at differing time points
xage1a_72_up_xage2_120_down_intersection <- data.frame(intersect(xage1a_72_upregulated_gene_list$ensembl_gene_id, 
                                                                 xage2_120_downregulated_gene_list$ensembl_gene_id))
colnames(xage1a_72_up_xage2_120_down_intersection)[1] <- "Ensembl_gene_ID"


xage1a_120_up_xage2_72_down_intersection <- data.frame(intersect(xage1a_120_upregulated_gene_list$ensembl_gene_id, 
                                                                 xage2_72_downregulated_gene_list$ensembl_gene_id))
colnames(xage1a_120_up_xage2_72_down_intersection)[1] <- "Ensembl_gene_ID"


xage1a_72_down_xage2_120_up_intersection <- data.frame(intersect(xage1a_72_downregulated_gene_list$ensembl_gene_id, 
                                                                 xage2_120_upregulated_gene_list$ensembl_gene_id))
colnames(xage1a_72_down_xage2_120_up_intersection)[1] <- "Ensembl_gene_ID"


xage1a_120_down_xage2_72_up_intersection <- data.frame(intersect(xage1a_120_downregulated_gene_list$ensembl_gene_id, 
                                                                 xage2_72_upregulated_gene_list$ensembl_gene_id))
colnames(xage1a_120_down_xage2_72_up_intersection)[1] <- "Ensembl_gene_ID"


# obtaining the gene names and descriptions for the intersecting gene lists
xage1a_72_up_xage2_120_down_intersecting_genes <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", 
                                                                                  "uniprot_gn_symbol", "description"),  
                                                                   filters = "ensembl_gene_id", 
                                                                   values = xage1a_72_up_xage2_120_down_intersection$Ensembl_gene_ID, 
                                                                   mart = ensembl.connected))


xage1a_120_up_xage2_72_down_intersecting_genes <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", 
                                                                                  "uniprot_gn_symbol", "description"),  
                                                                   filters = "ensembl_gene_id", 
                                                                   values = xage1a_120_up_xage2_72_down_intersection$Ensembl_gene_ID, 
                                                                   mart = ensembl.connected))


xage1a_72_down_xage2_120_up_intersecting_genes <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", 
                                                                                  "uniprot_gn_symbol", "description"),  
                                                                   filters = "ensembl_gene_id", 
                                                                   values = xage1a_72_down_xage2_120_up_intersection$Ensembl_gene_ID, 
                                                                   mart = ensembl.connected))

xage1a_120_down_xage2_72_up_intersecting_genes <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", 
                                                                                  "uniprot_gn_symbol", "description"),  
                                                                   filters = "ensembl_gene_id", 
                                                                   values = xage1a_120_down_xage2_72_up_intersection$Ensembl_gene_ID, 
                                                                   mart = ensembl.connected))


# uploading the new data tables to research folder on desktop
write_xlsx(xage1a_72_up_xage2_120_down_intersecting_genes, "xage1a_72_up_xage2_120_down_intersecting_genes.xlsx")
write_xlsx(xage1a_120_up_xage2_72_down_intersecting_genes, "xage1a_120_up_xage2_72_down_intersecting_genes.xlsx")
write_xlsx(xage1a_72_down_xage2_120_up_intersecting_genes, "xage1a_72_down_xage2_120_up_intersecting_genes.xlsx")
write_xlsx(xage1a_120_down_xage2_72_up_intersecting_genes, "xage1a_120_down_xage2_72_up_intersecting_genes.xlsx")



######################################################################################################################
# INTERSECTING COMMON DIFFERENTIALLY EXPRESSED GENES BETWEEN XAGES WITH CONTRASTING TIME POINTS
######################################################################################################################

# isolating differentially expressed genes with different regulatory effects on xage1a and xage2
xage1a_72_xage2_120_up_intersection <- data.frame(intersect(xage1a_72_upregulated_gene_list$ensembl_gene_id, 
                                                            xage2_120_upregulated_gene_list$ensembl_gene_id))
colnames(xage1a_72_xage2_120_up_intersection)[1] <- "Ensembl_gene_ID"


xage1a_72_xage2_120_down_intersection <- data.frame(intersect(xage1a_72_downregulated_gene_list$ensembl_gene_id, 
                                                              xage2_120_downregulated_gene_list$ensembl_gene_id))
colnames(xage1a_72_xage2_120_down_intersection)[1] <- "Ensembl_gene_ID"


xage1a_120_xage2_72_up_intersection <- data.frame(intersect(xage1a_120_upregulated_gene_list$ensembl_gene_id, 
                                                            xage2_72_upregulated_gene_list$ensembl_gene_id))
colnames(xage1a_120_xage2_72_up_intersection)[1] <- "Ensembl_gene_ID"


xage1a_120_xage2_72_down_intersection <- data.frame(intersect(xage1a_120_downregulated_gene_list$ensembl_gene_id, 
                                                              xage2_72_downregulated_gene_list$ensembl_gene_id))
colnames(xage1a_120_xage2_72_down_intersection)[1] <- "Ensembl_gene_ID"


# obtaining the gene names and descriptions for the intersecting gene lists
xage1a_72_xage2_120_up_intersecting_genes <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", 
                                                                             "uniprot_gn_symbol", "description"),  
                                                              filters = "ensembl_gene_id", 
                                                              values = xage1a_72_xage2_120_up_intersection$Ensembl_gene_ID, 
                                                              mart = ensembl.connected))

xage1a_72_xage2_120_down_intersecting_genes <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", 
                                                                               "uniprot_gn_symbol", "description"),  
                                                                filters = "ensembl_gene_id", 
                                                                values = xage1a_72_xage2_120_down_intersection$Ensembl_gene_ID, 
                                                                mart = ensembl.connected))

xage1a_120_xage2_72_up_intersecting_genes <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", 
                                                                             "uniprot_gn_symbol", "description"),  
                                                              filters = "ensembl_gene_id", 
                                                              values = xage1a_120_xage2_72_up_intersection$Ensembl_gene_ID, 
                                                              mart = ensembl.connected))

xage1a_120_xage2_72_down_intersecting_genes <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", 
                                                                               "uniprot_gn_symbol", "description"),  
                                                                filters = "ensembl_gene_id", 
                                                                values = xage1a_120_xage2_72_down_intersection$Ensembl_gene_ID, 
                                                                mart = ensembl.connected))

# uploading the new data tables to research folder on desktop
write_xlsx(xage1a_72_xage2_120_up_intersecting_genes, "xage1a_72_xage2_120_up_intersecting_genes.xlsx")
write_xlsx(xage1a_72_xage2_120_down_intersecting_genes, "xage1a_72_xage2_120_down_intersecting_genes.xlsx")
write_xlsx(xage1a_120_xage2_72_up_intersecting_genes, "xage1a_120_xage2_72_up_intersecting_genes.xlsx")
write_xlsx(xage1a_120_xage2_72_down_intersecting_genes, "xage1a_120_xage2_72_down_intersecting_genes.xlsx")


######################################################################################################################
#                                         INTRA GENE INTERSECTIONS
######################################################################################################################

######################################################################################################################
#       INTERSECTING COMMON DIFFERENTIALLY EXPRESSED GENES FOR EACH XAGE WITH SAME EFFECT AT DIFFERENT TIME POINTS
######################################################################################################################

# isolating differentially expressed genes at both time points in xage1a and xage2 respectively
xage1a_upregulated_intersection <- data.frame(intersect(xage1a_72_upregulated_gene_list$ensembl_gene_id, 
                                                        xage1a_120_upregulated_gene_list$ensembl_gene_id))
colnames(xage1a_upregulated_intersection)[1] <- "Ensembl_gene_ID"


xage1a_downregulated_intersection <- data.frame(intersect(xage1a_72_downregulated_gene_list$ensembl_gene_id, 
                                                          xage1a_120_downregulated_gene_list$ensembl_gene_id))
colnames(xage1a_downregulated_intersection)[1] <- "Ensembl_gene_ID"


xage2_upregulated_intersection <- data.frame(intersect(xage2_72_upregulated_gene_list$ensembl_gene_id, 
                                                       xage2_120_upregulated_gene_list$ensembl_gene_id))
colnames(xage2_upregulated_intersection)[1] <- "Ensembl_gene_ID"


xage2_downregulated_intersection <- data.frame(intersect(xage2_72_downregulated_gene_list$ensembl_gene_id, 
                                                         xage2_120_downregulated_gene_list$ensembl_gene_id))
colnames(xage2_downregulated_intersection)[1] <- "Ensembl_gene_ID"


# obtaining the gene names and descriptions for the intersecting gene lists
xage1a_upregulated_intersecting_genes <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", 
                                                                         "uniprot_gn_symbol", "description"),  
                                                          filters = "ensembl_gene_id", 
                                                          values = xage1a_upregulated_intersection$Ensembl_gene_ID, 
                                                          mart = ensembl.connected))

xage1a_downregulated_intersecting_genes <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", 
                                                                           "uniprot_gn_symbol", "description"),  
                                                            filters = "ensembl_gene_id", 
                                                            values = xage1a_downregulated_intersection$Ensembl_gene_ID, 
                                                            mart = ensembl.connected))

xage2_upregulated_intersecting_genes <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", 
                                                                        "uniprot_gn_symbol", "description"),  
                                                         filters = "ensembl_gene_id", 
                                                         values = xage2_upregulated_intersection$Ensembl_gene_ID, 
                                                         mart = ensembl.connected))

xage2_downregulated_intersecting_genes <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", 
                                                                          "uniprot_gn_symbol", "description"),  
                                                           filters = "ensembl_gene_id", 
                                                           values = xage2_downregulated_intersection$Ensembl_gene_ID, 
                                                           mart = ensembl.connected))


# uploading the new data tables to research folder on desktop
write_xlsx(xage1a_upregulated_intersecting_genes, "xage1a_upregulated_intersecting_genes.xlsx")

write_xlsx(xage1a_downregulated_intersecting_genes, "xage1a_downregulated_intersecting_genes.xlsx")

write_xlsx(xage2_upregulated_intersecting_genes, "xage2_upregulated_intersecting_genes.xlsx")

write_xlsx(xage2_downregulated_intersecting_genes, "xage2_downregulated_intersecting_genes.xlsx")


######################################################################################################################
#       INTERSECTING COMMON DIFFERENTIALLY EXPRESSED GENES FOR EACH XAGE WITH DIFFERENT EFFECT AT SAME TIME POINTS
######################################################################################################################

# isolating differentially expressed genes at both time points in xage1a and xage2 respectively
xage1a_72_intersection <- data.frame(intersect(xage1a_72_upregulated_gene_list$ensembl_gene_id, 
                                                        xage1a_72_downregulated_gene_list$ensembl_gene_id))
colnames(xage1a_72_intersection)[1] <- "Ensembl_gene_ID"


xage1a_120_intersection <- data.frame(intersect(xage1a_120_upregulated_gene_list$ensembl_gene_id, 
                                                          xage1a_120_downregulated_gene_list$ensembl_gene_id))
colnames(xage1a_120_intersection)[1] <- "Ensembl_gene_ID"


xage2_72_intersection <- data.frame(intersect(xage2_72_upregulated_gene_list$ensembl_gene_id, 
                                                       xage2_72_downregulated_gene_list$ensembl_gene_id))
colnames(xage2_72_intersection)[1] <- "Ensembl_gene_ID"


xage2_120_intersection <- data.frame(intersect(xage2_120_upregulated_gene_list$ensembl_gene_id, 
                                                         xage2_120_downregulated_gene_list$ensembl_gene_id))
colnames(xage2_120_intersection)[1] <- "Ensembl_gene_ID"

# obtaining the gene names and descriptions for the intersecting gene lists
xage1a_72_intersecting_genes <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", 
                                                                         "uniprot_gn_symbol", "description"),  
                                                          filters = "ensembl_gene_id", 
                                                          values = xage1a_72_intersection$Ensembl_gene_ID, 
                                                          mart = ensembl.connected))

xage1a_120_intersecting_genes <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", 
                                                                           "uniprot_gn_symbol", "description"),  
                                                            filters = "ensembl_gene_id", 
                                                            values = xage1a_120_intersection$Ensembl_gene_ID, 
                                                            mart = ensembl.connected))

xage2_72_intersecting_genes <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", 
                                                                        "uniprot_gn_symbol", "description"),  
                                                         filters = "ensembl_gene_id", 
                                                         values = xage2_72_intersection$Ensembl_gene_ID, 
                                                         mart = ensembl.connected))

xage2_120_intersecting_genes <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", 
                                                                          "uniprot_gn_symbol", "description"),  
                                                           filters = "ensembl_gene_id", 
                                                           values = xage2_120_intersection$Ensembl_gene_ID, 
                                                           mart = ensembl.connected))


# uploading the new data tables to research folder on desktop
write_xlsx(xage1a_72_intersecting_genes, "xage1a_72_intersecting_genes.xlsx")

write_xlsx(xage1a_120_intersecting_genes, "xage1a_120_intersecting_genes.xlsx")

write_xlsx(xage2_72_intersecting_genes, "xage2_72_intersecting_genes.xlsx")

write_xlsx(xage2_120_intersecting_genes, "xage2_120_intersecting_genes.xlsx")


######################################################################################################################
#   INTERSECTING COMMON DIFFERENTIALLY EXPRESSED GENES FOR EACH XAGE WITH DIFFERENT EFFECT AT DIFFERENT TIME POINTS
######################################################################################################################

# isolating differentially expressed genes at both time points in xage1a and xage2 respectively
xage1a_72_up_and_120_down_intersection <- data.frame(intersect(xage1a_72_upregulated_gene_list$ensembl_gene_id, 
                                               xage1a_120_downregulated_gene_list$ensembl_gene_id))
colnames(xage1a_72_up_and_120_down_intersection)[1] <- "Ensembl_gene_ID"


xage1a_120_up_and_72_down_intersection <- data.frame(intersect(xage1a_120_upregulated_gene_list$ensembl_gene_id, 
                                                xage1a_72_downregulated_gene_list$ensembl_gene_id))
colnames(xage1a_120_up_and_72_down_intersection)[1] <- "Ensembl_gene_ID"


xage2_72_up_and_120_down_intersection <- data.frame(intersect(xage2_72_upregulated_gene_list$ensembl_gene_id, 
                                              xage2_120_downregulated_gene_list$ensembl_gene_id))
colnames(xage2_72_up_and_120_down_intersection)[1] <- "Ensembl_gene_ID"


xage2_120_up_and_72_down_intersection <- data.frame(intersect(xage2_120_upregulated_gene_list$ensembl_gene_id, 
                                               xage2_72_downregulated_gene_list$ensembl_gene_id))
colnames(xage2_120_up_and_72_down_intersection)[1] <- "Ensembl_gene_ID"


# obtaining the gene names and descriptions for the intersecting gene lists
xage1a_72_up_and_120_down_intersecting_genes <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", 
                                                                "uniprot_gn_symbol", "description"),  
                                                 filters = "ensembl_gene_id", 
                                                 values = xage1a_72_up_and_120_down_intersection$Ensembl_gene_ID, 
                                                 mart = ensembl.connected))

xage1a_120_up_and_72_down_intersecting_genes <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", 
                                                                 "uniprot_gn_symbol", "description"),  
                                                  filters = "ensembl_gene_id", 
                                                  values = xage1a_120_up_and_72_down_intersection$Ensembl_gene_ID, 
                                                  mart = ensembl.connected))

xage2_72_up_and_120_down_intersecting_genes <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", 
                                                               "uniprot_gn_symbol", "description"),  
                                                filters = "ensembl_gene_id", 
                                                values = xage2_72_up_and_120_down_intersection$Ensembl_gene_ID, 
                                                mart = ensembl.connected))

xage2_120_up_and_72_down_intersecting_genes <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", 
                                                                "uniprot_gn_symbol", "description"),  
                                                 filters = "ensembl_gene_id", 
                                                 values = xage2_120_up_and_72_down_intersection$Ensembl_gene_ID, 
                                                 mart = ensembl.connected))


# uploading the new data tables to research folder on desktop
write_xlsx(xage1a_72_up_and_120_down_intersecting_genes, "xage1a_72_up_and_120_down_intersecting_genes.xlsx")

write_xlsx(xage1a_120_up_and_72_down_intersecting_genes, "xage1a_120_up_and_72_down_intersecting_genes.xlsx")

write_xlsx(xage2_72_up_and_120_down_intersecting_genes, "xage2_72_up_and_120_down_intersecting_genes.xlsx")

write_xlsx(xage2_120_up_and_72_down_intersecting_genes, "xage2_120_up_and_72_down_intersecting_genes.xlsx")



######################################################################################################################
#                      INVESTIGATING FOXP4 AS A POTENTIAL TRANSCRIPTION FACTOR OF XAGE
######################################################################################################################

######################################################################################################################
#                         IDENTIFYING THE PROMOTER REGIONS OF THE HUMAN GENOME
######################################################################################################################

# obtaining the promoter region from the human genome
# Read data file into R
human_genome <- read.table("/Users/Skye/Desktop/New_Xage/human_genome", header = FALSE, sep = "\t")

# number of rows in the file
nrow(human_genome)  # returned 266064

# assigning file to a variable, obtaining the promoter region and saving to a bed file
g <- human_genome

for ( k in 1:266064){ if(g[k,3] == "+") 
{ TSS = g[k,4]}else{ TSS = g[k,5]}; 
  write( paste(c(g[k,2], TSS-1000, TSS+1000, g[k,1]), collapse = "\t"),
         file = "human_genome_promoter_region.bed", append = TRUE, sep = "\n") }


######################################################################################################################
#                IDENTIFYING THE FOXP4 BINDING SITES IN THE PROMOTER REGIONS OF THE GENOME
######################################################################################################################

# using terminal (linux), compare human genome promoter file to FOXP4 ChipSeq ranked peaks file
# word file "XAGE Analysis in Terminal" contains the steps and code used


######################################################################################################################
#                           CONVERTING ENSEMBL GENES TO ENSEMBL TRANSCRIPTS
######################################################################################################################

# reading required files into R
# xage1a differentially expressed genes
xage1a_72_upregulated_gene_list <- read_excel("xage1a_72_upregulated_gene_list.xlsx")
xage1a_120_upregulated_gene_list <- read_excel("xage1a_120_upregulated_gene_list.xlsx")
xage1a_72_downregulated_gene_list <- read_excel("xage1a_72_downregulated_gene_list.xlsx")
xage1a_120_downregulated_gene_list <- read_excel("xage1a_120_downregulated_gene_list.xlsx")

# xage2 differentially expressed genes
xage2_72_upregulated_gene_list <- read_excel("xage2_72_upregulated_gene_list.xlsx")
xage2_120_upregulated_gene_list <- read_excel("xage2_120_upregulated_gene_list.xlsx")
xage2_72_downregulated_gene_list <- read_excel("xage2_72_downregulated_gene_list.xlsx")
xage2_120_downregulated_gene_list <- read_excel("xage2_120_downregulated_gene_list.xlsx")


# converting ens genes to ens transcripts and obtaining other transcript information using biomart
listEnsembl()
ensembl <- useEnsembl(biomart = "genes")
datasets <- listDatasets(ensembl)
ensembl.connected <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

attributes <- listAttributes(ensembl.connected)
filters <- listFilters(ensembl.connected)


# xage1a
xage1a_72_upregulated_transcripts <- data.frame(getBM(attributes = c("ensembl_transcript_id_version", "chromosome_name", "transcript_start", "transcript_end",
                                                                     "strand"),  
                                                      filters = "ensembl_gene_id", 
                                                      values = xage1a_72_upregulated_gene_list$ensembl_gene_id, 
                                                      mart = ensembl.connected))

xage1a_72_downregulated_transcripts <- data.frame(getBM(attributes = c("ensembl_transcript_id_version", "chromosome_name", "transcript_start", "transcript_end",
                                                                       "strand"),  
                                                        filters = "ensembl_gene_id", 
                                                        values = xage1a_72_downregulated_gene_list$ensembl_gene_id, 
                                                        mart = ensembl.connected))

xage1a_120_upregulated_transcripts <- data.frame(getBM(attributes = c("ensembl_transcript_id_version", "chromosome_name", "transcript_start", "transcript_end",
                                                                      "strand"),  
                                                       filters = "ensembl_gene_id", 
                                                       values = xage1a_120_upregulated_gene_list$ensembl_gene_id, 
                                                       mart = ensembl.connected))

xage1a_120_downregulated_transcripts <- data.frame(getBM(attributes = c("ensembl_transcript_id_version", "chromosome_name", "transcript_start", "transcript_end",
                                                                        "strand"),  
                                                         filters = "ensembl_gene_id", 
                                                         values = xage1a_120_downregulated_gene_list$ensembl_gene_id, 
                                                         mart = ensembl.connected))

# xage2
xage2_72_upregulated_transcripts <- data.frame(getBM(attributes = c("ensembl_transcript_id_version", "chromosome_name", "transcript_start", "transcript_end", 
                                                                    "strand"),  
                                                     filters = "ensembl_gene_id", 
                                                     values = xage2_72_upregulated_gene_list$ensembl_gene_id, 
                                                     mart = ensembl.connected))

xage2_72_downregulated_transcripts <- data.frame(getBM(attributes = c("ensembl_transcript_id_version", "chromosome_name", "transcript_start", "transcript_end",
                                                                      "strand"),  
                                                       filters = "ensembl_gene_id", 
                                                       values = xage2_72_downregulated_gene_list$ensembl_gene_id, 
                                                       mart = ensembl.connected))

xage2_120_upregulated_transcripts <- data.frame(getBM(attributes = c("ensembl_transcript_id_version", "chromosome_name", "transcript_start", "transcript_end",
                                                                     "strand"),  
                                                      filters = "ensembl_gene_id", 
                                                      values = xage2_120_upregulated_gene_list$ensembl_gene_id, 
                                                      mart = ensembl.connected))

xage2_120_downregulated_transcripts <- data.frame(getBM(attributes = c("ensembl_transcript_id_version", "chromosome_name", "transcript_start", "transcript_end",
                                                                       "strand"),  
                                                        filters = "ensembl_gene_id", 
                                                        values = xage2_120_downregulated_gene_list$ensembl_gene_id, 
                                                        mart = ensembl.connected))


# saving transcript files to desktop
# xage1a
write_xlsx(xage1a_72_upregulated_transcripts, "xage1a_72_upregulated_transcripts.xlsx")
write_xlsx(xage1a_72_downregulated_transcripts, "xage1a_72_downregulated_transcripts.xlsx")
write_xlsx(xage1a_120_upregulated_transcripts, "xage1a_120_upregulated_transcripts.xlsx")
write_xlsx(xage1a_120_downregulated_transcripts, "xage1a_120_downregulated_transcripts.xlsx")

# xage2
write_xlsx(xage2_72_upregulated_transcripts, "xage2_72_upregulated_transcripts.xlsx")
write_xlsx(xage2_72_downregulated_transcripts, "xage2_72_downregulated_transcripts.xlsx")
write_xlsx(xage2_120_upregulated_transcripts, "xage2_120_upregulated_transcripts.xlsx")
write_xlsx(xage2_120_downregulated_transcripts, "xage2_120_downregulated_transcripts.xlsx")


######################################################################################################################
#                   CREATING TRANSCRIPT DATA FILES FOR TRANSCRIPTION FACTOR (BEDTOOLS) ANALYSIS
######################################################################################################################

# number of rows in the transcript files
# xage1a
nrow(xage1a_72_upregulated_transcripts)       # 104
nrow(xage1a_72_downregulated_transcripts)     # 213
nrow(xage1a_120_upregulated_transcripts)      # 95
nrow(xage1a_120_downregulated_transcripts)    # 38

# xage2 
nrow(xage2_72_upregulated_transcripts)        # 266
nrow(xage2_72_downregulated_transcripts)      # 230
nrow(xage2_120_upregulated_transcripts)       # 3220
nrow(xage2_120_downregulated_transcripts)     # 535

# assigning files to a variable, obtaining the promoter region and saving to a bed file
# test this by creating a text file and observing the columns to ensure data is as expected

# xage1a
t <- xage1a_72_upregulated_transcripts
for ( k in 1:104){ if(t[k,5] == "1") 
{ TSS = t[k,3]}else{ TSS = t[k,4]}; 
  write( paste(c(t[k,2], TSS-1000, TSS+1000, t[k,1]), collapse = "\t"),
         file = "xage1a_72_upregulated_promoter_region.bed", append = TRUE, sep = "\n") }

t <- xage1a_72_downregulated_transcripts
for ( k in 1:213){ if(t[k,5] == "1") 
{ TSS = t[k,3]}else{ TSS = t[k,4]}; 
  write( paste(c(t[k,2], TSS-1000, TSS+1000, t[k,1]), collapse = "\t"),
         file = "xage1a_72_downregulated_promoter_region.bed", append = TRUE, sep = "\n") }

t <- xage1a_120_upregulated_transcripts
for ( k in 1:95){ if(t[k,5] == "1") 
{ TSS = t[k,3]}else{ TSS = t[k,4]}; 
  write( paste(c(t[k,2], TSS-1000, TSS+1000, t[k,1]), collapse = "\t"),
         file = "xage1a_120_upregulated_promoter_region.bed", append = TRUE, sep = "\n") }

t <- xage1a_120_downregulated_transcripts
for ( k in 1:38){ if(t[k,5] == "1") 
{ TSS = t[k,3]}else{ TSS = t[k,4]}; 
  write( paste(c(t[k,2], TSS-1000, TSS+1000, t[k,1]), collapse = "\t"),
         file = "xage1a_120_downregulated_promoter_region.bed", append = TRUE, sep = "\n") }


# xage2
t <- xage2_72_upregulated_transcripts
for ( k in 1:266){ if(t[k,5] == "1") 
{ TSS = t[k,3]}else{ TSS = t[k,4]}; 
  write( paste(c(t[k,2], TSS-1000, TSS+1000, t[k,1]), collapse = "\t"),
         file = "xage2_72_upregulated_promoter_region.bed", append = TRUE, sep = "\n") }

t <- xage2_72_downregulated_transcripts
for ( k in 1:230){ if(t[k,5] == "1") 
{ TSS = t[k,3]}else{ TSS = t[k,4]}; 
  write( paste(c(t[k,2], TSS-1000, TSS+1000, t[k,1]), collapse = "\t"),
         file = "xage2_72_downregulated_promoter_region.bed", append = TRUE, sep = "\n") }

t <- xage2_120_upregulated_transcripts
for ( k in 1:3220){ if(t[k,5] == "1") 
{ TSS = t[k,3]}else{ TSS = t[k,4]}; 
  write( paste(c(t[k,2], TSS-1000, TSS+1000, t[k,1]), collapse = "\t"),
         file = "xage2_120_upregulated_promoter_region.bed", append = TRUE, sep = "\n") }

t <- xage2_120_downregulated_transcripts
for ( k in 1:535){ if(t[k,5] == "1") 
{ TSS = t[k,3]}else{ TSS = t[k,4]}; 
  write( paste(c(t[k,2], TSS-1000, TSS+1000, t[k,1]), collapse = "\t"),
         file = "xage2_120_downregulated_promoter_region.bed", append = TRUE, sep = "\n") }



######################################################################################################################
#                IDENTIFYING THE FOXP4 BINDING SITES IN THE DIFFERENTIALLY EXPRESSED GENES
######################################################################################################################

# using bedtools to determine if the promoter regions of differentially expressed genes align with ChipSeq data of FOXP4
# using terminal (linux), compare FOXP4 binding sites file to the xage2 promoter region files
# word file "XAGE Analysis in Terminal" contains the steps and code used


######################################################################################################################
#                         CONVERTING BEDTOOLS RESULTS FROM TRANSCRIPTS TO GENE IDS
######################################################################################################################

# DIFFERENTIALLY EXPRESSED GENES THAT BIND FOXP4
# reading files into R (all files are xage2)
xage2_72_upregulated_foxp4_bound_transcripts <- read.table("xage2_72_upregulated_foxp4_bound_transcripts.txt", header = FALSE, sep = "\t")
colnames(xage2_72_upregulated_foxp4_bound_transcripts)[1:4] <- c("Chromosome_Name", "Transcription_Start", "Transcription_End", "Ensembl_Transcript_ID")

xage2_72_downregulated_foxp4_bound_transcripts <- read.table("xage2_72_downregulated_foxp4_bound_transcripts.txt", header = FALSE, sep = "\t")
colnames(xage2_72_downregulated_foxp4_bound_transcripts)[1:4] <- c("Chromosome_Name", "Transcription_Start", "Transcription_End", "Ensembl_Transcript_ID")

xage2_120_upregulated_foxp4_bound_transcripts <- read.table("xage2_120_upregulated_foxp4_bound_transcripts.txt", header = FALSE, sep = "\t")
colnames(xage2_120_upregulated_foxp4_bound_transcripts)[1:4] <- c("Chromosome_Name", "Transcription_Start", "Transcription_End", "Ensembl_Transcript_ID")

xage2_120_downregulated_foxp4_bound_transcripts <- read.table("xage2_120_downregulated_foxp4_bound_transcripts.txt", header = FALSE, sep = "\t")
colnames(xage2_120_downregulated_foxp4_bound_transcripts)[1:4] <- c("Chromosome_Name", "Transcription_Start", "Transcription_End", "Ensembl_Transcript_ID")

# converting transcript ids into gene data tables
xage2_72_upregulated_foxp4_bound_gene_list <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description"),  
                                                               filters = "ensembl_transcript_id_version", 
                                                               values = xage2_72_upregulated_foxp4_bound_transcripts$Ensembl_Transcript_ID, 
                                                               mart = ensembl.connected))

xage2_72_downregulated_foxp4_bound_gene_list <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description"),  
                                                                 filters = "ensembl_transcript_id_version", 
                                                                 values = xage2_72_downregulated_foxp4_bound_transcripts$Ensembl_Transcript_ID, 
                                                                 mart = ensembl.connected))

xage2_120_upregulated_foxp4_bound_gene_list <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description"),  
                                                                filters = "ensembl_transcript_id_version", 
                                                                values = xage2_120_upregulated_foxp4_bound_transcripts$Ensembl_Transcript_ID, 
                                                                mart = ensembl.connected))

xage2_120_downregulated_foxp4_bound_gene_list <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description"),  
                                                                  filters = "ensembl_transcript_id_version", 
                                                                  values = xage2_120_downregulated_foxp4_bound_transcripts$Ensembl_Transcript_ID, 
                                                                  mart = ensembl.connected))

# uploading the new data tables to research folder on desktop
write_xlsx(xage2_72_upregulated_foxp4_bound_gene_list, "xage2_72_upregulated_foxp4_bound_gene_list.xlsx")
write_xlsx(xage2_72_downregulated_foxp4_bound_gene_list, "xage2_72_downregulated_foxp4_bound_gene_list.xlsx")
write_xlsx(xage2_120_upregulated_foxp4_bound_gene_list, "xage2_120_upregulated_foxp4_bound_gene_list.xlsx")
write_xlsx(xage2_120_downregulated_foxp4_bound_gene_list, "xage2_120_downregulated_foxp4_bound_gene_list.xlsx")



# DIFFERENTIALLY EXPRESSED GENES THAT DO NOT BIND FOXP4
# reading files into R (all files are xage2)
xage2_72_upregulated_foxp4_unbound_transcripts <- read.table("xage2_72_upregulated_foxp4_unbound_transcripts.txt", header = FALSE, sep = "\t")
colnames(xage2_72_upregulated_foxp4_unbound_transcripts)[1:4] <- c("Chromosome_Name", "Transcription_Start", "Transcription_End", "Ensembl_Transcript_ID")

xage2_72_downregulated_foxp4_unbound_transcripts <- read.table("xage2_72_downregulated_foxp4_unbound_transcripts.txt", header = FALSE, sep = "\t")
colnames(xage2_72_downregulated_foxp4_unbound_transcripts)[1:4] <- c("Chromosome_Name", "Transcription_Start", "Transcription_End", "Ensembl_Transcript_ID")

xage2_120_upregulated_foxp4_unbound_transcripts <- read.table("xage2_120_upregulated_foxp4_unbound_transcripts.txt", header = FALSE, sep = "\t")
colnames(xage2_120_upregulated_foxp4_unbound_transcripts)[1:4] <- c("Chromosome_Name", "Transcription_Start", "Transcription_End", "Ensembl_Transcript_ID")

xage2_120_downregulated_foxp4_unbound_transcripts <- read.table("xage2_120_downregulated_foxp4_unbound_transcripts.txt", header = FALSE, sep = "\t")
colnames(xage2_120_downregulated_foxp4_unbound_transcripts)[1:4] <- c("Chromosome_Name", "Transcription_Start", "Transcription_End", "Ensembl_Transcript_ID")

# converting transcript ids into gene data tables
xage2_72_upregulated_foxp4_unbound_gene_list <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description"),  
                                                                 filters = "ensembl_transcript_id_version", 
                                                                 values = xage2_72_upregulated_foxp4_unbound_transcripts$Ensembl_Transcript_ID, 
                                                                 mart = ensembl.connected))

xage2_72_downregulated_foxp4_unbound_gene_list <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description"),  
                                                                   filters = "ensembl_transcript_id_version", 
                                                                   values = xage2_72_downregulated_foxp4_unbound_transcripts$Ensembl_Transcript_ID, 
                                                                   mart = ensembl.connected))

xage2_120_upregulated_foxp4_unbound_gene_list <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description"),  
                                                                  filters = "ensembl_transcript_id_version", 
                                                                  values = xage2_120_upregulated_foxp4_unbound_transcripts$Ensembl_Transcript_ID, 
                                                                  mart = ensembl.connected))

xage2_120_downregulated_foxp4_unbound_gene_list <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description"),  
                                                                    filters = "ensembl_transcript_id_version", 
                                                                    values = xage2_120_downregulated_foxp4_unbound_transcripts$Ensembl_Transcript_ID, 
                                                                    mart = ensembl.connected))

# uploading the new data tables to research folder on desktop
write_xlsx(xage2_72_upregulated_foxp4_unbound_gene_list, "xage2_72_upregulated_foxp4_unbound_gene_list.xlsx")
write_xlsx(xage2_72_downregulated_foxp4_unbound_gene_list, "xage2_72_downregulated_foxp4_unbound_gene_list.xlsx")
write_xlsx(xage2_120_upregulated_foxp4_unbound_gene_list, "xage2_120_upregulated_foxp4_unbound_gene_list.xlsx")
write_xlsx(xage2_120_downregulated_foxp4_unbound_gene_list, "xage2_120_downregulated_foxp4_unbound_gene_list.xlsx")



# NON-DIFFERENTIALLY EXPRESSED GENES THAT DO NOT BIND FOXP4
# reading files into R (all files are xage2)
xage2_72_upregulated_nonbinding_nonDE_transcripts <- read.table("xage2_72_upregulated_nonbinding_nonDE_transcripts.txt", header = FALSE, sep = "\t")
colnames(xage2_72_upregulated_nonbinding_nonDE_transcripts)[1:4] <- c("Chromosome_Name", "Transcription_Start", "Transcription_End", "Ensembl_Transcript_ID")

xage2_72_downregulated_nonbinding_nonDE_transcripts <- read.table("xage2_72_downregulated_nonbinding_nonDE_transcripts.txt", header = FALSE, sep = "\t")
colnames(xage2_72_downregulated_nonbinding_nonDE_transcripts)[1:4] <- c("Chromosome_Name", "Transcription_Start", "Transcription_End", "Ensembl_Transcript_ID")

xage2_120_upregulated_nonbinding_nonDE_transcripts <- read.table("xage2_120_upregulated_nonbinding_nonDE_transcripts.txt", header = FALSE, sep = "\t")
colnames(xage2_120_upregulated_nonbinding_nonDE_transcripts)[1:4] <- c("Chromosome_Name", "Transcription_Start", "Transcription_End", "Ensembl_Transcript_ID")

xage2_120_downregulated_nonbinding_nonDE_transcripts <- read.table("xage2_120_downregulated_nonbinding_nonDE_transcripts.txt", header = FALSE, sep = "\t")
colnames(xage2_120_downregulated_nonbinding_nonDE_transcripts)[1:4] <- c("Chromosome_Name", "Transcription_Start", "Transcription_End", "Ensembl_Transcript_ID")

# converting transcript ids into gene data tables
xage2_72_upregulated_nonbinding_nonDE_gene_list <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description"),  
                                                                    filters = "ensembl_transcript_id_version", 
                                                                    values = xage2_72_upregulated_nonbinding_nonDE_transcripts$Ensembl_Transcript_ID, 
                                                                    mart = ensembl.connected))

xage2_72_downregulated_nonbinding_nonDE_gene_list <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description"),  
                                                                      filters = "ensembl_transcript_id_version", 
                                                                      values = xage2_72_downregulated_nonbinding_nonDE_transcripts$Ensembl_Transcript_ID, 
                                                                      mart = ensembl.connected))

xage2_120_upregulated_nonbinding_nonDE_gene_list <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description"),  
                                                                     filters = "ensembl_transcript_id_version", 
                                                                     values = xage2_120_upregulated_nonbinding_nonDE_transcripts$Ensembl_Transcript_ID, 
                                                                     mart = ensembl.connected))

xage2_120_downregulated_nonbinding_nonDE_gene_list <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description"),  
                                                                       filters = "ensembl_transcript_id_version", 
                                                                       values = xage2_120_downregulated_nonbinding_nonDE_transcripts$Ensembl_Transcript_ID, 
                                                                       mart = ensembl.connected))

# uploading the new data tables to research folder on desktop
write_xlsx(xage2_72_upregulated_nonbinding_nonDE_gene_list, "xage2_72_upregulated_nonbinding_nonDE_gene_list.xlsx")
write_xlsx(xage2_72_downregulated_nonbinding_nonDE_gene_list, "xage2_72_downregulated_nonbinding_nonDE_gene_list.xlsx")
write_xlsx(xage2_120_upregulated_nonbinding_nonDE_gene_list, "xage2_120_upregulated_nonbinding_nonDE_gene_list.xlsx")
write_xlsx(xage2_120_downregulated_nonbinding_nonDE_gene_list, "xage2_120_downregulated_nonbinding_nonDE_gene_list.xlsx")



# NON-DIFFERENTIALLY EXPRESSED GENES THAT BIND FOXP4
# reading files into R (all files are xage2)
xage2_72_upregulated_nonDE_foxp4_bound_transcripts <- read.table("xage2_72_upregulated_nonDE_foxp4_bound_transcripts.txt", header = FALSE, sep = "\t")
colnames(xage2_72_upregulated_nonDE_foxp4_bound_transcripts)[1:4] <- c("Chromosome_Name", "Transcription_Start", "Transcription_End", "Ensembl_Transcript_ID")

xage2_72_downregulated_nonDE_foxp4_bound_transcripts <- read.table("xage2_72_downregulated_nonDE_foxp4_bound_transcripts.txt", header = FALSE, sep = "\t")
colnames(xage2_72_downregulated_nonDE_foxp4_bound_transcripts)[1:4] <- c("Chromosome_Name", "Transcription_Start", "Transcription_End", "Ensembl_Transcript_ID")

xage2_120_upregulated_nonDE_foxp4_bound_transcripts <- read.table("xage2_120_upregulated_nonDE_foxp4_bound_transcripts.txt", header = FALSE, sep = "\t")
colnames(xage2_120_upregulated_nonDE_foxp4_bound_transcripts)[1:4] <- c("Chromosome_Name", "Transcription_Start", "Transcription_End", "Ensembl_Transcript_ID")

xage2_120_downregulated_nonDE_foxp4_bound_transcripts <- read.table("xage2_120_downregulated_nonDE_foxp4_bound_transcripts.txt", header = FALSE, sep = "\t")
colnames(xage2_120_downregulated_nonDE_foxp4_bound_transcripts)[1:4] <- c("Chromosome_Name", "Transcription_Start", "Transcription_End", "Ensembl_Transcript_ID")

# converting transcript ids into gene data tables
xage2_72_upregulated_nonDE_foxp4_bound_gene_list <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description"),  
                                                                     filters = "ensembl_transcript_id_version", 
                                                                     values = xage2_72_upregulated_nonDE_foxp4_bound_transcripts$Ensembl_Transcript_ID, 
                                                                     mart = ensembl.connected))

xage2_72_downregulated_nonDE_foxp4_bound_gene_list <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description"),  
                                                                       filters = "ensembl_transcript_id_version", 
                                                                       values = xage2_72_downregulated_nonDE_foxp4_bound_transcripts$Ensembl_Transcript_ID, 
                                                                       mart = ensembl.connected))

xage2_120_upregulated_nonDE_foxp4_bound_gene_list <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description"),  
                                                                      filters = "ensembl_transcript_id_version", 
                                                                      values = xage2_120_upregulated_nonDE_foxp4_bound_transcripts$Ensembl_Transcript_ID, 
                                                                      mart = ensembl.connected))

xage2_120_downregulated_nonDE_foxp4_bound_gene_list <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description"),  
                                                                        filters = "ensembl_transcript_id_version", 
                                                                        values = xage2_120_downregulated_nonDE_foxp4_bound_transcripts$Ensembl_Transcript_ID, 
                                                                        mart = ensembl.connected))

# uploading the new data tables to research folder on desktop
write_xlsx(xage2_72_upregulated_nonDE_foxp4_bound_gene_list, "xage2_72_upregulated_nonDE_foxp4_bound_gene_list.xlsx")
write_xlsx(xage2_72_downregulated_nonDE_foxp4_bound_gene_list, "xage2_72_downregulated_nonDE_foxp4_bound_gene_list.xlsx")
write_xlsx(xage2_120_upregulated_nonDE_foxp4_bound_gene_list, "xage2_120_upregulated_nonDE_foxp4_bound_gene_list.xlsx")
write_xlsx(xage2_120_downregulated_nonDE_foxp4_bound_gene_list, "xage2_120_downregulated_nonDE_foxp4_bound_gene_list.xlsx")


######################################################################################################################
#                                            FISHERS EXACT TEST
######################################################################################################################

# converting transcripts from human genome file to genes
human_genome_gene_list <- data.frame(getBM(attributes = c("ensembl_gene_id", "external_gene_name"),  
                                           filters = "ensembl_transcript_id_version", 
                                           values = human_genome$V1, 
                                           mart = ensembl.connected))

# uploading the new data tables to research folder on desktop
write_xlsx(human_genome_gene_list, "human_genome_gene_list.xlsx")

# obtaining the number of genes in human genome
nrow(unique(human_genome_gene_list))                                # 68020

# obtaining the number of differentially expressed genes
nrow(xage2_72_upregulated_gene_list)                                # 93
nrow(xage2_72_downregulated_gene_list)                              # 72
nrow(xage2_120_upregulated_gene_list)                               # 662
nrow(xage2_120_downregulated_gene_list)                             # 140

# obtaining the number of genes in each condition
# differentially expressed genes bound by foxp4
nrow(unique(xage2_72_upregulated_foxp4_bound_gene_list))            # 20
nrow(unique(xage2_72_downregulated_foxp4_bound_gene_list))          # 7
nrow(unique(xage2_120_upregulated_foxp4_bound_gene_list))           # 100
nrow(unique(xage2_120_downregulated_foxp4_bound_gene_list))         # 44

# differentially expressed genes not bound by foxp4
nrow(unique(xage2_72_upregulated_foxp4_unbound_gene_list))          # 74
nrow(unique(xage2_72_downregulated_foxp4_unbound_gene_list))        # 66
nrow(unique(xage2_120_upregulated_foxp4_unbound_gene_list))         # 602
nrow(unique(xage2_120_downregulated_foxp4_unbound_gene_list))       # 115

# non-differentially expressed genes not bound by foxp4
nrow(unique(xage2_72_upregulated_nonbinding_nonDE_gene_list))       # 25593
nrow(unique(xage2_72_downregulated_nonbinding_nonDE_gene_list))     # 25587
nrow(unique(xage2_120_upregulated_nonbinding_nonDE_gene_list))      # 25326
nrow(unique(xage2_120_downregulated_nonbinding_nonDE_gene_list))    # 25557

# non-differentially expressed genes that bind foxp4
nrow(unique(xage2_72_upregulated_nonDE_foxp4_bound_gene_list))      # 8967
nrow(unique(xage2_72_downregulated_nonDE_foxp4_bound_gene_list))    # 9025
nrow(unique(xage2_120_upregulated_nonDE_foxp4_bound_gene_list))     # 8861
nrow(unique(xage2_120_downregulated_nonDE_foxp4_bound_gene_list))   # 8929

nrow(unique(rbind(xage2_72_upregulated_nonDE_foxp4_bound_gene_list, xage2_72_downregulated_nonDE_foxp4_bound_gene_list)))
nrow(unique(rbind(xage2_72_upregulated_nonbinding_nonDE_gene_list, xage2_72_downregulated_nonbinding_nonDE_gene_list)))

nrow(unique(rbind(xage2_120_upregulated_nonDE_foxp4_bound_gene_list, xage2_120_downregulated_nonDE_foxp4_bound_gene_list)))
nrow(unique(rbind(xage2_120_upregulated_nonbinding_nonDE_gene_list, xage2_120_downregulated_nonbinding_nonDE_gene_list)))

# creating dataframes for fishers test
xage2_72_fishers = matrix(c(27, 140, 9037, 25619), nrow = 2, dimnames = list(c("FOXP4 Bound", "FOXP4 Unbound"), c("Differentially Expressed", "Non-differentially Expressed")))
xage2_72_fishers
# p-value = 0.003368

xage2_120_fishers = matrix(c(144, 717, 8920, 25042), nrow = 2, dimnames = list(c("FOXP4 Bound", "FOXP4 Unbound"), c("Differentially Expressed", "Non-differentially Expressed")))
xage2_120_fishers
fisher.test(xage2_120_fishers)
# p-value = 5.468e-11