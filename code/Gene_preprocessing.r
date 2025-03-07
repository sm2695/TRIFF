
#### TRIFF Project #####
# Code to curate all genes associated to the risk of developing schizophrenia. The following sources are considered:
#A. common variants
# 1. GWAS (Trubetskoy et al. 2021)
# 2. HiC-defined SCZ risk genes from GWS (https://doi.org/10.1016/j.schres.2019.03.007)
# 3. MPRA-HiC (https://doi.org/10.1016/j.xgen.2023.100404)
#B. rare variants
# 4. CNV
# 5. SCHEMA (Singh et al. 2022)
#C. differentially expressed genes from postmortem comparisons
# 5. Postmortem-Batiuk
# 6. Postmortem-Ruzicka

source("R/utils.R")
rm(list = ls())
options(stringsAsFactors = FALSE)

if (!requireNamespace("rentrez", quietly = TRUE)) install.packages("rentrez")
if (!requireNamespace("XML", quietly = TRUE)) install.packages("XML")
if (!requireNamespace("biomaRt", quietly = TRUE)) install.packages("biomaRt")
if (!requireNamespace("httr", quietly = TRUE)) install.packages("httr")
if (!requireNamespace("jsonlite", quietly = TRUE)) install.packages("jsonlite")
if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")


# Load the required libraries
library(httr)
library(jsonlite)
library(rentrez)
library(XML)
library(readxl)
library(openxlsx)
library(biomaRt)

#### Get Common variant genes from the latest GWAS #######
gwas <- read.xlsx("data/Trubetskoy_SCZ_GWAS.xlsx")
gwas_genes <- unique(gwas$Symbol.ID[which(gwas$FINEMAPk3.5 == 1)]) ## FINEMAP credible SNPS fromTrubetskoy et al. 2021
length(gwas_genes)
#625

#gwas_genes <- ifelse(grepl("\\.", gwas_genes), sub("\\..*", "", gwas_genes), gwas_genes)

# Fix bad names
gwas_genes[grepl("ago", gwas_genes)] <- c("AGO1", "AGO3", "AGO4")
#gwas_genes[grepl("ENSG", gwas_genes)] ## Novel RNA genes. take care of these later if they persist

# Read the supp file https://doi.org/10.1016/j.schres.2019.03.007
hic_data <- read.table("data/HiC_defined_SCZ_risk_genes_from_GWS.txt", header = TRUE, sep = "\t")

# 
# Assuming we need to extract a specific column, e.g., 'Gene'
hic_genes <- unique(hic_data[[1]])
length(hic_genes)
#455

mpra <- read.xlsx("data/MPRA-HiC_SCZ.xlsx")
colnames(mpra) <- mpra[1, ]
mpra <- mpra[-1, ]
mpra_genes <- unique(mpra$HGNC_symbol)
length(mpra_genes)
#272

###### Get rare variant genes ######

schema <- c("SET1DA", "CUL1", "XPO7", "TRIO", "CACNA1G", "SP4", "GRIA3", "GRIN2A", "HERC1", "RB1CC1")
length(schema)
#10

all_cnv<- readRDS("data/Genes_with_CNV.rds")

# include 22q11.2del, PWS/AS dup, 3q29 del, NRXN1 del(2p16.3), 16p11.2 dup, 1q21.1. del, WBS dup, 15q13.3del, 16p12.1del, 1q21.1dup, 16p13.11dup, and 15q11.2del

keep <- c("22q11.21_DUP/DEL","2p16.3_DUP/DEL","7q11.23_DUP/DEL", "3q29_DUP/DEL", "16p11.2_DUP/DEL-A","16p11.2_DUP/DEL-B" ,"1q21.1-q21.2_DUP/DEL","15q11.2-q13.3_DUP", "15q13.2-q13.3_DEL","16p12.2_DUP/DEL","16p13.11_DUP/DEL")

cnv_genes <- unique(all_cnv[which(all_cnv$CNV %in% keep),1])
length(cnv_genes)
#199


###### Get differentially expressed genes from Postmortem comparisons ######

# Load the data from Batiuk et al.
kkh_deg <- read.xlsx("data/Batiuk_DEGenes.xlsx", sheet = 1)
kkh_deg$padj <- as.numeric(kkh_deg$padj)
kkh_genes <- unique(kkh_deg$Gene[which(kkh_deg$padj < 0.05)])
# Remove genes starting with "MT-" from kkh_genes
kkh_genes <- kkh_genes[!grepl("^MT-", kkh_genes)]
length(kkh_genes)
#53

# Load the data from Euzicka et al.
ruz_deg <- read.xlsx("data/Ruzicka_DEGenes.xlsx", sheet = 1)
ruz_deg$Meta_adj.P.Val <- as.numeric(ruz_deg$Meta_adj.P.Val)
ruz_genes <- unique(ruz_deg$gene[which(ruz_deg$Meta_adj.P.Val < 0.05)])
length(ruz_genes)
#224

# Load selected genes from the twins organoid scz study
twins_genes <- read.xlsx("data/Sellgren_dev_PDEgenes.xlsx")
twins_genes <- unique(twins_genes$Gene)
length(twins_genes)
#168

######## Combine all gene vectors into a single data frame  ############
all_genes <- data.frame(
  Gene = c(gwas_genes, hic_genes,mpra_genes, cnv_genes, schema, kkh_genes, ruz_genes, twins_genes),
  Source = c(
    rep("GWAS", length(gwas_genes)),
    rep("HiC-defined", length(hic_genes)),
    rep("MPRA-HiC", length(mpra_genes)),
    rep("CNV", length(cnv_genes)),
    rep("SCHEMA" , length(schema)),
    rep("Postmortem-Batiuk", length(kkh_genes)),
    rep("Postmortem-Ruzicka", length(ruz_genes)),
    rep("Dev-Organoid", length(twins_genes))
  )
)

#View(all_genes)
# Update the Source column to include CNV information
all_genes$Source<- ifelse(all_genes$Source == "CNV", paste0(all_genes$Source,"_", all_cnv$CNV[match(all_genes$Gene, all_cnv[[1]])]), all_genes$Source)
    
# Remove duplicates and append sources
all_genes <- aggregate(Source ~ Gene, data = all_genes, FUN = function(x) paste(unique(x), collapse = ","))

# Connect to Ensembl (Human & Mouse)
ensembl_human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")#, host = "http://useast.ensembl.org/")
ensembl_mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")#,host = "http://useast.ensembl.org/")


# Fetch gene biotype information and mouse homologs for all genes
biotype_info <- get_gene_info(all_genes$Gene)
#View(biotype_info)

# Merge the biotype information with the all_genes data frame
all_genes <- merge(all_genes, biotype_info, by.x = "Gene", by.y = "hgnc_symbol", all.x = TRUE)
#View(all_genes)


#Define the list of odd gene names
odd_gene_names <- all_genes[which(is.na(all_genes$external_gene_name)),1]
#Query for all genes and store the results
odd_gene_data_list <- lapply(odd_gene_names, query_exphewas)
#combine the list of results into a data frame
odd_gene_data_df <- do.call(rbind, odd_gene_data_list)
odd_gene_data_df <- as.data.frame(odd_gene_data_df)
#View(odd_gene_data_df)

# Match the names in all_genes$Gene with odd_gene_data_df$name and update gene_biotype and description
for (i in 1:nrow(all_genes)) {
    match_index <- which(odd_gene_data_df$name == all_genes$Gene[i])
    if (length(match_index) > 0) {
        all_genes$gene_biotype[i] <- odd_gene_data_df$biotype[match_index]
        all_genes$description[i] <- odd_gene_data_df$description[match_index]
        all_genes$ensembl_gene_id[i] <- odd_gene_data_df$ensembl_id[match_index]
    }
}

#View(all_genes)

# Create directory for BrainSpan if it doesn't exist
if (!dir.exists("results")) {
    dir.create("results", recursive = TRUE)
}
                       
# Write the second-last gene list to a file
all_genes[] <- lapply(all_genes, function(x) if (is.list(x)) sapply(x, toString) else x)
write.xlsx(all_genes, "results/TRIFF_genes_unupdated_tmp.xlsx")

### Run the gene names manually against https://www.genenames.org/tools/multi-symbol-checker/ and add a column with approved HGNC symbols

# Tried running it via httr here but the code fails consistently due the the website's API. 

# Make sure to save the file as an xlsx file before reading the file in. csv files give an error

# read in the results
symbol_check <- read.xlsx("results/hgnc-symbol-check.xlsx")

# Merge the symbol_check data with all_genes
all_genes2 <- merge(all_genes, symbol_check, by.x = "Gene", by.y = "Input", all.x = TRUE)
#View(all_genes2)

# Convert list columns to character
all_genes2[] <- lapply(all_genes2, function(x) if (is.list(x)) sapply(x, toString) else x)

# Write the final gene list to a file
write.xlsx(all_genes2, file = "results/TRIFF_SCZ_curated_list_of_genes_updated.xlsx", row.names = FALSE, na = "")

