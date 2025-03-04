


if (!requireNamespace("rentrez", quietly = TRUE)) install.packages("rentrez")
if (!requireNamespace("XML", quietly = TRUE)) install.packages("XML")
if (!requireNamespace("biomaRt", quietly = TRUE)) install.packages("biomaRt")
# Load required libraries
library(rentrez)
library(XML)
library(readxl)
library(biomaRt)

#### Get Common variant genes from the latest GWAS #######
gwas <- readxl::read_xls("/Volumes/projects/C3_Sellgren_lab/Lab Members/Susmita/Internal data/SCZ Meta-Analysis/TRIFF/data/Trubetskoy_Supplementary_Table_3.xls")

#  Extract column with the 'SNP IDs'
top_index <- gwas[[2]]

# Fetching and Parsing SNP Data from NCBI's Entrez database
snp_data <- entrez_fetch(db = "snp", id = top_index, rettype = "xml", parsed = FALSE)

# Convert fetched data into an XML document
xml_doc <- xmlParse(snp_data)  # Explicitly parse since parsed = FALSE

# Define XML namespace for NCBI SNP records
ns <- c(ns0 = "https://www.ncbi.nlm.nih.gov/SNP/docsum")

# Extract gene names
gwas_genes <- xpathSApply(xml_doc, "//ns0:GENE_E/ns0:NAME", xmlValue, namespaces = ns)
length(gwas_genes)
#240

## FINEMAP credible SNPS fromTrubetskoy et al. 2021
trubetskoy <- readRDS("/Volumes/projects/C3_Sellgren_lab/Lab Members/Susmita/Internal data/C1Q/Bipolar/Trubetskoy_scz_genes.rds")
gwas_genes <- unique(trubetskoy$all_scz_trubetskoy)

#gwas_genes <- ifelse(grepl("\\.", gwas_genes), sub("\\..*", "", gwas_genes), gwas_genes)

# Fix bad names

gwas_genes[grepl("ago", gwas_genes)] <- c("AGO1", "AGO3","AGO4")

gwas_genes[grepl("ENSG", gwas_genes)] ## Novel RNA genes. take care of these later after file generation

# Read the supp file https://doi.org/10.1016/j.schres.2019.03.007
hic_data <- read.table("/Volumes/projects/C3_Sellgren_lab/Lab Members/Susmita/Internal data/SCZ Meta-Analysis/TRIFF/data/HiC_defined_SCZ_risk_genes_from_GWS.txt", header = TRUE, sep = "\t")

# 
# Assuming we need to extract a specific column, e.g., 'Gene'
hic_genes <- unique(hic_data[[1]])
length(hic_genes)

mpra <- read.xlsx("/Volumes/projects/C3_Sellgren_lab/Lab Members/Susmita/Internal data/SCZ Meta-Analysis/TRIFF/data/MPRA-HiC_SCZ.xlsx")
colnames(mpra) <- mpra[1,]
mpra <- mpra[-1,]
mpra_genes <- unique(mpra$HGNC_symbol)
length(mpra_genes)


###### Get rare variant genes ######

schema <- c("SET1DA", "CUL1", "XPO7", "TRIO", "CACNA1G", "SP4", "GRIA3", "GRIN2A", "HERC1", "RB1CC1")

all_cnv<- readRDS("/Volumes/projects/C3_Sellgren_lab/Lab Members/Susmita/Internal data/CNV/2023/Genes_with_CNV.rds")

# include 22q11.2del, PWS/AS dup, 3q29 del, NRXN1 del(2p16.3), 16p11.2 dup, 1q21.1. del, WBS dup, 15q13.3del, 16p12.1del, 1q21.1dup, 16p13.11dup, and 15q11.2del

keep <- c("22q11.21_DUP/DEL","2p16.3_DUP/DEL","7q11.23_DUP/DEL", "3q29_DUP/DEL", "16p11.2_DUP/DEL-A","16p11.2_DUP/DEL-B" ,"1q21.1-q21.2_DUP/DEL","15q11.2-q13.3_DUP", "15q13.2-q13.3_DEL","16p12.2_DUP/DEL","16p13.11_DUP/DEL")

cnv_genes <- unique(all_cnv[which(all_cnv$CNV %in% keep),1])
length(cnv_genes)
#199


###### Get differentially expressed genes from Postmortem comparisons ######

# Load the data from Batiuk et al.
kkh_deg <- read.xlsx("/Volumes/projects/C3_Sellgren_lab/Lab Members/Susmita/Internal data/SCZ Meta-Analysis/TRIFF/data/Kostya_DEGenes.xlsx", sheet = 1)
colnames(kkh_deg) <- kkh_deg[1,]
kkh_deg <- kkh_deg[-1,]
kkh_deg$padj <- as.numeric(kkh_deg$padj)
kkh_genes <- unique(kkh_deg$Gene[which(kkh_deg$padj < 0.05)])
# Remove genes starting with "MT-" from kkh_genes
kkh_genes <- kkh_genes[!grepl("^MT-", kkh_genes)]
length(kkh_genes)


# Load the data from Euzicka et al.
ruz_deg <- read.xlsx("/Volumes/projects/C3_Sellgren_lab/Lab Members/Susmita/Internal data/SCZ Meta-Analysis/TRIFF/data/Ruzicka_DEGenes.xlsx", sheet = 1)
ruz_deg$Meta_adj.P.Val <- as.numeric(ruz_deg$Meta_adj.P.Val)
ruz_genes <- unique(ruz_deg$gene[which(ruz_deg$Meta_adj.P.Val < 0.05)])
length(ruz_genes)


# Combine all gene vectors into a single data frame
all_genes <- data.frame(
    Gene = c(gwas_genes, hic_genes,mpra_genes, cnv_genes, kkh_genes, ruz_genes),
    Source = c(
        rep("GWAS", length(gwas_genes)),
        rep("HiC-defined", length(hic_genes)),
        rep("MPRA-HiC", length(mpra_genes)),
        rep("CNV", length(cnv_genes)),
        rep("Postmortem-Batiuk", length(kkh_genes)),
        rep("Postmortem-Ruzicka", length(ruz_genes))
    )
)

View(all_genes)
# Update the Source column to include CNV information
all_genes$Source<- ifelse(all_genes$Source == "CNV", paste0(all_genes$Source,"_", all_cnv$CNV[match(all_genes$Gene, all_cnv[[1]])]), all_genes$Source)
    
# Remove duplicates and append sources
all_genes <- aggregate(Source ~ Gene, data = all_genes, FUN = function(x) paste(unique(x), collapse = ","))

# Connect to Ensembl (Human & Mouse)
ensembl_human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ensembl_mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl",host = "http://useast.ensembl.org/")

# Function to fetch gene biotype information from Ensembl using biomaRt 
get_gene_info <- function(gene_symbols) {
    # STEP 1: Get Human Gene Data
    human_attributes <- c("hgnc_symbol", "external_gene_name", "gene_biotype", "description", "ensembl_gene_id")  
    human_data <- getBM(
        attributes = human_attributes, 
        filters = "hgnc_symbol", 
        values = gene_symbols, 
        mart = ensembl_human
    )
    
    # STEP 2: Get Mouse Homologs using the Mouse Dataset
    mouse_attributes <- c("external_gene_name", "hsapiens_homolog_associated_gene_name")
    mouse_homologs <- getBM(
        attributes = mouse_attributes, 
        mart = ensembl_mouse
    )
    
    # Rename for merging
    colnames(mouse_homologs) <- c("mouse_gene_name", "hgnc_symbol")
    
    # STEP 3: Merge Human Data with Mouse Homologs
    final_data <- merge(human_data, mouse_homologs, by = "hgnc_symbol", all.x = TRUE)
    
    return(final_data)
}

# Fetch gene biotype information and mouse homologs for all genes
biotype_info <- get_gene_info(all_genes$Gene)
View(biotype_info)

# Merge the biotype information with the all_genes data frame
all_genes <- merge(all_genes, biotype_info, by.x = "Gene", by.y = "hgnc_symbol", all.x = TRUE)
View(all_genes)







if (!requireNamespace("httr", quietly = TRUE)) install.packages("httr")
if (!requireNamespace("jsonlite", quietly = TRUE)) install.packages("jsonlite")

# Load the required libraries
library(httr)
library(jsonlite)

# Define the list of gene names
odd_gene_names <- all_genes[which(is.na(all_genes$external_gene_name)),1]

# Function to query exphewas API for a single gene
query_exphewas <- function(gene_name) {
  url <- paste0("https://exphewas.statgen.org/v1/api/gene/name/", gene_name)
  response <- GET(url)  # Send GET request to the API
  
  # Check if the request was successful
  if (status_code(response) == 200) {
    content_data <- content(response, "text", encoding = "UTF-8")
    parsed_data <- fromJSON(content_data, flatten = TRUE)
    return(parsed_data)
  } else {
    warning(paste("Failed to fetch data for gene:", gene_name))
    return(NULL)
  }
}

# Query for all genes and store the results
odd_gene_data_list <- lapply(odd_gene_names, query_exphewas)

# combine the list of results into a data frame
odd_gene_data_df <- do.call(rbind, odd_gene_data_list)
odd_gene_data_df <- as.data.frame(odd_gene_data_df)
View(odd_gene_data_df)

# Match the names in all_genes$Gene with odd_gene_data_df$name and update gene_biotype and description
for (i in 1:nrow(all_genes)) {
    match_index <- which(odd_gene_data_df$name == all_genes$Gene[i])
    if (length(match_index) > 0) {
        all_genes$gene_biotype[i] <- odd_gene_data_df$biotype[match_index]
        all_genes$description[i] <- odd_gene_data_df$description[match_index]
        all_genes$ensembl_gene_id[i] <- odd_gene_data_df$ensembl_id[match_index]
    }
}

View(all_genes)

# Fix some leftover gene names
all_genes[grepl("ENSG", all_genes$Gene),6] <- all_genes[grepl("ENSG", all_genes$Gene),1]
all_genes[grepl("ENSG", all_genes$Gene),3] <- c("", "RN7SL199P", "RN7SL656P", "RDM1P1", "")
all_genes[grepl("ENSG", all_genes$Gene),1] <- c("ENSG0000026231","RN7SL199P","RN7SL656P", "RDM1P1","ENSG00000278770")


# Write the second-last gene list to a file
write.xlsx(all_genes, "/Volumes/projects/C3_Sellgren_lab/Lab Members/Susmita/Internal data/SCZ Meta-Analysis/TRIFF/results/TRIFF_SCZ_curated_list_of_genes_temp.xlsx", na = "")


## Run the gene names against https://www.genenames.org/tools/multi-symbol-checker/ and add a column with approved HGNC symbols

# read in the results
symbol_check <- read.xlsx("/Volumes/projects/C3_Sellgren_lab/Lab Members/Susmita/Internal data/SCZ Meta-Analysis/TRIFF/results/hgnc-symbol-check.xlsx")

head(symbol_check)
symbol_check <- symbol_check[-1,]

# Merge the symbol_check data with all_genes
all_genes2 <- merge(all_genes, symbol_check, by.x = "Gene", by.y = "Input", all.x = TRUE)
View(all_genes2)


# Convert list columns to character
all_genes2[] <- lapply(all_genes2, function(x) if (is.list(x)) sapply(x, toString) else x)

# Write the final gene list to a file
na.string=""
write.xlsx(all_genes2, file = "/Volumes/projects/C3_Sellgren_lab/Lab Members/Susmita/Internal data/SCZ Meta-Analysis/TRIFF/results/TRIFF_SCZ_curated_list_of_genes_updated.xlsx", row.names = FALSE, na = "")
