### Check expression in mouse brain tissue and get gene ontologies for each gene ###

# Load the required libraries
required_packages <- c("openxlsx",  "reshape2", 
                       "dplyr", "magrittr", "parallel")

# Install missing packages
for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg, dependencies = TRUE)
    }
}

if (!requireNamespace("clusterProfiler", quietly = TRUE)) {Biocmanager::install("clusterProfiler") }
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {Biocmanager::install("org.Hs.eg.db") }
if (!requireNamespace("clustifyr", quietly = TRUE)) {remotes::install_github("rnabioco/clustifyR") }


# Load the libraries

suppressWarnings(suppressMessages(library(openxlsx)))
suppressWarnings(suppressMessages(library(reshape2)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(magrittr)))
suppressWarnings(suppressMessages(library(clusterProfiler)))
suppressWarnings(suppressMessages(library(clustifyR)))
suppressWarnings(suppressMessages(library(org.Hs.eg.db)))
suppressWarnings(suppressMessages(library(parallel)))

source('R/utils.R')

# Load the latest results file
genes <- read.xlsx("results/TRIFF_risk_genes_info.xlsx")

## Add function/roles of each gene to help selection of genes by adding top GO terms to the genes dataframe
                 
genes$TopGO_BP <- NA
genes$TopGO_CC <- NA
genes$TopGO_MF <- NA

# Extract top GO terms for every gene in the dataframe. Use mclapply to parallelize the process. Takes a long time to run
go_terms_list <- mclapply(genes$Approved.symbol, function(symbol) {
    if (!is.na(symbol)) {
        get_top_go_terms(symbol)
    } else {
        list(BP = NA, CC = NA, MF = NA)
    }
}, mc.cores = detectCores() - 1) # ran this on the server with 256 cores; needs >500gb RAM for some reason. 

# Extract the results and assign them to the genes dataframe
genes$TopGO_BP <- sapply(go_terms_list, function(x) x$BP)
genes$TopGO_CC <- sapply(go_terms_list, function(x) x$CC)
genes$TopGO_MF <- sapply(go_terms_list, function(x) x$MF)

#write.xlsx(genes, "results/TRIFF_risk_genes_info_with_GO.xlsx")

                         
# Load the Linnarsson dataset
lin <- get_ucsc_reference(cb_url = "https://cells.ucsc.edu/?ds=mouse-dev-brain",cluster_col = "Class",if_log = FALSE)
rownames(lin) = gsub(".+[|]", "", rownames(lin))

# Check if genes$mouse_gene_name is expressed in any of the celltypes and add the cell types if expression is > 0.001
genes$MouseBrainExpressed <- NA
genes$MouseCTExpressed <- NA

for (i in seq_len(nrow(genes))) {
    if (genes$mouse_gene_name[i] %in% rownames(lin)) {
        gene_expression <- lin[genes$mouse_gene_name[i], ]
        if (any(gene_expression > 0.001)) {
            genes$MouseBrainExpressed[i] <- TRUE
            genes$MouseCTExpressed[i] <- paste(colnames(lin)[gene_expression > 0.001], collapse = ";")
        } else {
            genes$MouseBrainExpressed[i] <- FALSE
        }
    } else {
        genes$MouseBrainExpressed[i] <- FALSE
    }
}

# Write the results to a file
write.xlsx(genes, "results/TRIFF_risk_genes_info_with_Ms_and_Ontologies.xlsx")



# Count how many genes are BrainExpressed == TRUE, MouseBrainExpressed == TRUE, 
# and are either prenatal or postnatal enriched in genes$Human_Period_Enrichment
count_prenatal_enriched <- sum(genes$BrainExpressed == TRUE & 
                               genes$MouseBrainExpressed == TRUE & 
                               grepl("prenatal", genes$Human_Period_Enrichment, ignore.case = TRUE) & 
                               genes$MouseBrainExpressed == TRUE, na.rm = TRUE)
#556
count_postnatal_enriched <- sum(genes$BrainExpressed == TRUE & 
                                genes$MouseBrainExpressed == TRUE & 
                                grepl("postnatal", genes$Human_Period_Enrichment, ignore.case = TRUE) & 
                                genes$MouseBrainExpressed == TRUE, na.rm = TRUE)
#596

count_brain_expressed <- sum(genes$BrainExpressed == TRUE & genes$MouseBrainExpressed == TRUE, na.rm = TRUE) 
#1239

count_brain_not_mouse <- sum(genes$BrainExpressed == TRUE & genes$MouseBrainExpressed == FALSE, na.rm = TRUE)
#188

sum(genes$BrainExpressed == FALSE & genes$MouseBrainExpressed == TRUE, na.rm = TRUE)
#178
