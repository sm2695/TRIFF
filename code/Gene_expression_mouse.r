### Still working on this 

# Load the required libraries
suppressWarnings(suppressMessages(library(Seurat)))
suppressWarnings(suppressMessages(library(SeuratExtend)))
suppressWarnings(suppressMessages(library(openxlsx)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(reshape2)))
suppressWarnings(suppressMessages(library(viridis)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(magrittr)))
suppressWarnings(suppressMessages(library(clusterProfiler)))
suppressWarnings(suppressMessages(library(org.Hs.eg.db)))
suppressWarnings(suppressMessages(library(parallel)))

source('R/utils.R')

# Load the latest results file
genes <- read.xlsx("/Volumes/projects/C3_Sellgren_lab/Lab Members/Susmita/Internal data/SCZ Meta-Analysis_triff/TRIFF_main/results/TRIFF_risk_genes_info.xlsx")

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
}, mc.cores = detectCores() - 1) # run on ther server with >500gb RAM for some reason. 

# Extract the results and assign them to the genes dataframe
genes$TopGO_BP <- sapply(go_terms_list, function(x) x$BP)
genes$TopGO_CC <- sapply(go_terms_list, function(x) x$CC)
genes$TopGO_MF <- sapply(go_terms_list, function(x) x$MF)

write.xlsx(genes, "/Volumes/projects/C3_Sellgren_lab/Lab Members/Susmita/Internal data/SCZ Meta-Analysis_triff/TRIFF_main/results/TRIFF_risk_genes_info_with_GO.xlsx")


# Load the Linnarsson dataset
lin <- readRDS("/Volumes/projects/C3_Sellgren_lab/Lab Members/Susmita/External datasets/Single Cell Studies/Mouse datasets/Linnarsson_mousedev.rds")
lin

# Check if genes$mouse_gene_name is expressed in lin and add the colnames from tab if expression is > 0.1
genes_expressed <- genes$mouse_gene_name[genes$mouse_gene_name %in% rownames(tab)[rowSums(tab > 0.1) > 0]]
colnames_expressed <- colnames(tab)[apply(tab[genes_expressed, ], 2, function(x) any(x > 0.1))]

tab <- AggregateExpression(lin, group.by = c("Region","Class"))$RNA

tab[1:5,1:5]


filtered_tab <- tab[rownames(tab) %in% genes$mouse_gene_name, ]
filtered_tab[1:5, 1:5]

names_to_df <- colnames(filtered_tab)

names_df <- data.frame(
    region = sapply(strsplit(names_to_df, "_"), `[`, 1),
    celltype = sapply(strsplit(names_to_df, "_"), `[`, 2)
)

head(names_df)
names_df$region <- gsub("b|'","", names_df$region)
colnames(tab) <- paste(names_df$region, names_df$celltype, sep = "_")

genes$MouseBrainExpressed <- NA
genes$MouseCTExpressed <- NA
for(i in seq_len(nrow(genes))){
    if(genes$mouse_gene_name[i] %in% rownames(filtered_tab)){
        genes$MouseBrainExpressed[i] <- TRUE
        if(genes$MouseBrainExpressed[i]){
            genes$MouseCTExpressed[i] <- paste(colnames(filtered_tab)[filtered_tab[genes$mouse_gene_name[i], ] > 0.1], collapse = ";")
        }
    }
}
View(genes)

write.xlsx(genes, "/Volumes/projects/C3_Sellgren_lab/Lab Members/Susmita/Internal data/SCZ Meta-Analysis_triff/TRIFF_main/results/TRIFF_risk_genes_info_with_mouse_expression.xlsx")
