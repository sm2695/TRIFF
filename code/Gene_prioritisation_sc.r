if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")
if (!requireNamespace("EWCE", quietly = TRUE)) install.packages("EWCE")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("reshape2", quietly = TRUE)) install.packages("reshape2")
if (!requireNamespace("viridis", quietly = TRUE)) install.packages("viridis")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("magrittr", quietly = TRUE)) install.packages("magrittr")

suppressWarnings(suppressMessages(library(openxlsx)))
suppressWarnings(suppressMessages(library(EWCE)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(reshape2)))
suppressWarnings(suppressMessages(library(viridis)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(magrittr)))

## Read in results from Brainspan 

genes<- read.xlsx("results/TRIFF_risk_gene_BrainSpan_info.xlsx")

## Download the CTD data for Human Fetal scRNAseq dataset from https://figshare.com/s/6cba3ae11e6c17335f4f
download.file("https://figshare.com/ndownloader/files/34065392", destfile = "data/ctd_FetalDT2023.rda")
download.file("https://figshare.com/ndownloader/files/34065395", destfile = "data/ctd_FetalCT2023.rda")

load("data/ctd_FetalDT2023.rda")

# Get pre-computed mean expression and specificity values for all genes.
exp <- ctd[[1]]$mean_exp
spec <- ctd[[1]]$specificity

exp <- exp[, colnames(exp) != "cerebellum"]
spec <- spec[, colnames(spec) != "cerebellum"]

# Define the cutoff (Mean + 1 SD)
nonzero_values <- exp[exp > 0]
mean_expr <- mean(nonzero_values)
std_expr <- sd(nonzero_values)
cutoff <- mean_expr + std_expr


# Check which genes are not in rownames of exp
genes_not_in_exp <- genes$Gene[!genes$Gene %in% rownames(exp)]
# Check if any of the genes not in exp have their Approved.Symbols in exp instead
approved_symbols <- genes$Approved.symbol[!genes$Gene %in% rownames(exp)]
genes_not_in_exp_with_approved_symbols <- approved_symbols[approved_symbols %in% rownames(exp)]
genes_not_in_exp_with_approved_symbols

# Initialize the new column with NA
genes$sc_Fetal_Region_Expressed <- NA
cutoff <- 0.01
for (i in seq_len(nrow(genes))) {
    gene <- genes$Gene[i]
    approved_symbol <- genes$Approved.symbol[i]
    
    if (!is.na(genes$Human_Period_Enrichment[i]) && genes$Human_Period_Enrichment[i] == "Prenatal") {
        expressed_regions <- c()
        
        if (gene %in% rownames(exp)) {
            expressed_regions <- colnames(exp)[exp[gene, ] > cutoff]
        } else if (approved_symbol %in% rownames(exp)) {
            expressed_regions <- colnames(exp)[exp[approved_symbol, ] > cutoff]
        }
        
        if (length(expressed_regions) > 0) {
            if (is.na(genes$sc_Fetal_Region_Expressed[i])) {
                genes$sc_Fetal_Region_Expressed[i] <- paste(expressed_regions, collapse = ", ")
            } else {
                genes$sc_Fetal_Region_Expressed[i] <- paste(genes$sc_Fetal_Region_Expressed[i], paste(expressed_regions, collapse = ", "), sep = ", ")
            }
        }
    }
}


# Print the number of genes with non-NA values in the last column
sum(!is.na(genes$sc_Fetal_Region_Expressed))
#630

#View(genes)


## Specificity

nonzero_values <- spec[spec > 0] # Ignore zeros if necessary
mean_spec <- mean(nonzero_values)
sd_spec <- sd(nonzero_values)
cutoff_spec <- mean_spec + sd_spec

# Initialize the new column with NA
genes$sc_Fetal_Region_Specificity <- NA

for (i in seq_len(nrow(genes))) {
    gene <- genes$Gene[i]
    approved_symbol <- genes$Approved.symbol[i]
    
    if (!is.na(genes$Human_Period_Enrichment[i]) && genes$Human_Period_Enrichment[i] == "Prenatal") {
        specific_regions <- c()
        
        if (gene %in% rownames(spec)) {
            specific_regions <- colnames(spec)[spec[gene, ] > cutoff_spec]
        } else if (approved_symbol %in% rownames(spec)) {
            specific_regions <- colnames(spec)[spec[approved_symbol, ] > cutoff_spec]
        }
        
        if (length(specific_regions) > 0) {
            if (is.na(genes$sc_Fetal_Region_Specificity[i])) {
                genes$sc_Fetal_Region_Specificity[i] <- paste(specific_regions, collapse = ", ")
            } else {
                genes$sc_Fetal_Region_Specificity[i] <- paste(genes$sc_Fetal_Region_Specificity[i], paste(specific_regions, collapse = ", "), sep = ", ")
            }
        }
    }
}

sum(!is.na(genes$sc_Fetal_Region_Specificity))
#24
#View(genes)


load("data/ctd_FetalCT2023.rda")

#head(ctd[[2]]$annot)

exp1 <- ctd[[2]]$mean_exp
spec1 <- ctd[[2]]$specificity

# remove cerebellum
exp1 <- exp1[, !grepl("cerebellum_", colnames(exp1))]
spec1 <- spec1[, !grepl("cerebellum_", colnames(spec1))]

#exp1 <- exp1[, colnames(exp1) != "Endothelial-Mural"]
#spec1 <- spec1[, colnames(spec1) != "Endothelial-Mural"]

# Define the cutoff (Mean + 1 SD)
nonzero_values1 <- exp1[exp1 > 0]
mean_expr1 <- mean(nonzero_values1)
std_expr1 <- sd(nonzero_values1)
cutoff1 <- mean_expr1 + std_expr1
print(paste("Expression cutoff:", cutoff1))

# going with 0.01 as cutoff since the other is too high. #point for discussion
cutoff1 <- 0.01

# Initialize the new column with NA
genes$sc_Fetal_CellType_Expressed <- NA
cutoff1 <- 0.01
for (i in seq_len(nrow(genes))) {
    gene <- genes$Gene[i]
    approved_symbol <- genes$Approved.symbol[i]
    
    if (!is.na(genes$Human_Period_Enrichment[i]) && genes$Human_Period_Enrichment[i] == "Prenatal") {
        expressed_regions <- c()
        
        if (gene %in% rownames(exp1)) {
            expressed_regions <- colnames(exp1)[exp1[gene, ] > cutoff1]
        } else if (approved_symbol %in% rownames(exp1)) {
            expressed_regions <- colnames(exp1)[exp1[approved_symbol, ] > cutoff1]
        }
        
        if (length(expressed_regions) > 0) {
            if (is.na(genes$sc_Fetal_CellType_Expressed[i])) {
                genes$sc_Fetal_CellType_Expressed[i] <- paste(expressed_regions, collapse = ", ")
            } else {
                genes$sc_Fetal_CellType_Expressed[i] <- paste(genes$sc_Fetal_CellType_Expressed[i], paste(expressed_regions, collapse = ", "), sep = ", ")
            }
        }
    }
}

sum(!is.na(genes$sc_Fetal_CellType_Expressed))
#652

#View(genes)


## Specificity

nonzero_values1 <- spec1[spec1 > 0] # Ignore zeros if necessary
mean_spec1 <- mean(nonzero_values1)
sd_spec1 <- sd(nonzero_values1)
cutoff_spec1 <- mean_spec1 + sd_spec1

# Initialize the new column with NA
genes$sc_Fetal_CellType_Specificity <- NA

for (i in seq_len(nrow(genes))) {
    gene <- genes$Gene[i]
    approved_symbol <- genes$Approved.symbol[i]
    
    if (!is.na(genes$Human_Period_Enrichment[i]) && genes$Human_Period_Enrichment[i] == "Prenatal") {
        specific_regions <- c()
        
        if (gene %in% rownames(spec1)) {
            specific_regions <- colnames(spec1)[spec1[gene, ] > cutoff_spec1]
        } else if (approved_symbol %in% rownames(spec1)) {
            specific_regions <- colnames(spec1)[spec[approved_symbol, ] > cutoff_spec1]
        }
        
        if (length(specific_regions) > 0) {
            if (is.na(genes$sc_Fetal_CellType_Specificity[i])) {
                genes$sc_Fetal_CellType_Specificity[i] <- paste(specific_regions, collapse = ", ")
            } else {
                genes$sc_Fetal_CellType_Specificity[i] <- paste(genes$sc_Fetal_CellType_Specificity[i], paste(specific_regions, collapse = ", "), sep = ", ")
            }
        }
    }
}
sum(!is.na(genes$sc_Fetal_CellType_Specificity))
#123
View(genes)

write.xlsx(genes, "results/TRIFF_risk_gene_SC_Fetal_info.xlsx")

### Summary of regions and cell types expressed

# Summarize the regions and cell types expressed
regions_summary <- table(unlist(strsplit(genes$sc_Fetal_Region_Expressed, ", ")))
celltypes_summary <- table(unlist(strsplit(genes$sc_Fetal_CellType_Expressed, ", ")))

regions_df <- as.data.frame(regions_summary)
colnames(regions_df) <- c("Region", "Count")
celltypes_df <- as.data.frame(celltypes_summary)
colnames(celltypes_df) <- c("CellType", "Count")

ggplot(regions_df, aes(x = reorder(Region, Count), y = Count, fill = Count)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_viridis_c() +
    labs(title = "Summary of Regions Expressed", x = "Region", y = "Count") +
    theme_minimal()
ggsave("results/Summary_regions_expressed.pdf")

# Plot the cell types expressed
ggplot(celltypes_df, aes(x = reorder(CellType, Count), y = Count, fill = Count)) +
    geom_bar(stat = "identity") +
     coord_flip() +
    viridis::scale_fill_viridis(option="C") +
    labs(title = "Summary of Cell Types Expressed", x = "Cell Type", y = "Count") +
    theme_minimal()
ggsave("results/Summary_Celltypes_expressed.pdf")



#Summarize the specificity regions and cell types
regions_spec_summary <- table(unlist(strsplit(genes$sc_Fetal_Region_Specificity, ", ")))
celltypes_spec_summary <- table(unlist(strsplit(genes$sc_Fetal_CellType_Specificity, ", ")))

regions_spec_df <- as.data.frame(regions_spec_summary)
colnames(regions_spec_df) <- c("Region", "Count")

celltypes_spec_df <- as.data.frame(celltypes_spec_summary)
colnames(celltypes_spec_df) <- c("CellType", "Count")

ggplot(regions_spec_df, aes(x = "", y = Count, fill = Region)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = pal3[c(4,6,2)]) +
    labs(title = "Specificity of Regions Expressed") +
    theme_minimal() + 
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(), panel.grid = element_blank(), axis.text.x = element_blank()) +
    geom_text(aes(label = scales::percent(Count / sum(Count))), position = position_stack(vjust = 0.5))
ggsave("results/Summary_Specificity_regions_expressed.pdf")

# Plot the cell types specificity in a pie chart
ggplot(celltypes_spec_df, aes(x = "", y = Count, fill = CellType)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    labs(title = "Specificity of Cell Types Expressed") +
    theme_minimal() +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(), panel.grid = element_blank(), axis.text.x = element_blank())
ggsave("results/Summary_Specificity_Celltypes_expressed.pdf")


pal <- c("#DFE1BC", "#90C6A7", "#3D646A", "#795862", 
             "#141313", "#E2E1E4", "#A49FB8", "#6A83A9") 
pal2 <- c("#A6A6A6", "#CDAA7D", "#EEE685", "#FFBBFF", "#D1EEEE", "#6495ED", "#6959CD")
pal3<- c("#0E1118", "#DC383A", "#A90C1B", "#FED568", "#C48233", "#84C4EF", "#1D4564")
