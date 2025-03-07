#### TRIFF Project #####

## Code to explore expression of genes in the brain (adapted from mmkim1210) ##

# First round of exploration will be done using the BrainSpan data - bulk RNA-seq data from the developing human brain

source("R/utils.R")
rm(list = ls())
options(stringsAsFactors = FALSE)


if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
if (!requireNamespace("cowplot", quietly = TRUE)) install.packages("cowplot")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")
if (!requireNamespace("parallel", quietly = TRUE)) install.packages("parallel")
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(cowplot)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(readxl)))
suppressWarnings(suppressMessages(library(openxlsx)))
suppressWarnings(suppressMessages(library(parallel)))


# Create directory for BrainSpan if it doesn't exist
if (!dir.exists("data/BrainSpan")) {
    dir.create("data/BrainSpan", recursive = TRUE)
}

system("curl http://development.psychencode.org/files/raw_data/mRNA-seq_Sample%20metadata.xlsx --output data/BrainSpan/metadata_subject.xlsx")
system("curl http://development.psychencode.org/files/processed_data/RNA-seq/mRNA-seq_QC.xlsx --output data/BrainSpan/metadata_sample.xlsx")
system("curl http://development.psychencode.org/files/processed_data/RNA-seq/mRNA-seq_hg38.gencode21.wholeGene.geneComposite.STAR.nochrM.gene.RPKM.normalized.CQNCombat.txt --output data/BrainSpan/data_normalized.txt")


datMeta_subject <- read_excel("data/BrainSpan/metadata_subject.xlsx", skip = 3)
datMeta_subject <- slice(datMeta_subject, -1)
#datMeta_subject %>% print(n = 5, width = Inf)

datMeta_sample <- read_excel("data/BrainSpan/metadata_sample.xlsx", skip = 3)
#datMeta_sample %>% print(n = 5, width = Inf)
datMeta_sample <- slice(datMeta_sample, -1) %>% 
  select(2:3) %>% 
  unite(ID, c(Braincode, Regioncode), sep = ".", remove = FALSE)

datMeta <- left_join(datMeta_sample, datMeta_subject, by = "Braincode")  # 607 samples

## Subset data to include selected brain regions (mentioned in the project plan)
# Check which regiona are included in the data
#table(datMeta$Regioncode)

#   A1C    AMY    CBC    CGE    DFC    DTH    HIP    IPC    ITC    LGE    M1C 
#    36     37     33      2     40      2     39     39     37      2     33 
#M1CS1C     MD    MFC    MGE    MSC     OC    OFC     PC    S1C    STC    STR 
#     2     33     39      2      4      2     37      2     34     37     34 
#    TC    URL    V1C    VFC 
#     1      2     38     40 


datMeta <- datMeta %>% 
  mutate(Region = ifelse(Regioncode %in% c("DFC", "MFC", "PC", "VFC", "M1C", "OFC", "FC"), "Frontal Cortex",
                        ifelse(Regioncode %in% c("A1C", "ITC", "STC", "TC"), "Temporal Cortex",
                               ifelse(Regioncode %in% c("IPC", "S1C", "PC"), "Parietal Cortex",
                                      ifelse(Regioncode %in% c("OC", "V1C"), "Visual Cortex",
                                             ifelse(Regioncode %in% c("CB", "CBC", "URL"), "Cerebellum",
                                                    ifelse(Regioncode %in% c("DTH", "MD", "DIE"), "Thalamus",
                                                           ifelse(Regioncode %in% c("CGE", "LGE" ,"MGE", "STR"), "Striatum",
                                                                  ifelse(Regioncode == "HIP", "Hippocampus", 
                                                                         ifelse(Regioncode == "AMY","Amygdala", NA))))))))))

# remove cerebellum
datMeta <- filter(datMeta, Region != "Cerebellum")


datMeta <- mutate(datMeta, Period = ifelse(Days < 40 * 7, "Prenatal", "Postnatal"))
datMeta$Period <- factor(datMeta$Period, levels = c("Prenatal", "Postnatal"))

#datMeta %>% print(n = nrow(.), width = Inf)


expr <- as.data.frame(fread("data/BrainSpan/data_normalized.txt"))
rownames(expr) <- expr$Geneid
expr <- select(expr, -1)  # 60,155 x 607

annot <- data.frame(Geneid = rownames(expr))
annot[, c("ENSG", "gene")] <- str_split_fixed(annot$Geneid,"[|]", 2)

ind <- match(colnames(expr), datMeta$ID)
#sum(is.na(ind))
datMeta <- datMeta[ind, ]


### Read in disease-risk genes ###
gene_main_df <- read.xlsx("results/TRIFF_SCZ_curated_list_of_genes_updated.xlsx")

#### Begin Expression Exploration here ####

#Create a data frame to store gene expression information 
gene_expression_info <- gene_main_df %>% mutate(
    Dataset = "BrainSpan",
    BrainExpressed = NA,
    BrainExpressionLevel = NA,
    AvgExpressionPrenatal = NA,
    AvgExpressionPostnatal = NA,
    Significance_Pval = NA,
    Significance_Bestimate = NA,
    Human_Period_Enrichment = NA)

#Create directory for individual gene expression plots if it doesn't exist
if (!dir.exists("results/Individual Gene Exp Plots")) {
        dir.create("results/Individual Gene Exp Plots", recursive = TRUE)
    }


n_cores <- 6 # update if running on the server

#Populate the data frame with average expression values
for (i in seq_len(nrow(gene_expression_info))) {
    gene <- gene_expression_info$Gene[i]
    gene_expression_info$BrainExpressed[i] <- is_gene_expressed(gene)
    ind <- which(annot$gene == gene)
    gene_data <- as.numeric(unlist(expr[ind, ]))
    gene_expression_info$BrainExpressionLevel[i] <- mean(gene_data)
        
        if (gene_expression_info$BrainExpressed[i]) {
        
        # Calculate average expression for Prenatal and Postnatal periods
        gene_data_df <- data.frame(RPKM = gene_data, datMeta)
        prenatal_data <- filter(gene_data_df, Period == "Prenatal")
        postnatal_data <- filter(gene_data_df, Period == "Postnatal")
        gene_expression_info$AvgExpressionPrenatal[i] <- mean(prenatal_data$RPKM, na.rm = TRUE)
        gene_expression_info$AvgExpressionPostnatal[i] <- mean(postnatal_data$RPKM, na.rm = TRUE)

        # Calculate p-value for prenatal/postnatal enrichment using linear models
        sig <- mclapply(gene, calculate_significance, mc.cores = n_cores)
        gene_expression_info$Significance_Pval[i] <- sig[[1]]$p_value
        gene_expression_info$Significance_Bestimate[i] <- sig[[1]]$effect_size
        if (gene_expression_info$Significance_Pval[i] < 0.05) {
            gene_expression_info$Human_Period_Enrichment[i] <- ifelse(gene_expression_info$Significance_Bestimate[i] < 0, "Prenatal", "Postnatal")
        } else {
            gene_expression_info$Human_Period_Enrichment[i] <- NA
        }
        }
        # save the plot
        ggsave(paste0("results/Individual Gene Exp Plots/",gene, "_BrainSpan_Human.pdf"), plot = plot_expr(gene), width = 8, height = 7, dpi = 600)
        
}


#View(gene_expression_info)


# Identify genes in gene_main_df that are not in annot
missing_genes <- gene_main_df$Gene[!gene_main_df$Approved.symbol %in% annot$gene]
missing_df<- gene_main_df[which(gene_main_df$Gene %in% missing_genes),]
View(gene_main_df[which(gene_main_df$Gene %in% missing_genes),])

# check if the updated Approved symbols of these genes are present in annot
missing_approved <- annot$gene[annot$gene %in% missing_df$Approved.symbol]
missing_approved#
#"ELN"      "CHD5"     "CHD6"     "MPHOSPH6" "GALNTL6"

## Fill in missing gene information for these 3 genes

for (gene in missing_approved) { 
    row_index <- which(gene_expression_info$Approved.symbol == gene)
    if (length(row_index) == 1) {
        gene_expression_info$BrainExpressed[row_index] <- is_gene_expressed(gene)
        gene_data <- as.numeric(unlist(expr[which(annot$gene == gene), ]))
        gene_expression_info$BrainExpressionLevel[row_index] <- mean(gene_data)
        gene_data_df <- data.frame(RPKM = gene_data, datMeta)
        gene_expression_info$AvgExpressionPrenatal[row_index] <- mean(filter(gene_data_df, Period == "Prenatal")$RPKM, na.rm = TRUE)
        gene_expression_info$AvgExpressionPostnatal[row_index] <- mean(filter(gene_data_df, Period == "Postnatal")$RPKM, na.rm = TRUE)
        sig <- calculate_significance(gene)
        gene_expression_info$Significance_Pval[row_index] <- sig[[1]]
        gene_expression_info$Significance_Bestimate[row_index] <- sig[[2]]
        if (gene_expression_info$Significance_Pval[row_index] < 0.05) {
            gene_expression_info$Human_Period_Enrichment[row_index] <- ifelse(gene_expression_info$Significance_Bestimate[row_index] < 0, "Prenatal", "Postnatal")
        } else {
            gene_expression_info$Human_Period_Enrichment[row_index] <- NA
        }
        ggsave(paste0("results/Individual Gene Exp Plots/", gene, "_BrainSpan_Human.pdf"), plot = plot_expr(gene), width = 8, height = 7, dpi = 600)
    }
}

# Check if any columns are a list in gene_expression_info otherwise it doesn't save it properly
list_columns <- sapply(gene_expression_info, is.list)
if (any(list_columns)) {
    print("The following columns are lists:")
    print(names(gene_expression_info)[list_columns])
} else {
    print("No columns are lists.")
}

gene_expression_info$Significance_Pval <- as.numeric(gene_expression_info$Significance_Pval)

write.xlsx(gene_expression_info, "results/TRIFF_risk_gene_BrainSpan_info.xlsx")

#gene_expression_info <- read.xlsx("results/TRIFF_risk_gene_BrainSpan_info.xlsx")


# Summarize and Plot the percentage of genes that are brain expressed

#get some colors
pal <- c("#DFE1BC", "#90C6A7", "#3D646A", "#795862", 
             "#141313", "#E2E1E4", "#A49FB8", "#6A83A9") 
pal2 <- c("#A6A6A6", "#CDAA7D", "#EEE685", "#FFBBFF", "#D1EEEE", "#6495ED", "#6959CD")
pal3<- c("#0E1118", "#DC383A", "#A90C1B", "#FED568", "#C48233", "#84C4EF", "#1D4564")


brain_expressed_summary <- gene_expression_info %>%
    mutate(BrainExpressed = ifelse(is.na(BrainExpressed), FALSE, BrainExpressed)) %>%
    group_by(BrainExpressed) %>%
    summarise(Count = n()) %>%
    mutate(Percentage = Count / sum(Count) * 100)

#Plot pie chart for brain expressed genes
brain_expressed_pie <- ggplot(brain_expressed_summary, aes(x = "", y = Percentage, fill = BrainExpressed)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    theme_void() +
    scale_fill_manual(values = c("#A6A6A6", "#6A83A9")) +
    geom_text(aes(label = paste0(round(Percentage, 1), "%")), position = position_stack(vjust = 0.5)) +
    labs(title = "Percentage of Genes that are Brain Expressed")

#Summarize the human period enrichment for brain expressed genes
enrichment_summary <- gene_expression_info %>%
    filter(BrainExpressed == TRUE) %>%
    group_by(Human_Period_Enrichment) %>%
    summarise(Count = n()) %>%
    mutate(Percentage = Count / sum(Count) * 100)

#Plot pie chart for human period enrichment
enrichment_pie <- ggplot(enrichment_summary, aes(x = "", y = Percentage, fill = Human_Period_Enrichment)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    theme_void() +
    scale_fill_manual(values = pal3[c( 3 , 5, 7)]) +
    geom_text(aes(label = paste0(round(Percentage, 1), "%")), position = position_stack(vjust = 0.5)) +
    labs(title = "Human Period Enrichment for Brain Expressed Genes")

cowplot::plot_grid(brain_expressed_pie, enrichment_pie, nrow = 1)
ggsave("results/Summary_BrainSpan_results.pdf", width = 10, height = 5, dpi = 600)



#Remove temporary files from results directory
files_to_remove <- c("results/TRIFF_genes_unupdated_tmp.xlsx", "results/hgnc-symbol-check.csv", "results/hgnc-symbol-check.xlsx" )
file.remove(files_to_remove)
