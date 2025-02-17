#### TRIFF Project #####

## Code to explore expression of genes in the brain (adapted from mmkim1210) ##

# First round of exploration will be done using the BrainSpan data - bulk RNA-seq data from the developing human brain

source("R/utils.R")
rm(list = ls())
options(stringsAsFactors = FALSE)


suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(cowplot)))
library(readxl)

system("curl http://development.psychencode.org/files/raw_data/mRNA-seq_Sample%20metadata.xlsx --output metadata_subject.xlsx")
system("curl http://development.psychencode.org/files/processed_data/RNA-seq/mRNA-seq_QC.xlsx --output metadata_sample.xlsx")
system("curl http://development.psychencode.org/files/processed_data/RNA-seq/mRNA-seq_hg38.gencode21.wholeGene.geneComposite.STAR.nochrM.gene.RPKM.normalized.CQNCombat.txt --output data_normalized.txt")


datMeta_subject <- read_excel("metadata_subject.xlsx", skip = 3)
datMeta_subject <- slice(datMeta_subject, -1)
#datMeta_subject %>% print(n = 5, width = Inf)

datMeta_sample <- read_excel("metadata_sample.xlsx", skip = 3)
#datMeta_sample %>% print(n = 5, width = Inf)
datMeta_sample <- slice(datMeta_sample, -1) %>% 
  select(2:3) %>% 
  unite(ID, c(Braincode, Regioncode), sep = ".", remove = FALSE)

datMeta <- left_join(datMeta_sample, datMeta_subject, by = "Braincode")  # 607 samples




## Subset data to include selected brain regions (mentioned in the project plan)
# Check which regiona are included in the data
table(datMeta$Regioncode)

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


expr <- as.data.frame(fread("data_normalized.txt"))
rownames(expr) <- expr$Geneid
expr <- select(expr, -1)  # 60,155 x 607

annot <- data.frame(Geneid = rownames(expr))
annot[, c("ENSG", "gene")] <- str_split_fixed(annot$Geneid,"[|]", 2)

ind <- match(colnames(expr), datMeta$ID)
#sum(is.na(ind))
datMeta <- datMeta[ind, ]


### Read in disease-risk genes ###

disease_list1<- readRDS("/Volumes/projects/C3_Sellgren_lab/Lab Members/Susmita/Internal data/C1Q/data/Geneset_with_disease_riskgenes.rds")
# Lets start with the latest SCZ gwas from 2022
trubetskoy <- readRDS("/Volumes/projects/C3_Sellgren_lab/Lab Members/Susmita/Internal data/C1Q/Bipolar/Trubetskoy_scz_genes.rds")
hic <- disease_list1$SCZ_HiC$`High-confident Hi-C defined Schizophrenia risk genes`
schema <- c("SET1DA", "CUL1", "XPO7", "TRIO", "CACNA1G", "SP4", "GRIA3", "GRIN2A", "HERC1", "RB1CC1")

### According to our discussions, for the rest we basically go with every CNV with odds ratio > 5.0 from Owen et al., 2023. 
#load all 127 cnvs
all_cnv<- readRDS("/Volumes/projects/C3_Sellgren_lab/Lab Members/Susmita/Internal data/CNV/2023/Genes_with_CNV.rds")
# this includes 22q11.2del, PWS/AS dup, 3q29 del, NRXN1 del, 16p11.2 dup and 1q21.1. del

keep <- c("22q11.21_DUP/DEL","7q11.23_DUP/DEL", "3q29_DUP/DEL", "16p11.2_DUP/DEL-A","16p11.2_DUP/DEL-B" ,"1q21.1-q21.2_DUP/DEL")

# Dont forget NRXN1

cnvs <- list("22q11.2del" = c(), "PWS/ASdup" = c(), "3q29del"= c(), "NRXN1del"= c(), "16p11.2dup" =c(), "1q21.1del" = c())


#### Begin Expression Exploration here ####

# Create a data frame to store gene expression information
gene_expression_info <- data.frame(
    Dataset = "BrainSpan",
    Source = "GWAS (Trubetskoy et al., 2022)",
    Gene = trubetskoy$all_scz_trubetskoy,
    BrainExpressed = NA,
    BrainExpressionLevel = NA,
    AvgExpressionPrenatal = NA,
    AvgExpressionPostnatal = NA,
    Significance_Pval = NA,
    Significance_Bestimate = NA,
    Human_Period_Enrichment = NA)



# Populate the data frame with average expression values
for (i in 1:nrow(gene_expression_info)) {
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
        sig <- calculate_significance(gene)
        gene_expression_info$Significance_Pval[i] <- sig[[1]]
        gene_expression_info$Significance_Bestimate[i] <- sig[[2]]
        if (gene_expression_info$Significance_Pval[i] < 0.05) {
            gene_expression_info$Human_Period_Enrichment[i] <- ifelse(gene_expression_info$Significance_Bestimate[i] < 0, "Prenatal", "Postnatal")
        } else {
            gene_expression_info$Human_Period_Enrichment[i] <- NA
        }
        }
        # save the plot
        ggsave(paste0("./results/Individual Gene Exp Plots/",gene, "_BrainSpan_Human.pdf"), plot = plot_expr(gene), width = 8, height = 7, dpi = 600)
        
}

# Print the gene expression information
View(gene_expression_info)

openxlsx::write.xlsx(gene_expression_info, "./results/TRIFF_risk_gene_info.xlsx")


