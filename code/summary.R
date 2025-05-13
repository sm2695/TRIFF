### Summaries of the TRIFF gene sets ####

library(openxlsx)
library(dplyr)

df <- read.xlsx("TRIFF_genes_info_Ontologies_Hs&Ms.xlsx")
head(df)

table(df$BrainExpressed)

table(df$BrainExpressed, df$Human_Period_Enrichment)

table(df$Human_Period_Enrichment, df$MouseBrainExpressed)

library(VennDiagram)

brain_expressed <- which(df$BrainExpressed == TRUE)
human_enriched <- which(!is.na(df$Human_Period_Enrichment) & df$Human_Period_Enrichment == "Prenatal")
mouse_expressed <- which(df$MouseBrainExpressed == TRUE)
human_enriched_mouse_expressed <- which(!is.na(df$Human_Period_Enrichment) & df$Human_Period_Enrichment == "Prenatal" & df$MouseBrainExpressed == TRUE)





library(ggplot2)

# A. How many genes are expressed in the brain vs not expressed
brain_expr_counts <- table(df$BrainExpressed)

pie(brain_expr_counts,
    labels = paste0(names(brain_expr_counts), ": ", brain_expr_counts),
    col = c("darkred", "darkblue"),
    main = "A. Genes Expressed vs Not Expressed in Brain")


# B. Of expressed genes, how many are enriched in Prenatal
brain_expr_df <- subset(df, BrainExpressed == TRUE)
prenatal_counts <- table(brain_expr_df$Human_Period_Enrichment == "Prenatal")

# Bar plot B
pie(prenatal_counts,
    labels = paste0(names(prenatal_counts), ": ", prenatal_counts),
    col = c("purple", "turquoise"),
    main = "B. Prenatal Enrichment in Brain-Expressed Genes")

brain_prenatal_df <- subset(brain_expr_df, Human_Period_Enrichment == "Prenatal")
mouse_counts <- table(brain_prenatal_df$MouseBrainExpressed)
names(mouse_counts) <- c("Not Mouse Expressed", "Mouse Expressed")

pie(mouse_counts,
    labels = paste0(names(mouse_counts), ": ", mouse_counts),
    col = c("orange", "blue"),
    main = "C. Mouse Expression in Prenatal+Brain Genes")



df1 <- df

library(dplyr)
library(tidyr)
df$Source <- gsub("^CNV[^,]*", "CNV", df$Source)

# Split sources by commas
df_sources <- df %>%
  mutate(SourceList = strsplit(as.character(Source), ",")) %>%
  unnest(SourceList) %>%
  mutate(Value = TRUE) %>%
  pivot_wider(names_from = SourceList, values_from = Value, values_fill = FALSE)

# View column names of distinct sources (for debugging)
colnames(df_sources)
source_only_df <- df_sources %>% select(where(is.logical))

library(UpSetR)
upset(source_only_df, 
      nsets = ncol(source_only_df),     # Adjust number of sets based on columns
      nintersects = 30,                 # Limit the number of intersections shown
      order.by = "freq")    

library(EWCE)

load("/Users/susmita.malwade/ctd_Oligo.rda")


org_exp <- ctd[[1]]$mean_exp



genes <- df1
genes$sc_Organoid_CellType_Expressed <- NA
cutoff1 <- 0.01
for (i in seq_len(nrow(genes))) {
    gene <- genes$Gene[i]
    approved_symbol <- genes$Approved.symbol[i]
    
    if (!is.na(genes$Human_Period_Enrichment[i]) && genes$Human_Period_Enrichment[i] == "Prenatal") {
        expressed_regions <- c()
        
        if (gene %in% rownames(org_exp)) {
            expressed_regions <- colnames(org_exp)[org_exp[gene, ] > cutoff1]
        } else if (approved_symbol %in% rownames(org_exp)) {
            expressed_regions <- colnames(exp)[org_exp[approved_symbol, ] > cutoff1]
        }
        
        if (length(expressed_regions) > 0) {
            if (is.na(genes$sc_Organoid_CellType_Expressed[i])) {
                genes$sc_Organoid_CellType_Expressed[i] <- paste(expressed_regions, collapse = ", ")
            } else {
                genes$sc_Organoid_CellType_Expressed[i] <- paste(genes$sc_Organoid_CellType_Expressed[i], paste(expressed_regions, collapse = ", "), sep = ", ")
            }
        }
    }
}

sum(!is.na(genes$sc_Organoid_CellType_Expressed))

genes$sc_Organoid_Expressed <- ifelse(!is.na(genes$sc_Organoid_CellType_Expressed), TRUE, FALSE)

write.xlsx(genes, file = "TRIFF_genes_info_Ontologies_Hs&Ms&Org.xlsx", sheetName = "Sheet1", append = TRUE, row.names = FALSE)


brain_expr_df <- subset(genes, BrainExpressed == TRUE)
brain_prenatal_df <- subset(brain_expr_df, Human_Period_Enrichment == "Prenatal")
brain_prenatal_mouse <- subset(brain_prenatal_df, MouseBrainExpressed == TRUE)
#brain_prenatal_mouse_sc <- subset(brain_prenatal_df, sc_Organoid_Expressed == TRUE)
organoid_counts <- table(brain_prenatal_mouse$sc_Organoid_Expressed)
names(organoid_counts) <- c("Not Organoid Expressed", "Organoid Expressed")
#table(genes$sc_Organoid_Expressed, genes$MouseBrainExpressed)

png("Organoid_Expression.png", width = 800, height = 600)
pie(organoid_counts,
    labels = paste0(names(organoid_counts), ": ", organoid_counts),
    col = c("orange", "forestgreen"),
    main = "D. Organoid Expression in Prenatal Brain + Mouse Genes")
dev.off()


library(stringr)

# Step 1: Split entries and unlist
celltypes <- str_split(main$sc_Fetal_CellType_Expressed, ",\\s*") %>% unlist()

# Step 2: Trim whitespace
celltypes <- str_trim(celltypes)

# Step 3: Create a data frame of counts
celltype_df <- as.data.frame(table(celltypes))
colnames(celltype_df) <- c("CellType", "Count")

# Step 4: Plot barplot
ggplot(celltype_df, aes(x = reorder(CellType, -Count), y = Count)) +
  geom_bar(stat = "identity", fill = "darkgreen") +
  labs(title = "Expression Count per Fetal Cell Type", x = "Cell Type", y = "Gene Count") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.title.x = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5)
  )
ggsave("Fetal_CellType_Expression_Counts.pdf", width = 10, height = 6, dpi = 300)


# Step 1: Split entries and unlist
specificity <- str_split(main$sc_Fetal_CellType_Specificity, ",\\s*") %>% unlist()

# Step 2: Trim whitespacet
specificity <- str_trim(specificity)

# Step 3: Create a data frame of counts
specificity_df <- as.data.frame(table(specificity))
colnames(specificity_df) <- c("CellType", "Count")

# Step 4: Plot barplot
ggplot(specificity_df, aes(x = reorder(CellType, -Count), y = Count)) +
  geom_bar(stat = "identity", fill = "purple") +
  labs(title = "Specific Expression per Fetal Cell Type", x = "Cell Type", y = "Gene Count") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.title.x = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5)
  )
ggsave("Fetal_CellType_Specificity_Counts.pdf", width = 10, height = 6, dpi = 300)

