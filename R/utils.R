#' Check if a gene is expressed in the brain
#'
#' This function checks whether a given gene is expressed in the brain based on expression threshold.
#'
#' @param gene A character string specifying the gene name.
#' @return TRUE if the gene is expressed, FALSE otherwise.
#' @examples
#' is_gene_expressed("GENE1")
is_gene_expressed <- function(gene) {
  ind <- which(annot$gene == gene)
  
  if (length(ind) == 0) {
    return(FALSE)
  }
  
  if (length(ind) > 1) {
    gene_data <- as.numeric(rowMeans(expr[ind, ]))
  } else {
    gene_data <- as.numeric(expr[ind, ])
  }
  
  return(mean(gene_data) >= 0.5)
}

#' Plot Gene Expression
#'
#' This function generates a plot of gene expression over time, with different regions highlighted.
#'
#' @param gene A character string specifying the gene name.
#' @return A ggplot2 object visualizing gene expression.
#' @examples
#' plot_expr("GENE1")
plot_expr <- function(gene) {
  ind <- which(annot$gene == gene)
  if (length(ind) > 1) {
    ind <- ind[1]
  }
  
  gene_data <- data.frame(RPKM = as.numeric(expr[ind, ]), datMeta) %>% 
    filter(!is.na(Region))
  
  ggplot(gene_data, aes(x = Days, y = log2(.1 + RPKM), fill = Region, color = Region)) + 
    geom_point(aes(shape = Period)) + 
    scale_x_log10() +
    geom_smooth(method = "loess", formula = y ~ x, alpha = .1, span = 1, se = FALSE) + 
    labs(x = "Days post conception", title = annot$gene[ind]) + 
    theme_bw() +
    theme(plot.title = element_text(hjust = .5), panel.grid = element_blank()) +
    geom_vline(xintercept = 40 * 7, lty = 1) + 
    scale_color_manual(values = c(
      "Frontal Cortex" = "#E41A1C", "Temporal Cortex" = "#377EB8", "Parietal Cortex" = "#4DAF4A", 
      "Visual Cortex" = "#984EA3", "Thalamus" = "#FF7F00", "Striatum" = "#FFFF33", 
      "Hippocampus" = "#A65628", "Amygdala" = "#F781BF"
    )) +
    scale_fill_manual(values = c(
      "Frontal Cortex" = "#E41A1C", "Temporal Cortex" = "#377EB8", "Parietal Cortex" = "#4DAF4A", 
      "Visual Cortex" = "#984EA3", "Thalamus" = "#FF7F00", "Striatum" = "#FFFF33", 
      "Hippocampus" = "#A65628", "Amygdala" = "#F781BF"
    ))
}

#' Calculate Significance of Prenatal Expression
#'
#' This function calculates the p-value and effect size for prenatal expression of a given gene.
#'
#' @param gene A character string specifying the gene name.
#' @return A list containing the p-value and effect size.
#' @examples
#' calculate_significance("GENE1")
calculate_significance <- function(gene) {
  ind <- which(annot$gene == gene)
  
  if (length(ind) > 1) {
    ind <- ind[1]  # Take the first match to avoid duplicate mappings
  }
  
  gene_data <- data.frame(RPKM = as.numeric(expr[ind, ]), datMeta) %>% 
    filter(!is.na(Region))
  
  # Convert categorical variables to factors
  gene_data <- gene_data %>% mutate(
    Region = factor(Region),
    Sex = factor(Sex),
    Ethnicity = factor(Ethnicity),
    ID = factor(ID),
    Sequencing.Site = factor(Sequencing.Site)
  )
  
  # Perform linear model analysis
  model <- lm(log2(.1 + RPKM) ~ Period + Region + Sex + Ethnicity, data = gene_data)
  lm_test_result_p <- summary(model)$coefficients["PeriodPostnatal", "Pr(>|t|)"]
  lm_test_result_b <- summary(model)$coefficients["PeriodPostnatal", "Estimate"]
  
  return(list(p_value = lm_test_result_p, effect_size = lm_test_result_b))
}

#' Retrieve gene biotype information from Ensembl
#'
#' This function fetches gene biotype information from Ensembl using the biomaRt package.
#' It retrieves data for human genes and their corresponding mouse homologs.
#'
#' @param gene_symbols A character vector of HGNC gene symbols.
#' @return A data frame containing gene biotype information, descriptions, and mouse homologs.
#' @examples
#' get_gene_info(c("TP53", "BRCA1"))
#' 
#' @import biomaRt
#' @export
get_gene_info <- function(gene_symbols) {
  # Load required package
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop("Package 'biomaRt' is required but not installed. Install it using install.packages('biomaRt').")
  }
  
  # STEP 1: Retrieve Human Gene Data
  human_attributes <- c("hgnc_symbol", "external_gene_name", "gene_biotype", "description", "ensembl_gene_id")  
  human_data <- biomaRt::getBM(
    attributes = human_attributes, 
    filters = "hgnc_symbol", 
    values = gene_symbols, 
    mart = ensembl_human
  )
  
  # STEP 2: Retrieve Mouse Homologs using the Mouse Dataset
  mouse_attributes <- c("external_gene_name", "hsapiens_homolog_associated_gene_name")
  mouse_homologs <- biomaRt::getBM(
    attributes = mouse_attributes, 
    mart = ensembl_mouse
  )
  
  # Rename columns for merging
  colnames(mouse_homologs) <- c("mouse_gene_name", "hgnc_symbol")
  
  # STEP 3: Merge Human Data with Mouse Homologs
  final_data <- merge(human_data, mouse_homologs, by = "hgnc_symbol", all.x = TRUE)
  
  return(final_data)
}

#' Query Exphewas API for a single gene
#'
#' This function queries the Exphewas API to retrieve data for a given gene.
#'
#' @param gene_name A character string specifying the gene name.
#' @return A list containing parsed JSON data from the API response, or NULL if the request fails.
#' @examples
#' query_exphewas("TP53")
#'
#' @import httr
#' @import jsonlite
#' @export
query_exphewas <- function(gene_name) {
  # Construct API URL
  url <- paste0("https://exphewas.statgen.org/v1/api/gene/name/", gene_name)
  
  # Send GET request to the API
  response <- httr::GET(url)
  
  # Check if the request was successful
  if (httr::status_code(response) == 200) {
    content_data <- httr::content(response, "text", encoding = "UTF-8")
    parsed_data <- jsonlite::fromJSON(content_data, flatten = TRUE)
    return(parsed_data)
  } else {
    warning(paste("Failed to fetch data for gene:", gene_name))
    return(NULL)
  }
}


#' Function to get top GO terms for each gene
#'
#' This function extracts the GO:BP, GO:CC and GO:MF terms for a given gene.
#'
#' @param gene_name A character string specifying the gene name.
#' @return A list containing top GO terms, or NULL if the request fails.
#' @examples
#' get_top_go_terms("TP53")
#'
#' @import clusterProfiler
#' @import org.Hs.eg.db
#' @export
get_top_go_terms <- function(gene_name) {
    tryCatch({
        # Convert gene name to Entrez ID
        entrez_id <- bitr(gene_name, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
        if (length(entrez_id) == 0) return(list(BP = NA, CC = NA))
        
        # Perform GO enrichment analysis for Biological Process (BP)
        go_bp <- enrichGO(entrez_id, OrgDb = org.Hs.eg.db, ont = "BP", readable = TRUE)
        top_bp <- if (!is.null(go_bp)) paste(head(paste(go_bp@result$ID, go_bp@result$Description, sep = ": "), 10), collapse = "; ") else NA
        
        # Perform GO enrichment analysis for Cellular Component (CC)
        go_cc <- enrichGO(entrez_id, OrgDb = org.Hs.eg.db, ont = "CC", readable = TRUE)
        top_cc <- if (!is.null(go_cc)) paste(head(paste(go_cc@result$ID, go_cc@result$Description, sep = ": "), 3), collapse = "; ") else NA

        # Perform GO enrichment analysis for Molecular Function (MF)
        go_mf <- enrichGO(entrez_id, OrgDb = org.Hs.eg.db, ont = "MF", readable = TRUE)
        top_mf <- if (!is.null(go_mf)) paste(head(paste(go_mf@result$ID, go_mf@result$Description, sep = ": "), 3), collapse = "; ") else NA
        
        return(list(BP = top_bp, CC = top_cc, MF=top_mf))
    }, error = function(e) {
        return(list(BP = NA, CC = NA, MF = NA))
    })
}

