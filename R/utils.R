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
