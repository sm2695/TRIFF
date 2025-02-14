# Function to check if a gene is expressed in the brain
is_gene_expressed <- function(gene) {
    ind <- which(annot$gene == gene)
    if (length(ind) == 0) {
        return(FALSE)
    }
    if(length(ind) > 1) {
        gene_data <- as.numeric(rowMeans(expr[ind, ]))
    }else {
       gene_data <- as.numeric(expr[ind, ])
    }
    
    if (mean(gene_data) >= 0.5){
        return(TRUE)
    } else {
        return(FALSE)
    }
}



## Plot gene epression ###
plot_expr <- function(gene) {
    ind <- which(annot$gene == gene)
    if (length(ind) > 1) {
        ind <- ind[1]
    }
    data.frame(RPKM = as.numeric(expr[ind, ]), datMeta) %>% 
      filter(!is.na(Region)) %>% 
      ggplot(aes(x = Days, y = log2(.1 + RPKM), fill = Region, color = Region)) + 
        geom_point(aes(shape = Period)) + 
        scale_x_log10() +
        geom_smooth(method = "loess", formula = y ~ x, alpha = .1, span = 1, se = FALSE) + 
        labs(x = "Days post conception") + 
        ggtitle(annot$gene[ind]) +
        theme_bw() +
        theme(plot.title = element_text(hjust = .5), panel.grid = element_blank()) +
        geom_vline(xintercept = 40 * 7, lty = 1) + 
        scale_color_manual(values = c("Frontal Cortex" = "#E41A1C", "Temporal Cortex" = "#377EB8", "Parietal Cortex" = "#4DAF4A", "Visual Cortex" = "#984EA3", "Thalamus" = "#FF7F00", "Striatum" = "#FFFF33", "Hippocampus" = "#A65628", "Amygdala" = "#F781BF")) +
        scale_fill_manual(values = c("Frontal Cortex" = "#E41A1C", "Temporal Cortex" = "#377EB8", "Parietal Cortex" = "#4DAF4A", "Visual Cortex" = "#984EA3", "Thalamus" = "#FF7F00", "Striatum" = "#FFFF33", "Hippocampus" = "#A65628", "Amygdala" = "#F781BF"))
}



# Calculate p-value for prenatal expression of the gene
calculate_significance <- function(gene) {
    ind <- which(annot$gene == gene)
   
    if(length(ind) > 1) {
        # take the first hit for the gene. The second one comes up most likely due to a novel protein being mapped to it. Note that expression levels are indeed different though. Troubleshoot this. 
        gene_data <- data.frame(RPKM = as.numeric(expr[ind[1], ]), datMeta) %>% 
            filter(!is.na(Region))

        #gene_data <- data.frame(RPKM = as.numeric(rowMeans(expr[ind, ])), datMeta) %>% 
        #filter(!is.na(Region))
    } else {
            gene_data <- data.frame(RPKM = as.numeric(expr[ind, ]), datMeta) %>% 
            filter(!is.na(Region))
    }
    gene_data <- data.frame(RPKM = as.numeric(expr[ind[1], ]), datMeta) %>% 
            filter(!is.na(Region))
    
    # Perform a lm comparing expression in different periods, including region as a covariate
    gene_data$Region <- factor(gene_data$Region)
    gene_data$Sex <- factor(gene_data$Sex)
    gene_data$Ethnicity <- factor(gene_data$Ethnicity)
    #gene_data$ID <- factor(gene_data$ID)
    #gene_data$Sequencing.Site <- factor(gene_data$Sequencing.Site)
    gene_data$Region <- factor(gene_data$Region)
    
    #Get p-value and B-estimates for the models
    model <- lm(log2(.1 + RPKM) ~ Period + Region + Sex + Ethnicity, data = gene_data)
    lm_test_result_p <- summary(model)$coefficients["PeriodPostnatal", "Pr(>|t|)"]
    lm_test_result_b <- summary(model)$coefficients["PeriodPostnatal", "Estimate"]
    #t_test_result <- t.test(log2(.1 + gene_data$RPKM) ~ gene_data$Period)
    return(list(lm_test_result_p, lm_test_result_b))
}
