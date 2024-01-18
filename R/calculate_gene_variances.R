#' Calculate Gene Variances
#'
#' This function calculates the variances associated with different sources (residuals and covariates) for a given set of genes in single-nucleus RNA-seq data.
#'
#' @param gene_names A character vector specifying the names of genes for which variances need to be calculated.
#' @param df A data frame containing the relevant data, including expression values of genes and specified covariates.
#' @param covariates A character vector specifying the covariates (default: c("species", "celltype")) to include in the linear model.
#' @param verbose A logical value indicating whether to display a progress bar (default: FALSE).
#'
#' @return A data frame containing gene names and corresponding variances for residuals, covariate 1, and covariate 2.
#'
#' @details The function fits linear models for each gene using specified covariates and extracts the relevant ANOVA information to calculate variances.
#'
#' @examples
#' gene_names <- c("Gene1", "Gene2", "Gene3")
#' data_frame <- read.csv("snRNAseq_data.csv")
#' result <- calculate_gene_variances(gene_names, data_frame)
#'
#' @export
calculate_gene_variances <- function(gene_names, df, covariates = c("species", "celltype"), verbose = FALSE) {
    # Create empty vectors to store sum_sqs
    resid_variances <- vector("numeric", length(gene_names))
    X1_variances <- vector("numeric", length(gene_names))
    X2_variances <- vector("numeric", length(gene_names))

    # Create a progress bar if verbose is TRUE
    if (verbose) {
        pb <- txtProgressBar(min = 0, max = length(gene_names), style = 3)
    }

    # Loop through genes and fit linear models
    for (i in seq_along(gene_names)) {
        gene <- gene_names[i]
        formula <- as.formula(paste(gene, "~", covariates[1], "+", covariates[2]))
        model <- lm(formula, data = df)

        # Extract ANOVA table
        anova_table <- anova(model)

        # Extract relevant information
        sum_sq <- anova_table$`Sum Sq`

        # Store sum_sqs for each gene
        resid_variances[i] <- sum_sq[3]
        X2_variances[i] <- sum_sq[2]
        X1_variances[i] <- sum_sq[1]

        # Update the progress bar if verbose is TRUE
        if (verbose) {
            setTxtProgressBar(pb, i)
        }
    }

    # Close the progress bar if verbose is TRUE
    if (verbose) {
        close(pb)
    }

    # Create a data frame with sum_sqs
    sum_sqs_df <- data.frame(
        Gene = gene_names,
        Resid_Variance = resid_variances
    )

    colnames(sum_sqs_df)[2] <- covariates[1]
    colnames(sum_sqs_df)[3] <- covariates[2]

    return(sum_sqs_df)
}
