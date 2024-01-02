library(DESeq2)
source('de_analysis/utils.R')


run_deseq2_analysis <- function(arguments){

    conditions <- make_conditions_vector(arguments$conditions)

    preprocessed_data <- deseq2_preprocess(arguments$count_matrix,
                                           conditions,
                                           arguments$normal_condition,
                                           arguments$features,
                                           arguments$cutoff,
                                           arguments$features_id_colname)

    comparisons <- make_comparisons(conditions = conditions,
                                    normal_condition = arguments$normal_condition,
                                    contrast_condition = arguments$contrast_condition)
    # obtain DE results
    DEresults <- deseq2_analyse(preprocessed_data = preprocessed_data,
                                conditions = conditions,
                                comparisons = comparisons)

    return(DEresults)
}


deseq2_preprocess <- function(count_matrix, conditions, normal_condition, features, cutoff, features_id_colname){

    conditions <- factorize_conditions(conditions = conditions,
                                       normal_condition = normal_condition)
    filtered <- filter_lib(count_matrix = count_matrix, features = features, cutoff = cutoff,
            features_colname = features_id_colname)
    count_matrix <- filtered[[1]]

    preprocessed_data <- DESeqDataSetFromMatrix(countData = count_matrix,
                                                colData = data.frame(Treatment = conditions),
                                                design = formula(paste("~", "Treatment")))
    return(preprocessed_data)
}

deseq2_analyse <- function(preprocessed_data, conditions, comparisons) {

    #' Run DEResultColumnsLimma analysis
    #' @description
    #' Run deseq2 on Deseq data set
    #' @param preprocesssed_data Deseq data set
    #' @return
    #' Named list  with DataFrames containing DEResultColumnsLimma analysis results.
    #' Each name in the list is a contrast in the study

    design <- model_matrix(conditions = conditions)
    mod_mat <- design$design
    dds <- DESeq(preprocessed_data)
    DEresults <- calc_coefficients(dds = dds, mod_mat = mod_mat, comparisons = comparisons)
    return(DEresults)
}

calc_coefficients <- function(dds, mod_mat, comparisons) {
    #' calculate coefficients to get deseq results based on normal and contrast condition
    #' @param dds Deseq data set
    #' @param conditions vector with different conditions or treatments
    #' @param mod_mat matrix. Model design matrix
    #' @param comparisons vector with normal and contrast conditions
    #' @return Deseq results based on comparisons vector

    coefficients_list <- list()
    for (comparison in comparisons) {
        comparison <- unlist(strsplit(comparison, "-"))
        for (condition in comparison) {
            coefficients_results <- colMeans(mod_mat[dds[["Treatment"]] == condition,
                ])
            coefficients_list[[condition]] <- coefficients_results
        }
    }
    DEresults_coefficients <- list()
    for (comparison in comparisons) {
        comparison <- unlist(strsplit(comparison, "-"))
        condition1 <- comparison[1]
        condition2 <- comparison[2]
        DEresults_coefficients[[paste(condition1, "-", condition2)]] <- as.data.frame(results(dds,
            contrast = coefficients_list[[condition1]] - coefficients_list[[condition2]], ))
    }

    return(DEresults_coefficients)
}