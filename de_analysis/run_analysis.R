library(magrittr)
library(stringr)

source('de_analysis/limma.R')
source('de_analysis/deseq2.R')
source('de_analysis/utils.R')

run_da_analysis <- function(count_matrix, conditions, method = "limma", padjmethod = "BH",
                            features = NULL, features_id_colname = NULL, cutoff = NULL,
                            contrast_condition = regex("\\w"), normal_condition = NULL,
                            preprocess = TRUE, filename, save) {


    #' Differential Expression analysis
    #' @description
    #' Perform DEResultColumnsLimma or DA analysis for given microarrays or RNAseq experiment
    #' @param count_matrix matrix with sequence or species counts. Rows are features, columns are samples
    #' @param conditions data.frame. Contains data on phenotypes/diseases/treatments
    #' @param method character. Name of method to use. Options are "limma", "DEseq2".
    #' @param features data.frame. Contains feature data.
    #' @param features_id_colname character. Column if features with IDs corresponding to count_matrix row names.
    #' @param filename character. Study name which will be included in filenames of all results
    #' @param save logical. Whether saving results is needed

    #' @return DE results for each contrast in the study

    methods <- list(limma = run_limma_analysis,
                    DESeq2 = run_deseq2_analysis)
    common_arguments <- list(count_matrix = count_matrix,
                             conditions = conditions,
                             normal_condition = normal_condition,
                             contrast_condition = contrast_condition,
                             features = features,
                             cutoff = cutoff,
                             features_id_colname = features_id_colname)
    arguments <- list(limma = list(padjmethod = padjmethod,
                                   preprocess = preprocess),
                      DESeq2 = list())

    stopifnot(method %in% names(methods))
    DEresults <- methods[[method]](c(common_arguments, arguments[[method]]))

    if (save) {
        # save
        for (contrast in names(DEresults)) {
            save_de_results(DEresults = DEresults, contrast = contrast, filename = filename)
        }
    }
    return(DEresults)
}




