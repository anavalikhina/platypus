library(dplyr)

normalize_libs <- function(count_matrix) {
    #' Normalize libraries by size
    # initialize DGElist, counts, library size, normalization factors
    count_matrix <- edgeR::DGEList(count_matrix)

    # calculate normalization factors regarding library size
    count_matrix <- edgeR::calcNormFactors(count_matrix)

    return(count_matrix)
}

filter_lib <- function(count_matrix, features = NULL, features_colname = NULL, cutoff = NULL) {
    #' Filter low-expressed genes
    if (is.null(cutoff)) {
        cutoff <- 0
    }

    # max expression between samples must be > cutoff
    drop <- which(apply(edgeR::cpm(count_matrix), 1, max) < cutoff)
    if (length(drop) != 0) {
        count_matrix_filtered <- count_matrix[-drop, ]
    } else {
        count_matrix_filtered <- count_matrix
    }

    # filter genes and dge lists to have similar probes
    if (!is.null(features)){
        features <- features[features[, features_colname] %in% rownames(count_matrix_filtered$counts), ]
    }

    return(list(count_matrix_filtered = count_matrix_filtered, features = features))
}

make_conditions_vector <- function(conditions) {
    # create conditions vector
    if (class(conditions) != "character") {
        conditions <- conditions[, stringr::str_detect(string = tolower(colnames(conditions)),
            pattern = "factor|diagnosis|disease|treat|group")]
    }
    if ((class(conditions) != "character") & (class(conditions) != "factor")) {
        if (ncol(conditions) > 1) {
            conditions <- conditions[[1]]
        }
    }

    # modify conditions vector
    conditions <- sapply(tolower(conditions), stringr::str_replace_all, " ", "")

    return(conditions)
}

model_matrix <- function(conditions) {
    #' Design experiment
    #' @description
    #' Make model matrix from conditions vector
    #' @param conditions character. Vector with conditions
    #' @return named list with model design (matrix) and conditions (chracter)

    # build experiment design
    group <- factor(conditions, levels = unique(conditions))
    design <- model.matrix(~0 + group)
    colnames(design) <- levels(group)

    return(list(design = design, conditions = group))
}

make_comparisons <- function(conditions, normal_condition, contrast_condition) {

    contrast_condition <- tolower(contrast_condition)
    contrast_condition <- sapply(contrast_condition, stringr::str_replace_all, " ",
        "")

    # make contrasts
    conditions_unique <- unique(conditions)
    normal_condition <- levels(factorize_conditions(conditions = conditions,
                                                    normal_condition = normal_condition))[[1]]

    diseases <- conditions_unique[conditions_unique != normal_condition]
    diseases <- diseases[stringr::str_detect(string = stringr::fixed(diseases, ignore_case = TRUE),
        pattern = contrast_condition)]

    comparisons <- NULL
    for (disease in diseases) {
        comparisons <- c(comparisons, paste(disease, normal_condition, sep = "-"))
    }
    return(comparisons)
}


factorize_conditions <- function(conditions, normal_condition){
    require(stringr)

    conditions <- unlist(lapply(conditions, stringr::str_to_lower))
    conditions_unique <- unique(conditions)

    if (is.null(normal_condition)) {
        finder <- stringr::str_detect(conditions_unique, "normal|control|healthy")
        normal_condition <- conditions_unique[finder]

    } else if (!is.null(normal_condition)) {
        normal_condition <- tolower(normal_condition)
    }
    # WHAT IF THERE'S NO NORMAL/CONTROL/HEALTHY - SET NORMAL TO THE FIRST CONDITION IN ALPHABETICAL ORDER
    if (length(normal_condition) == 0) {
        conditions <- factor(conditions)
        warning(paste0("Setting normal condition to ", levels(conditions)[[1]]))
    } else {
      conditions <- relevel(factor(conditions), ref = normal_condition)
    }
    return(conditions)
}


select_comparison_data <- function(count_matrix, conditions, comparison){
    contrast_condition <- unlist(strsplit(comparison, "-"))[[1]]
    normal_condition <- unlist(strsplit(comparison, "-"))[[2]]
    idx <- conditions %in% c(normal_condition, contrast_condition)
    count_matrix_selected <- count_matrix[, idx]
    conditions_selected <- conditions[idx]
    return(list(count_matrix_selected = count_matrix_selected,
                conditions_selected = conditions_selected))
}

save_de_results <- function(DEresults, filename, contrast) {

    #' Save results
    #' @description
    #' Save results
    #' @param DEresults
    #' @param filename
    #' @param contrast

    results <- DEresults[[contrast]]
    if (!("feature" %in% colnames(results))){
        results[["feature"]] <- rownames(results)
        rownames(results) <- NULL
    }
    write.csv(as.matrix(results), paste0(paste("DE", contrast, filename,
        sep = "_"), ".csv"))
}
