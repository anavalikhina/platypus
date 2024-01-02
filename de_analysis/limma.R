library(limma)

source('de_analysis/utils.R')

run_limma_analysis <- function(arguments){

    count_matrix <- as.matrix(arguments$count_matrix)

    if (is.null(arguments$features)){
        features <- data.frame(features = row.names(count_matrix))
        features_id_colname <- "features"
    } else{
        features <- arguments$features
        features_id_colname <- arguments$features_id_colname
    }

    conditions <- make_conditions_vector(arguments$conditions)

    if (arguments$preprocess){
        preprocessed_data <- limma_preprocessing(count_matrix = count_matrix,
                                                 features = features,
                                                 features_id_colname = features_id_colname,
                                                 conditions = conditions,
                                                 cutoff = arguments$cutoff)
        EList <- preprocessed_data$EList
        features <- preprocessed_data$features
    }else{
        EList <- ExpressionSet(assayData = count_matrix,
                               phenoData = AnnotatedDataFrame(data.frame(arguments$conditions,
                                                                         row.names = colnames(count_matrix))),
                               featureData = AnnotatedDataFrame(data.frame(arguments$features,
                                                                           row.names = rownames(count_matrix))))
    }

    design <- model_matrix(conditions = conditions)
    conditions <- design$conditions
    design <- design$design

    bayes_fit <- limma_fit(AEset = EList,
                           design = design,
                           conditions = conditions,
                           contrast_condition = arguments$contrast_condition,
                           normal_condition = arguments$normal_condition,
                           features = features)

    # obtain DEResultColumnsLimma results
    DEresults <- limma_analyse(bayes_fit = bayes_fit, padjmethod = arguments$padjmethod)
    return(DEresults)
}


limma_preprocessing <- function(count_matrix, features, conditions, cutoff, features_id_colname){

    #' limma preprocessing
    #' @description
    #' limma data preprocessing pipeline.
    #' Normalize libraries by size and
    #' perform filtering for low expression and variation.
    #' Then apply voom on expression data to feed it into limma linear models in analysis pipe.

    count_matrix <- normalize_libs(count_matrix = count_matrix)
    filtered <- filter_lib(count_matrix = count_matrix, features = features, cutoff = cutoff,
            features_colname = features_id_colname)
    count_matrix <- filtered[[1]]
    features <- filtered[[2]]
    y <- voom_dge(count_matrix = count_matrix, conditions = conditions)

    preprocessed_data <- (list(EList = y, features = features))
    return (preprocessed_data)
}

voom_dge <- function(count_matrix, conditions) {
    #' Apply voom
    design <- model_matrix(conditions = conditions)
    design <- design$design

    y <- limma::voom(count_matrix, design, plot = FALSE)
    return(y)
}

limma_fit <- function(AEset, design, conditions, contrast_condition, normal_condition,
                      features) {

    #' Fit model
    #' @description
    #' Fit linear model, contrast matrix, and Bayes model
    #' @param AEset ExpressionSet or EList
    #' @param design matrix. Model design matrix
    #' @param conditions character. Vector with condition
    #' @param normal_condition character. Name of condition which will be in the denominator of logFC values.
    #' @param contrast_condition character. Name of condition which will be in the nominator of logFC values.
    #' @param annotations data.frame. Contains gene data in columns 'PROBEID', 'ENTREZID', 'SYMBOL'
    #' @return MArrayLM object. Bayes model fitted into expression data

    # fit linear model
    data.fit <- limma::lmFit(AEset, design)

    # make contrast matrix
    comparisons <- make_comparisons(conditions = conditions, normal_condition = normal_condition,
        contrast_condition = contrast_condition)
    contrast.matrix <- limma::makeContrasts(contrasts = comparisons, levels = design)

    # fit contrast
    data.fit.con <- limma::contrasts.fit(data.fit, contrast.matrix)

    # fit Bayes
    data.fit.eb <- limma::eBayes(data.fit.con)

    # annotate
    if (!is.null(features)) {
        data.fit.eb$genes <- features
    }

    return(data.fit.eb)
}

limma_analyse <- function(bayes_fit, padjmethod) {

    #' Run DEResultColumnsLimma analysis
    #' @description
    #' Run decision test on Bayes model fit
    #' @param bayes_fit MArrayLM. Bayes model fitted into expression data
    #' @return
    #' Named list  with DataFrames containing DEResultColumnsLimma analysis results.
    #' Each name in the list is a contrast in the study

    ngenes <- dim(bayes_fit$t)[1]
    DEresults <- list()
    for (coef in (1:dim(bayes_fit$contrasts)[2])) {
        stats <- limma::topTable(bayes_fit, adjust = padjmethod, coef = coef, sort.by = "logFC",
            number = ngenes)

        DEresults[[colnames(bayes_fit$contrasts)[coef]]] <- as.data.frame(stats)
    }

    return(DEresults)
}