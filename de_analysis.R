library("reticulate")
source("de_analysis/run_analysis.R")
source_python("de_analysis/volcano.py")

expression <- read.csv("input_files/brain_counts.csv", row.names = 1)

targets <- read.csv("input_files/brain_metadata.csv", row.names = 1)
conditions <- targets[, 'cell_ontology_class']


# ANALYSE
de_results <- run_da_analysis(count_matrix = t(as.matrix(expression)),
                              conditions = conditions,
                              normal_condition = "astrocyte",
                              preprocess = TRUE,
                              save = FALSE)
for (name in names(de_results)){
  plot_volcano(df = de_results[[name]], filename = name)
}
