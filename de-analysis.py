import pandas as pd
import rpy2.robjects as ro

from de_analysis.utils import convert_dataframe_to_python, convert_dataframe_to_r
from de_analysis.volcano import plot_volcano

r = ro.r
r.source("de_analysis/run_analysis.R")


expression_matrix = pd.read_csv("input_files/brain_counts.csv", index_col=0)
metadata = pd.read_csv("input_files/brain_metadata.csv", index_col=0)

results = r.run_da_analysis(
    count_matrix=convert_dataframe_to_r(expression_matrix.T),
    conditions=convert_dataframe_to_r(metadata['cell_ontology_class']),
    normal_condition='astrocyte',
    preprocess=True,
    method="limma",
    filename="results",
    save=True
)

da_results_limma = []
for name, df in results.items():
    da_results_limma.append(convert_dataframe_to_python(df))
    plot_volcano(convert_dataframe_to_python(df), filename=name)

