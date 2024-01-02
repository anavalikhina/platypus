from sklearn.cluster import KMeans
import seaborn as sns
from matplotlib import pyplot as plt
import pandas as pd
import scanpy as sc
import anndata as ad


expression_matrix = pd.read_csv("brain_counts.csv", index_col=0)
metadata = pd.read_csv("brain_metadata.csv", index_col=0)

# making annotated dataset
adata = ad.AnnData(expression_matrix, obs=metadata)

# calculating UMAP coordinates
sc.pp.neighbors(adata)
sc.tl.umap(adata, min_dist=0.7, spread=3.0, random_state=1, n_components=2)
umap_coordinates = pd.DataFrame(adata.obsm['X_umap'], columns=['component1', 'component2'])

# running K-Means clusterization on the UMAP coords
kmeans = KMeans(n_clusters=7, random_state=0).fit(umap_coordinates)
umap_coordinates['labels'] = kmeans.labels_
umap_coordinates['class'] = metadata['cell_ontology_class'].values

# visualizing results
plt.figure()
sns.scatterplot(data=umap_coordinates, x='component1', y='component2', style='labels', hue='class')
plt.show()


