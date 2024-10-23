# PCA y guardar resultados por lote y al final el resultado completo

import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

# Cargar la matriz de genotipos previamente guardada
genotype_matrix = pd.read_csv('/content/drive/My Drive/Tesi/genotype_matrix.csv')

# Aplicar PCA
pca = PCA(n_components=10)
pca_result = pca.fit_transform(genotype_matrix)

# Guardar los componentes principales
pca_components = pd.DataFrame(pca_result, columns=[f'PC{i}' for i in range(1, 11)])
pca_components.to_csv('/content/drive/My Drive/Tesi/pca_components.csv', index=False)

# Visualizar los dos primeros componentes
plt.figure(figsize=(8, 6))
plt.scatter(pca_components['PC1'], pca_components['PC2'], alpha=0.5)
plt.title('PCA - Componente Principal 1 vs Componente Principal 2')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.grid(True)
plt.show()
