# Calcular PRS y guardar resultados por lote y al final el resultado completo

import pandas as pd

# Cargar la matriz de genotipos previamente guardada
genotype_matrix = pd.read_csv('/content/drive/My Drive/Tesi/genotype_matrix.csv')

# Cargar los datos GWAS de FinnGen
gwas_file = '/content/drive/My Drive/Tesi/zip_mental/Animic/summary_stats_finngen_R11_F5_ANXIETY'
gwas_data = pd.read_csv(gwas_file, sep='\t')

# Calcular PRS utilizando la matriz de genotipos y los betas del GWAS
genotype_matrix['weighted_beta'] = genotype_matrix.mul(gwas_data['beta'], axis=0)
prs_scores = genotype_matrix.sum(axis=1)

# Guardar los PRS calculados
prs_scores.to_csv('/content/drive/My Drive/Tesi/prs_scores.csv', index=False)

print("PRS calculados y guardados.")
