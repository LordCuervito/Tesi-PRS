# PCA

```python
import os
import pandas as pd
from sklearn.decomposition import PCA
from google.colab import drive

# Montar Google Drive
drive.mount('/content/drive', force_remount=True)

# Ruta a la carpeta de VCFs y datos GWAS
vcf_folder = '/content/drive/My Drive/Tesi/OpenSNP/1-5k'
gwas_file = '/content/drive/My Drive/Tesi/zip_mental/Animic/summary_stats_finngen_R11_F5_ANXIETY'

# Cargar los datos GWAS
gwas_data = pd.read_csv(gwas_file, sep='\t', usecols=['rsids', 'beta'], dtype={'rsids': str, 'beta': float})
gwas_data['rsids'] = gwas_data['rsids'].str.lower().str.strip()

# Obtener todos los archivos .vcf
vcf_files = [f for f in os.listdir(vcf_folder) if f.endswith('.vcf')]

# Dividir los archivos en lotes de 100
lote_size = 100
lotes_vcf = [vcf_files[i:i + lote_size] for i in range(0, len(vcf_files), lote_size)]

# Crear una lista para almacenar todos los genotipos
genotype_matrix = []

def procesar_archivo_vcf(user_file_path):
    try:
        # Cargar datos del archivo VCF
        user_data = pd.read_csv(user_file_path, sep='\t', comment='#', header=None, names=['chrom', 'pos', 'ID', 'ref', 'alt', 'qual', 'filter', 'info', 'format', 'genotype'])
        print(f"Procesando archivo: {os.path.basename(user_file_path)}")
        print(f"Número de columnas en {os.path.basename(user_file_path)}: {user_data.shape[1]}")

        # Convertir los genotipos a valores numéricos
        user_data['genotype_numeric'] = user_data['genotype'].replace({'0/0': 0, '0/1': 1, '1/1': 2, './.': None})
        user_data['genotype_numeric'] = pd.to_numeric(user_data['genotype_numeric'], errors='coerce')

        # Filtrar SNPs que están en el GWAS
        user_data_filtered = user_data[user_data['ID'].isin(gwas_data['rsids'])]

        # Añadir los genotipos numéricos a la matriz de genotipos
        genotype_matrix.append(user_data_filtered['genotype_numeric'].values)
        return True
    except Exception as e:
        print(f"Error procesando {os.path.basename(user_file_path)}: {e}")
        return False

# Procesar los archivos en lotes de 100
for idx, lote in enumerate(lotes_vcf):
    print(f"\nProcesando lote {idx+1}/{len(lotes_vcf)} con {len(lote)} archivos...")
    archivos_fallidos = []

    for user_file in lote:
        user_file_path = os.path.join(vcf_folder, user_file)
        if not procesar_archivo_vcf(user_file_path):
            archivos_fallidos.append(user_file)

    # Reintentar procesar los archivos que fallaron
    if archivos_fallidos:
        print(f"\nReintentando los {len(archivos_fallidos)} archivos fallidos del lote {idx+1}/{len(lotes_vcf)}...")
        for user_file in archivos_fallidos:
            user_file_path = os.path.join(vcf_folder, user_file)
            procesar_archivo_vcf(user_file_path)

# Convertir la matriz de genotipos a un DataFrame
genotype_matrix = pd.DataFrame(genotype_matrix).T.dropna()

# Aplicar PCA
pca = PCA(n_components=10)  # Puedes ajustar el número de componentes principales
pca_result = pca.fit_transform(genotype_matrix)

# Guardar los componentes principales para utilizarlos en el cálculo de PRS
pca_components = pd.DataFrame(pca_result, columns=[f'PC{i}' for i in range(1, 11)])

print("\nPCA completado. Los componentes principales han sido calculados.")
