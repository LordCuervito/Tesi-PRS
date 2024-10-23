#ARCHIVOS VCF OpenSNP

import os
import pandas as pd
from google.colab import drive

# Montar Google Drive
drive.mount('/content/drive', force_remount=True)

# Ruta a la carpeta de VCFs y los datos de FinnGen
vcf_folder = '/content/drive/My Drive/Tesi/OpenSNP/1-5k'
finngen_folder = '/content/drive/My Drive/Tesi/zip_mental/Animic'
gwas_file = os.path.join(finngen_folder, 'summary_stats_finngen_R11_F5_ANXIETY')

# Cargar los datos de FinnGen
gwas_data = pd.read_csv(gwas_file, sep='\t', low_memory=False)
gwas_data['rsids'] = gwas_data['rsids'].str.lower().str.strip()

# Obtener todos los archivos .vcf
vcf_files = [f for f in os.listdir(vcf_folder) if f.endswith('.vcf')]

# Dividir los archivos en lotes de 100
lote_size = 100
lotes_vcf = [vcf_files[i:i + lote_size] for i in range(0, len(vcf_files), lote_size)]

genotype_matrix = []  # Lista para guardar los genotipos alineados

def procesar_archivo_vcf(user_file_path):
    try:
        user_data = pd.read_csv(user_file_path, sep='\t', comment='#', header=None,
                                names=['chrom', 'pos', 'ID', 'ref', 'alt', 'qual', 'filter', 'info', 'format', 'genotype'])
        user_data['ID'] = user_data['ID'].str.lower().str.strip()
        common_snps = set(user_data['ID']).intersection(set(gwas_data['rsids']))

        if len(common_snps) == 0:
            print(f"Sin SNPs comunes en {os.path.basename(user_file_path)}.")
            return None

        # Filtrar solo los SNPs comunes y alinear
        user_genotypes = user_data[user_data['ID'].isin(common_snps)]
        gwas_snps = gwas_data[gwas_data['rsids'].isin(common_snps)]
        merged_data = pd.merge(user_genotypes, gwas_snps, left_on='ID', right_on='rsids')

        # Alineación de alelos
        mismatched_alleles = (merged_data['ref_x'] != merged_data['ref_y'])
        merged_data.loc[mismatched_alleles, 'genotype'] = merged_data.loc[mismatched_alleles, 'genotype'].replace({'0/1': '1/0', '1/0': '0/1'})

        # Convertir genotipos a valores numéricos
        merged_data['genotype_numeric'] = merged_data['genotype'].replace({'0/0': 0, '0/1': 1, '1/1': 2, './.': None})
        merged_data['genotype_numeric'] = pd.to_numeric(merged_data['genotype_numeric'], errors='coerce')

        # Guardar los genotipos numéricos alineados
        genotype_matrix.append(merged_data['genotype_numeric'].values)
        return True

    except Exception as e:
        print(f"Error procesando {os.path.basename(user_file_path)}: {e}")
        return False

# Procesar los archivos en lotes
for idx, lote in enumerate(lotes_vcf):
    print(f"\nProcesando lote {idx+1}/{len(lotes_vcf)} con {len(lote)} archivos...")
    for user_file in lote:
        user_file_path = os.path.join(vcf_folder, user_file)
        procesar_archivo_vcf(user_file_path)

# Guardar la matriz de genotipos procesada para usarla posteriormente
genotype_matrix_df = pd.DataFrame(genotype_matrix).T.dropna()
genotype_matrix_df.to_csv('/content/drive/My Drive/Tesi/genotype_matrix.csv', index=False)

print("Matriz de genotipos guardada.")


