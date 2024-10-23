# Calcular PRS
import os
import pandas as pd
from google.colab import drive

# Montar Google Drive
drive.mount('/content/drive', force_remount=True)

# Ruta a la carpeta de VCFs y datos GWAS
vcf_folder = '/content/drive/My Drive/Tesi/OpenSNP/1-5k'
gwas_file = '/content/drive/My Drive/Tesi/zip_mental/Animic/summary_stats_finngen_R11_F5_ANXIETY'  # Ajusta la ruta si es necesario

# Cargar los datos GWAS
gwas_data = pd.read_csv(gwas_file, sep='\t', usecols=['rsids', 'beta'], dtype={'rsids': str, 'beta': float})
gwas_data['rsids'] = gwas_data['rsids'].str.lower().str.strip()

# Obtener todos los archivos .vcf
vcf_files = [f for f in os.listdir(vcf_folder) if f.endswith('.vcf')]

# Dividir los archivos en lotes de 100
lote_size = 100
lotes_vcf = [vcf_files[i:i + lote_size] for i in range(0, len(vcf_files), lote_size)]

def calcular_prs(user_genotypes, gwas_data):
    # Filtrar SNPs comunes entre el genotipo del usuario y el GWAS
    common_snps = set(user_genotypes['ID']).intersection(set(gwas_data['rsids']))

    # Filtrar solo los SNPs comunes para el análisis
    user_genotypes = user_genotypes[user_genotypes['ID'].isin(common_snps)]
    gwas_snps = gwas_data[gwas_data['rsids'].isin(common_snps)]

    # Unir por la columna 'ID' o 'rsids' para emparejar los SNPs
    merged_data = pd.merge(user_genotypes, gwas_snps, left_on='ID', right_on='rsids')

    # Asegurarnos de que los valores de genotipo sean numéricos
    merged_data['genotype_numeric'] = merged_data['genotype'].replace({'0/0': 0, '0/1': 1, '1/1': 2, './.': None})

    # Conversión a numérico, y manejar valores no numéricos
    merged_data['genotype_numeric'] = pd.to_numeric(merged_data['genotype_numeric'], errors='coerce')  # Convertir a numérico
    merged_data['beta'] = pd.to_numeric(merged_data['beta'], errors='coerce')  # Asegurarse de que los valores de beta sean numéricos

    # Verificar si hay valores inválidos
    if merged_data['genotype_numeric'].isnull().any() or merged_data['beta'].isnull().any():
        print(f"Valores inválidos detectados en genotipo o beta. Los SNPs afectados serán excluidos.")

    # Filtrar filas donde los valores de genotipo o beta no sean válidos
    merged_data.dropna(subset=['genotype_numeric', 'beta'], inplace=True)

    if merged_data.empty:
        print("No hay datos válidos después del filtrado.")
        return 0

    # Calcular el efecto ponderado de cada SNP
    merged_data['weighted_beta'] = merged_data['genotype_numeric'] * merged_data['beta']

    # Calcular el PRS sumando los betas ponderados
    prs = merged_data['weighted_beta'].sum()

    return prs

def procesar_archivo_vcf(user_file_path):
    try:
        # Cargar datos del archivo VCF
        user_data = pd.read_csv(user_file_path, sep='\t', comment='#', header=None, names=['chrom', 'pos', 'ID', 'ref', 'alt', 'qual', 'filter', 'info', 'format', 'genotype'])
        print(f"Procesando archivo: {os.path.basename(user_file_path)}")
        print(f"Número de columnas en {os.path.basename(user_file_path)}: {user_data.shape[1]}")

        # Validar que las columnas sean las esperadas
        if 'ID' not in user_data.columns or 'genotype' not in user_data.columns:
            print(f"Archivo {os.path.basename(user_file_path)} no tiene las columnas esperadas.")
            return False

        # Calcular PRS para este archivo
        prs = calcular_prs(user_data, gwas_data)  # Llamada al cálculo de PRS real
        print(f"PRS calculado para {os.path.basename(user_file_path)}: {prs}")
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

print("\nProcesamiento de todos los lotes completado.")
