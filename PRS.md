#Calcular PRS
import os
import pandas as pd
from google.colab import drive

# Montar Google Drive
drive.mount('/content/drive', force_remount=True)

# Ruta a la carpeta de VCFs y los datos de FinnGen (sin compresión, según la información proporcionada)
vcf_folder = '/content/drive/My Drive/Tesi/OpenSNP/1-5k'
finngen_folder = '/content/drive/My Drive/Tesi/zip_mental/Animic'
gwas_file = os.path.join(finngen_folder, 'summary_stats_finngen_R11_F5_ANXIETY')  # Archivo GWAS sin compresión

# Cargar los datos de FinnGen
gwas_data = pd.read_csv(gwas_file, sep='\t', low_memory=False)
gwas_data['rsids'] = gwas_data['rsids'].str.lower().str.strip()  # Limpiar IDs de SNP

# Obtener todos los archivos .vcf
vcf_files = [f for f in os.listdir(vcf_folder) if f.endswith('.vcf')]

# Dividir los archivos en lotes de 100 para procesarlos en partes
lote_size = 100
lotes_vcf = [vcf_files[i:i + lote_size] for i in range(0, len(vcf_files), lote_size)]

def procesar_archivo_vcf(user_file_path):
    try:
        # Cargar datos del archivo VCF
        user_data = pd.read_csv(user_file_path, sep='\t', comment='#', header=None, names=['chrom', 'pos', 'ID', 'ref', 'alt', 'qual', 'filter', 'info', 'format', 'genotype'])
        print(f"Procesando archivo: {os.path.basename(user_file_path)}")
        print(f"Número de columnas en {os.path.basename(user_file_path)}: {user_data.shape[1]}")
        
        # Filtrar SNPs comunes
        user_data['ID'] = user_data['ID'].str.lower().str.strip()
        common_snps = set(user_data['ID']).intersection(set(gwas_data['rsids']))
        
        if len(common_snps) == 0:
            print(f"Sin SNPs comunes en {os.path.basename(user_file_path)}.")
            return 0
        
        # Filtrar solo los SNPs comunes
        user_genotypes = user_data[user_data['ID'].isin(common_snps)]
        gwas_snps = gwas_data[gwas_data['rsids'].isin(common_snps)]
        
        # Asegurarse de que estén alineados los SNPs por ID
        merged_data = pd.merge(user_genotypes, gwas_snps, left_on='ID', right_on='rsids')

        # Alineación de alelos
        mismatched_alleles = (merged_data['ref_x'] != merged_data['ref_y'])
        merged_data.loc[mismatched_alleles, 'genotype'] = merged_data.loc[mismatched_alleles, 'genotype'].replace({'0/1': '1/0', '1/0': '0/1'})
        
        # Convertir genotipos a valores numéricos
        merged_data['genotype_numeric'] = merged_data['genotype'].replace({'0/0': 0, '0/1': 1, '1/1': 2, './.': None})
        merged_data['genotype_numeric'] = pd.to_numeric(merged_data['genotype_numeric'], errors='coerce')
        
        # Calcular el PRS
        merged_data['weighted_beta'] = merged_data['genotype_numeric'] * merged_data['beta']
        prs = merged_data['weighted_beta'].sum()
        
        print(f"PRS calculado para {os.path.basename(user_file_path)}: {prs}")
        return prs

    except Exception as e:
        print(f"Error procesando {os.path.basename(user_file_path)}: {e}")
        return None

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
