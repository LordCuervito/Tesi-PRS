import os
import pandas as pd

# Ruta a las carpetas de tus archivos
vcf_folder = r'C:\Users\Cuervito\Tesi\DataSets\OpenSNP'
gwas_file = r'C:\Users\Cuervito\Tesi\DataSets\FinnGen\summary_stats_finngen_R11_F5_ANXIETY-006_'


# Cargar los datos de FinnGen
gwas_data = pd.read_csv(gwas_file, sep='\t', low_memory=False)
gwas_data['rsids'] = gwas_data['rsids'].str.lower().str.strip()

# Obtener todos los archivos .vcf
vcf_files = [f for f in os.listdir(vcf_folder) if f.endswith('.vcf')]

# Crear una lista para almacenar los genotipos alineados
genotype_matrix = []

def procesar_archivo_vcf(user_file_path):
    try:
        # Cargar datos del archivo VCF
        user_data = pd.read_csv(user_file_path, sep='\t', comment='#', header=None,
                                names=['chrom', 'pos', 'ID', 'ref', 'alt', 'qual', 'filter', 'info', 'format', 'genotype'])
        user_data['ID'] = user_data['ID'].str.lower().str.strip()
        
        # Encontrar SNPs comunes con el archivo de GWAS
        common_snps = set(user_data['ID']).intersection(set(gwas_data['rsids']))
        if len(common_snps) == 0:
            print(f"Sin SNPs comunes en {os.path.basename(user_file_path)}.")
            return None

        # Filtrar y alinear los SNPs comunes
        user_genotypes = user_data[user_data['ID'].isin(common_snps)]
        gwas_snps = gwas_data[gwas_data['rsids'].isin(common_snps)]
        merged_data = pd.merge(user_genotypes, gwas_snps, left_on='ID', right_on='rsids')

        # Alineación de alelos y conversión de genotipos a numéricos
        mismatched_alleles = (merged_data['ref_x'] != merged_data['ref_y'])
        merged_data.loc[mismatched_alleles, 'genotype'] = merged_data.loc[mismatched_alleles, 'genotype'].replace({'0/1': '1/0', '1/0': '0/1'})
        merged_data['genotype_numeric'] = merged_data['genotype'].replace({'0/0': 0, '0/1': 1, '1/1': 2, './.': None})
        merged_data['genotype_numeric'] = pd.to_numeric(merged_data['genotype_numeric'], errors='coerce')

        # Agregar los genotipos numéricos a la matriz de genotipos
        genotype_matrix.append(merged_data['genotype_numeric'].values)
        return True

    except Exception as e:
        print(f"Error procesando {os.path.basename(user_file_path)}: {e}")
        return False

# Procesar todos los archivos VCF
for idx, user_file in enumerate(vcf_files):
    user_file_path = os.path.join(vcf_folder, user_file)
    procesar_archivo_vcf(user_file_path)
    if (idx + 1) % 100 == 0:  # Guardar cada 100 archivos procesados para evitar pérdida de datos
        genotype_matrix_df = pd.DataFrame(genotype_matrix).T.dropna()
        genotype_matrix_df.to_csv(f'C:\\Users\\Cuervito\\Tesi\\genotype_matrix_partial_{(idx+1)//100}.csv', index=False)

# Guardar la matriz completa
genotype_matrix_df = pd.DataFrame(genotype_matrix).T.dropna()
genotype_matrix_df.to_csv(r'C:\Users\Cuervito\Tesi\genotype_matrix.csv', index=False)

print("Matriz de genotipos guardada.")
