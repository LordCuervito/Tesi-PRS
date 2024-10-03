# Filtrar los SNPs del archivo GWAS que coinciden con los SNPs del usuario
    user_snps = set(user_data['SNP'])  # Ajustar según el formato de tu archivo de usuario
    gwas_snps = gwas_data[gwas_data['SNP'].isin(user_snps)]

# Verificar los SNPs alineados
    print(gwas_snps.head())
