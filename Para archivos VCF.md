#importar archivos/librerias
    
    from google.colab import files

    from google.colab import drive
    drive.mount('/content/drive')

    import pandas as pd
    import os

# Abrir el archivo y leer las primeras líneas
    with open(user_file, 'r') as f:
        for i in range(10):
            print(f.readline())

    !pip install cyvcf2

    from cyvcf2 import VCF
    import os

# Ruta a la carpeta donde subiste los archivos VCF en Google Drive
    vcf_folder = '/content/drive/My Drive/Tesi/OpenSNP/vcf/'

# Listar todos los archivos que contienen "vcf" en el nombre y terminan en ".txt"
    vcf_files = [f for f in os.listdir(vcf_folder) if 'vcf' in f and f.endswith('.txt')]

# Mostrar cuántos archivos VCF se encontraron
    print(f'Se encontraron {len(vcf_files)} archivos VCF en formato .txt.')

# Procesar cada archivo VCF, ignorando aquellos que no se puedan abrir
    valid_vcf_files = []
        for vcf_file in vcf_files:
            vcf_path = os.path.join(vcf_folder, vcf_file)

    try:
        vcf = VCF(vcf_path)
# Mostrar la primera variante de cada archivo como ejemplo
        for variant in vcf:
            print(f"Archivo {vcf_file} procesado correctamente. Primera variante: {variant}")
            break  # Mostrar solo la primera variante

        valid_vcf_files.append(vcf_file)  # Guardar los archivos válidos para análisis posterior

    except Exception as e:
        print(f"Error abriendo el archivo {vcf_file}: {e}")
        continue  # Saltar archivos que no se pueden abrir

# Mostrar cuántos archivos fueron válidos
    print(f'Se procesaron correctamente {len(valid_vcf_files)} archivos VCF.')






