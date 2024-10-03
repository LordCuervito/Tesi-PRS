#importar de drive y librerias

    from google.colab import drive
    drive.mount('/content/drive')

    import os
    import pandas as pd

# Ruta a la carpeta
    gwas_folder = '/content/drive/My Drive/Tesi/zip_mental/Animic/'

# Listar archivos en la carpeta
    files = os.listdir(gwas_folder)
    print(files)

# Ruta al archivo sin extensión
    gwas_file = '/content/drive/My Drive/Tesi/zip_mental/Animic/summary_stats_finngen_R11_F5_ANXIETY'

# Intentar cargar el archivo directamente
    gwas_data = pd.read_csv(gwas_file, sep='\t')  # Ajusta el separador si es necesario (por defecto es tabulación)

# Verificar las primeras filas para asegurarnos de que se cargó correctamente
    gwas_data.head()
