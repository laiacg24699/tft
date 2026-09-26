import pandas as pd
from sklearn.model_selection import train_test_split

class filter:
    # CONFIGURACIÓN
    acmg_classes = [
        "Pathogenic",
        "Likely pathogenic",
        "Uncertain significance",
        "Likely benign",
        "Benign"
    ]
    required_columns = [
        "Chromosome",
        "Start",
        "Stop",
        "ReferenceAllele",
        "AlternateAllele"
    ]
    valid_types = [
        "single nucleotide variant",
    ]
    invalid_values = ["-", ".", ""]
    invalid_values_vcf = ["-", ".", "", "na"]
    chunksize = 10000
    n_per_class = 10000
    ################################################
    # Funciones:
    # APLICA FILTROS
    def clean(self, chunk):
        #Filtrar variantes dejando solo aquellas cuya ClinicalSignificance sea benigna, probablemente benigna, de significado incierto, probablemente patogénica o patogénica. 
        mask_acmg = chunk["ClinicalSignificance"].isin(self.acmg_classes)
        #Filtrar dejando solo variantes que sean SNV:
        mask_type = chunk["Type"].isin(self.valid_types)
        #Filtrar variantes dejando solo aquellas que estén en assembly GRCh38:
        mask_assembly = chunk["Assembly"] == "GRCh38"
        #Filtrar variantes cuyos alelos de referencia y alternativo, inicio y final, o cromosoma sean nulos.
        mask_not_null = chunk[self.required_columns].notna().all(axis=1)
        #Filtrar variantes con posiciones VCF, alelos de referencia VCF o alelos alternaticos VCF con valores nulos, o bien iguales a "-", ".", "" o "na".
        mask_valid_alleles = (
            chunk["PositionVCF"].notna() &
            chunk["ReferenceAlleleVCF"].notna() &
            chunk["AlternateAlleleVCF"].notna() &
            ~chunk["PositionVCF"].astype(str).str.lower().isin(self.invalid_values_vcf) &
            ~chunk["ReferenceAlleleVCF"].astype(str).str.lower().isin(self.invalid_values_vcf) &
            ~chunk["AlternateAlleleVCF"].astype(str).str.lower().isin(self.invalid_values_vcf)
        )
        #Filtrar variantes cuyas columnas de alelo de referencia y alternativo, inicio y final, y del cromosoma tengan alguno de los valores "-", "." o "".
        mask_valid_values = ~chunk[self.required_columns].isin(self.invalid_values).any(axis=1)
        #Leer valor numérico de start:
        start = pd.to_numeric(
            chunk["Start"],
            errors="coerce"
        )
        #Leer valor numérico de stop:
        stop = pd.to_numeric(
            chunk["Stop"],
            errors="coerce"
        )
        #Filtrar variantes con valores de start y stop mayores a 0 y con start igual o inferior a stop, para filtrar datos incoherentes:
        mask_coordinates = (
            (start > 0) &
            (stop > 0) &
            (start <= stop)
        )
        #Devolver chunk con todos los filtros aplicados: 
        return (
            mask_acmg &
            mask_type &
            mask_assembly &
            mask_not_null &
            mask_valid_alleles &
            mask_valid_values &
            mask_coordinates
        )
    # ELIMINA VARIANTES CON DATOS CONFLICTIOS
    def consolidation(self, df):
        #Agrupar filas del dataframe por AlleleId, seleccionar la columna "ClinicalSignificance" de cada alelo y contar cuantos valores distintos hay.
        n_classes = (
            df
            .groupby("#AlleleID")["ClinicalSignificance"]
            .nunique()
        )
        # Obtener los AlleleIDs de los alelos cuyas variantes tengan un único valor de "ClinicalSignificance":
        non_conflicting_ids = n_classes[
            n_classes == 1
        ].index
        # Generar un dataframe con las variantes cuyo AllelleID se encuentre en la lista "non_conflicting_ids".
        consolidated_df = df[
            df["#AlleleID"].isin(non_conflicting_ids)
        ].copy()
        # Borrar variantes duplicadas (que tengan el mismo AlleleID):
        consolidated_df = consolidated_df.drop_duplicates(
            subset="#AlleleID"
        ).reset_index(drop=True)
        # Devolver dataframe filtrado:
        return consolidated_df
    # DEVUELVE UN DATAFRAME CON LAS VARIANTES AGRUPADAS POR CLÍNICAL SIGNIFICANCE Y CON NÚMERO ESPECIFICADO DE VARIANTES EN CADA GRUPO
    def prepare_df(self, df):
        return (
            df
            .groupby("ClinicalSignificance", group_keys=False) # Agrupar las variantes por "ClinicalSignificance".
            .sample(
                n=min(self.n_per_class, len(df)), # Para cada grupo seleccionar "n_per_class" valores o bien todos los elementos del grupo (tantos elementos como la longitud del grupo) en el caso de no llegar a ese valor.
                random_state=42 # Aplicar una semilla para la selección aleatoria de variantes (usar esta semilla permite obtener siempre las mismas variantes aleatorias).
            )
            .reset_index(drop=True) #Resetear los índices del dataframe para el nuevo dataframe obtenido, que incorpora las variantes seleccionadas anteriormente de cada grupo.
        )
    # OBTIENE EL DATAFRAME FILTRADO LISTO PARA GUARDAR
    def get_consolidated_dataframe(self, input_file):
        # Crear lista de datos filtrados:
        filtered_chunks = []
        # Leer archivo por tramos:
        for chunk in pd.read_csv(input_file, sep="\t", compression="gzip", chunksize=self.chunksize):
            # Crear variable con las condiciones de filtraje:
            mask = self.clean(chunk)
            # Filtrar chunk utilizando las condiciones de la variable "mask":
            filtered_chunk = chunk.loc[mask]
            # Añadir los datos del chunk que han pasado el filtro a la lista de resultados:
            filtered_chunks.append(filtered_chunk)
        # Generar un dataframe concatenando los datos de la lista de resultados.
        filtered_df = pd.concat(
            filtered_chunks,
            ignore_index=True
        )
        # Aplicar la función de consolidación
        consolidated_df = self.consolidation(filtered_df)
        return self.prepare_df(consolidated_df)
    #################################################
    # Función principal:
    def filter_data(self, input_file, csv_file):
        #Obtener dataframe filtrado:
        consolidated_df = self.get_consolidated_dataframe(input_file)
        #guardar dataframe filtrado. 
        consolidated_df.to_csv(csv_file, index=False)