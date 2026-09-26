import pandas as pd
import requests
import time

class EnsemblData:
    # CONFIGURACIÓN
    Ensembl_URL = "https://rest.ensembl.org/vep/human/region"
    BATCH_SIZE = 100
    MAX_RETRIES = 3
    TIMEOUT = 120
    # Guardar el progreso cada X lotes
    SAVE_EVERY = 10
    ################################################
    # Funciones:
    # NORMALIZAR RSID
    def normalize_rsid(self, rsid):
        if pd.isna(rsid):
            return None
        rsid = str(rsid).strip()
        if rsid in ["", "-", ".", "nan", "None", "-1"]:
            return None
        if rsid.startswith("rs"):
            return rsid
        if rsid.isdigit():
            return "rs" + rsid
        return None
    # CREAR VARIANTE PARA Ensembl REST API:
    def create_variant_string(self, row):
        chromosome = str(row["Chromosome"])
        position = str(row["PositionVCF"])
        ref = str(row["ReferenceAlleleVCF"])
        alt = str(row["AlternateAlleleVCF"])
        rsid = self.normalize_rsid(
            row["RS# (dbSNP)"]
        )
        # Si existe RSID
        if rsid is not None:
            return (
                f"{chromosome} {position} {rsid} {ref} {alt}"
            )
        # Si no existe RSID, utilizar coordenadas
        return (
            f"{chromosome} {position} . {ref} {alt}"
        )
    # CONSULTAR Ensembl REST API
    def query_ensembl_data(self, variantes):
        headers = {
            "Content-Type": "application/json",
            "Accept": "application/json"
        }

        payload = {
            "variants": variantes
        }
        params = {
            "CADD": "snv" 
        }
        try_number = 1
        while True:
            try:
                response = requests.post(
                    self.Ensembl_URL,
                    headers=headers,
                    params=params,
                    json=payload,
                    timeout=self.TIMEOUT
                )
                # RESPUESTA CORRECTA
                if response.status_code == 200:
                    return response.json(), True
                # ERRORES TEMPORALES DEL SERVIDOR
                if response.status_code in [500,502,503,504]:
                    if try_number >= self.MAX_RETRIES:
                        time.sleep(300)
                        try_number = 1
                        continue
                    wait = 10 * try_number
                    time.sleep(wait)
                    try_number += 1
                    continue
                # OTROS ERRORES
                return [], False
            # TIMEOUT Y OTROS ERRORES DE CONEXIÓN
            except requests.exceptions.Timeout or requests.exceptions.RequestException:
                if try_number >= self.MAX_RETRIES:
                    time.sleep(300)
                    try_number = 1
                else:
                    wait = 10 * try_number
                    time.sleep(wait)
                    try_number += 1
    # EXTRAER CADD
    def get_cadd(self, record, alt):
        transcript_consequences = record.get(
            "transcript_consequences",
            []
        )
        valores = []
        for transcript in transcript_consequences:
            if transcript.get("variant_allele") != alt:
                continue
            cadd = transcript.get("cadd_phred")
            if cadd is not None:
                try:
                    valores.append(float(cadd))
                except (ValueError, TypeError):
                    pass

        if valores:
            return max(valores)

        return None
    # EXTRAER gnomAD_AF
    def get_gnomad_af(self, record, alt):
        colocated_variants = record.get("colocated_variants", [])

        for variant in colocated_variants:
            frequencies = variant.get("frequencies", {})

            if alt in frequencies:
                datos = frequencies[alt]

                # Primero intentamos gnomAD Genome
                gnomadg = datos.get("gnomadg")

                if gnomadg is not None:
                    return float(gnomadg)

                # Si no existe, usamos gnomAD Exome
                gnomade = datos.get("gnomade")

                if gnomade is not None:
                    return float(gnomade)

        return None
    ################################################
    # Función principal:
    def getCaddGnomead(self, input_file, output_file):
        df = pd.read_csv(input_file)
        # AÑADIR COLUMNAS DE FEATURES
        df["CADD"] = pd.NA
        df["gnomAD_AF"] = pd.NA
        # PROCESAMIENTO POR LOTES
        for start in range(0,len(df),self.BATCH_SIZE):
            end = min(start + self.BATCH_SIZE,len(df))
            batch = df.iloc[start:end]
            # Crear variantes
            variants = []
            for _, row in batch.iterrows():
                variants.append(self.create_variant_string(row))
            # Consultar Ensembl REST API
            response, success = self.query_ensembl_data(variants)
            # Procesar resultados
            for row_index, record in zip(batch.index,response):
                original_alt = str(df.loc[row_index, "AlternateAlleleVCF"])
                # CADD
                cadd = self.get_cadd(record,original_alt)
                # Si Ensembl REST API encuentra CADD:
                # sustituimos NA por el número.
                # Si no encuentra CADD:
                # se mantiene NA.
                if cadd is not None:
                    df.loc[row_index,"CADD"] = cadd
                # gnomAD_AF
                gnomad_af = self.get_gnomad_af(record, original_alt)
                if gnomad_af is not None:
                    df.loc[row_index, "gnomAD_AF"] = gnomad_af
                # Pausa entre lotes
                time.sleep(2)
                # GUARDADO FINAL
                df.to_csv(output_file, index=False, na_rep="NA")