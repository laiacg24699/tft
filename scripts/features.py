import pandas as pd
import myvariant
import re
import requests
import ast
import json

class features:

    def get_variant_features(self, variant_str, genome_build="grch38"):
        # Añadir información para consultar la base de datos d'Ensembl
        server = "https://rest.ensembl.org"
        ext = f"/vep/human/region/{variant_str}/1?"
        
        headers = {
            "Content-Type": "application/json",
            "Accept": "application/json"
        }
        
        # Activar complementos necesarios (dbNSFP incluye REVEL y PhyloP)
        params = {
            "CADD": 1,
            "gnomADe": 1,      # gnomAD Exomes
            "gnomADg": 1,      # gnomAD Genomes
            "dbNSFP": "REVEL_score,phyloP100way_vertebrate"
        }
        # Enivar petición:
        response = requests.get(server + ext, headers=headers, params=params)
        # Petición fallida.
        if not response.ok:
            response.raise_for_status()
            return None
        # Si la respuesta es exitosa, se deuelve la respuesta como objeto json.
        return response.json()

    def normalize_rsid(self, value):
        # Si rsid no tiene valor, devolver pd.NA (valor ausente).
        if pd.isna(value):
            return pd.NA
        # Eliminar espacios al principio y al final del rsid.
        value = str(value).strip()
        # Si el valor coincide con alguno de los valores de la lista, devolver pd.NA (valor ausente).
        if value in ["", "-", "na", ".", "NA", "-1", "nan"]:
            return pd.NA
        # Si el valor empieza con "rs", devolver el rsid, puesto que tiene el formato correcto
        if value.startswith("rs"):
            return value
        # Si no se han cumplido las condiciones anteriores:
        try: #Convertir valor a float. Convertir el float a un número entero y devolver un sting con el prefijo "rs" más ese número entero.
            return f"rs{int(float(value))}"
        except (ValueError, TypeError): # En caso de fallar la conversión, devolver pd.NA.
            return pd.NA

    def extract_ref_alt(self, variant_id):
        # Si el variant id expresado en rsid no es un string, devolver valores nulos para los alelos de referencia y alternativo.
        if not isinstance(variant_id, str):
            return None, None
        # Buscar estructura ":g." + número + "A", "C", "G" o "T" (grupo 1) + ">" + "A", "C", "G" o "T" (grupo 2):
        match = re.search(
            r":g\.\d+([ACGT]+)>([ACGT]+)",
            variant_id
        )
        # Si se encuentra la estructura, devoler el caracter del grupo 1 y 2 (son los alelos de referencia y alternativo respectivamente): 
        if match:
            return match.group(1), match.group(2)
        # Si no se encuentra la estructura, devolver valores nulos para los alelos de referencia y alternativo:
        return None, None

    def arrange_revel(self, revel):
        # Si revel es una lista o una tupla:
        if isinstance(revel, (list, tuple)):
            # Elimina valores nulos y NaN de revel, y convierte el resto de valores a float y los añade a la lista "valid_floats".
            valid_floats = [
                float(x) for x in revel if x is not None and str(x) != "nan"
            ]
            # Devuelve el valor de revel máximo de la lista anterior, o devuelve un valor nulo si la lista "valid_floats" está vacía: 
            return max(valid_floats) if valid_floats else None

        # Si revel es una cadena de texto:
        elif isinstance(revel, str):
            try:
                # Eliminar los espacios al inicio y al final de la cadena de texto del revel y convertir la cadena resultante a un objeto de python
                parsed = ast.literal_eval(revel.strip())
                # Si objeto resultante es una lista, una tupla o un set:
                if isinstance(parsed, (list, tuple, set)):
                    # Elimina valores nulos y NaN, y convierte el resto de valores a float y los añade a la lista "valid_floats", que contendrá un conjunto de revel.
                    valid_floats = [
                        float(x)
                        for x in parsed
                        if x is not None and str(x) != "nan"
                    ]
                    # Devuelve el valor de revel máximo de la lista anterior, o devuelve un valor nulo si la lista "valid_floats" está vacía: 
                    return max(valid_floats) if valid_floats else None
                # Si el objeto resultante no es una lista, una tupla ni una cadena de texto:
                else:
                    # Devolver el objeto resultante, que corresponde al valor de revel, convertido a float:
                    return float(parsed)
            # Si el código no puede procesarse:
            except (ValueError, SyntaxError):
                # Devolver valor nulo:
                return None

        # Si revel es un float o un número entero:
        elif isinstance(revel, (float, int)):
            # Devolver el revel convertido a float:
            return float(revel)
        # Si revel no es una lista, una tupla, una cadena de texto, un número entero ni un float:
        else:
            # Devolver valor nulo:
            return None

    def fill_annotations(self, df, matching_rows, record):
        for idx in matching_rows.index:
            # Recuperar alelo de referencia:
            ref = df.at[idx, "ReferenceAlleleVCF"]
            # Recuperar alelo alternativo:
            alt = df.at[idx, "AlternateAlleleVCF"]
            # Obtener diccionario de anotaciones con los valores de CADD, gnomAD_AF, REVEL y PhyloP:
            annotations = self.get_annotation(record, ref, alt)
            # Si el valor "CADD" del diccionario de anotaciones no es nulo:
            if annotations["CADD"] is not None:
                # Rellenar la columna "CADD" de la fila con índice "idx" del dataframe con el valor "CADD" del diccionario de anotaciones:
                df.at[idx, "CADD"] = annotations["CADD"]
            # Si el valor "gnomAD_AF" del diccionario de anotaciones no es nulo:
            if annotations["gnomAD_AF"] is not None:
                # Rellenar la columna "gnomAD_AF" de la fila con índice "idx" del dataframe con el valor "gnomAD_AF" del diccionario de anotaciones:
                df.at[idx, "gnomAD_AF"] = annotations[
                    "gnomAD_AF"
                ]
            # Si el valor "REVEL"" del diccionario de anotaciones no es nulo:
            if annotations["REVEL"] is not None:
                # Rellenar la columna "REVEL" de la fila con índice "idx" del dataframe con el valor "REVEL" del diccionario de anotaciones, aplicando una transformacion a ese alor antes:
                df.at[idx, "REVEL"] = self.arrange_revel(annotations["REVEL"])
            # Si el valor "PhyloP"" del diccionario de anotaciones no es nulo:
            if annotations["PhyloP"] is not None:
                # Rellenar la columna "PhyloP" de la fila con índice "idx" del dataframe con el valor "PhyloP" del diccionario de anotaciones:
                df.at[idx, "PhyloP"] = annotations["PhyloP"]

    def get_annotation(self, record, ref, alt):
        # Extracción de CADD y de PhyloP
        # Extraer variant id expresado en rsid de los datos de Myvariant::
        variant_id = record.get("_id", "")
        # Extraer alelo alternativo y de referencia del ariant id para validar con los obtenidos desde Myvariant:
        record_ref, record_alt = self.extract_ref_alt(variant_id)
        # Si los alelos de referancia y alternativo extraídos del variant id no son nulos:
        if ref and alt:
            # Comprovar alguno de los alelos extraídos del variant id no corresponde con el alelo correspondiente extraído de Myvariant:
            if (
                record_ref != str(ref).upper()
                or record_alt != str(alt).upper()
            ):
                # Devolver CADD, gnomAD_AF, REVEL y PhyloP nulos.
                return {
                    "CADD": None,
                    "gnomAD_AF": None,
                    "REVEL": None,
                    "PhyloP": None,
                }
        # Llenar diccionario de anotaciones con CADD, gnomAD_AF, REVEL y PhyloP nulos.
        annotations = {
            "CADD": None,
            "gnomAD_AF": None,
            "REVEL": None,
            "PhyloP": None,
        }

        #Extraer CADD de los datos de Myvariant:
        cadd = record.get("cadd")

        # Definir un diccionario para almacenar valores de CADD obtenido de Myariant:
        cadd_dict = {}
        # Si el CADD extraido de Myvariant es una lista no vacía, almacenar primer elemento de la lista en el diccionario cadd_dict:
        if isinstance(cadd, list) and len(cadd) > 0:
            cadd_dict = cadd[0]
        # Si el CADD extraido de Myvariant es un diccionario, almacenar su valor directamente en el diccionario cadd_dict:
        elif isinstance(cadd, dict):
            cadd_dict = cadd
        # Si se ha podido almacenar información en el diccionario cadd_dict:
        if cadd_dict:
            # Almacenar valor "phred" del diccionario cadd_dict en la columna CADD del diccionario de anotaciones.
            annotations["CADD"] = cadd_dict.get("phred")
            # Buscar phyplop dentro del diccionario cadd_dict:
            cadd_phylop = cadd_dict.get("phylop", {})
            # Si el phyplop encontrado es un diccionario:
            if isinstance(cadd_phylop, dict):
                # Obtener valor de la propiedad "vertebrate" del diccionario phyplop, o el valor de la propiedad "mammalian" en su defecto.
                annotations["PhyloP"] = cadd_phylop.get(
                    "vertebrate"
                ) or cadd_phylop.get("mammalian")

        # Extracción de REVEL, y de PhyloP a partir del dbnsfp (si el método anterior ha fallado).
        # Obtener dbnsfp:
        dbnsfp = record.get("dbnsfp")
        # Generar diccionario para almacenar valores de dbnsfp extraidos de MyVariant:
        dbnsfp_dict = {}
        # Si el dbnsfp extraído es una lista no vacía:
        if isinstance(dbnsfp, list) and len(dbnsfp) > 0:
            # Almacenar el primer elemento de la lista en el diccionario "dbnsfp_dict":
            dbnsfp_dict = dbnsfp[0]
        # Si el dbnsfp extraído es un diccionario:
        elif isinstance(dbnsfp, dict):
            # Almacenar el dbnsfp extraído en el diccionario "dbnsfp_dict":
            dbnsfp_dict = dbnsfp
        # Si el diccionario "dbnsfp_dict" no es nulo:
        if dbnsfp_dict:
            # Obtener el revel del diccionario "dbnsfp_dict":
            revel = dbnsfp_dict.get("revel")
            # Si revel es un diccionario:
            if isinstance(revel, dict):
                # Almacenar en la propiedad "REVEL" del diccionario de anotaciones el valor de la propiedad "score" del diccionario "dbnsfp_dict":
                annotations["REVEL"] = revel.get("score")
            # Si revel es un número entero, un float o una cadena de texto:
            elif isinstance(revel, (int, float, str)):
                # Almacenar en la propiedad "REVEL" del diccionario de anotaciones el valor de revel transformado a float:
                annotations["REVEL"] = float(revel)

            # Si la propiedad "PhyloP" del diccionario de anotaciones es nula, es decir, no se ha podido encontrar el PhyloP buscando en MyVariant mediante el rsid:
            if annotations["PhyloP"] is None:
                # Buscar el PhyloP en el diccionario "dbnsfp_dict" usando los isguientes campos:
                for p_key in [
                    "phylop100way_vertebrate",
                    "phyloP100way_vertebrate",
                    "phylop30way_mammalian",
                    "phyloP30way_mammalian",
                ]:
                    # Intentar obtener PhyloP a partir de cada campo definido en "p_key":
                    p_val = dbnsfp_dict.get(p_key)
                    # Si el valor encontrado no es nulo:
                    if p_val is not None:
                        # Si el valor encontrado es un diccionario:
                        if isinstance(p_val, dict):
                            # Almacenar en la propiedad "PhyloP" del diccionario de anotaciones el valor de la propiedad "score" del diccionario encontrado:
                            annotations["PhyloP"] = p_val.get("score")
                        # Si el valor encontrado es una lista no vacía:
                        elif isinstance(p_val, list) and len(p_val) > 0:
                            # Si el primer valor de la lista es un diccionario, almacenar el valor de la propiedad "score" de ese diccionario en la propiedad "Phylop" del diccionario de anotaciones. Si no es un diccionario, almacenar directamente el primer valor de la lista n la propiedad "Phylop" del diccionario de anotaciones.
                            annotations["PhyloP"] = (
                                p_val[0].get("score")
                                if isinstance(p_val[0], dict)
                                else p_val[0]
                            )
                        # Si el valor encontrado es un nñumero entero, un float o una cadena de texto:
                        elif isinstance(p_val, (int, float, str)):
                            # Almacenar en la propiedad "PhyloP" del diccionario de anotaciones el valor encontrado transformado a float:
                            annotations["PhyloP"] = float(p_val)
                        # Si el valor encontrado no es un diccionario, una lista, un número entero, un float ni una cadena de texto, pasar al siguiente p_key:
                        break

        # Extracción de gnomad_af:
        # Buscar gnomad_genome de los datos de Myvariant, o gnomad_exome en su defecto:
        gnomad = record.get("gnomad_genome") or record.get("gnomad_exome")
        # Crear un diccionario para almacenar valores de gnomad obtenidos de Myvariant:
        gnomad_dict = {}
        # Si el valor de gnomad es una lista no vacía:
        if isinstance(gnomad, list) and len(gnomad) > 0:
            # Almacenar en el diccionario "gnomad_dict" el primer elemento de la lista:
            gnomad_dict = gnomad[0]
        # Si el valor de gnomad es un diccionario:
        elif isinstance(gnomad, dict):
            # Almacenar en el diccionario "gnomad_dict" el diccionario obtenido:
            gnomad_dict = gnomad
        # Si el diccionario "gnomad_dict" no es nulo:
        if gnomad_dict:
            # Obtener valor de la propiedad "af" del diccionario:
            af = gnomad_dict.get("af")
            # Si el valor "af" es un diccionario: 
            if isinstance(af, dict):
                # Almacenar valor de la propiedad "af" del diccionario anterior en la propiedad "gnomAD_AF" del diccionario de anotaciones:
                annotations["gnomAD_AF"] = af.get("af")
            # Si el valor "af" no es un diccionario: 
            else:
                # Almacenar valor "af" en la propiedad "gnomAD_AF" del diccionario de anotaciones:
                annotations["gnomAD_AF"] = af
        # Devolver diccionario de anotaciones:
        return annotations
    def find_correct_allele(self, records, ref, alt):
        """
        Dels resultats associats a un rsID,
        busca el registre que tingui exactament
        el mateix REF i ALT que ClinVar.
        """

        if not isinstance(records, list):
            records = [records]

        ref = str(ref).upper()
        alt = str(alt).upper()

        for record in records:

            variant_id = record.get("_id")

            record_ref, record_alt = self.extract_ref_alt(
                variant_id
            )

            if (
                record_ref == ref
                and record_alt == alt
            ):
                return record

        return None
    
    def add_features_myvariant(self, df):
        # Procesar datos a partir del rsid:
        # Normalizar rsids:
        df["dbSNP"] = df["RS# (dbSNP)"].apply(self.normalize_rsid)
        # Generar el hvs id en una nueva columna llamada "hgvs_id" juntando en un string "chr" más el número de cromosoma más ":g." más la posición VCF más el alelo de referencia más ">" más el alelo alternativo.
        df["hgvs_id"] = (
            "chr"
            + df["Chromosome"].astype(str)
            + ":g."
            + df["PositionVCF"].astype(str)
            + df["ReferenceAlleleVCF"].astype(str)
            + ">"
            + df["AlternateAlleleVCF"].astype(str)
        )
        # Generar columnas CADD, gnomAD_AF, REVEL y PhyloP vacias:
        df["CADD"] = pd.NA
        df["gnomAD_AF"] = pd.NA
        df["REVEL"] = pd.NA
        df["PhyloP"] = pd.NA
        # Cargar cliente de Myvariant:
        mv = myvariant.MyVariantInfo()
        # Almacenar campos a consultar con el módulo de Myariant:
        fields = "cadd,gnomad_genome.af,gnomad_exome.af,dbnsfp"
        # Procesar datos con rsid:
        # Obtener índices de las variantes con rsid:
        indexes_with_rsid = df[df["dbSNP"].notna()].index
        rsids = (
            df.loc[indexes_with_rsid, "dbSNP"] # Obtener rsids de la columna "dbSNP" de las variantes con rsid a partir de los índices.
            .drop_duplicates() # Eliminar rsids dublicados.
            .tolist() #Añadir rsids a una lista
        )
        # Procesar rsids en lotes de 100.
        batch_size = 100
        for i in range(0, len(rsids), batch_size):
            # Recopilar rsids del lote:
            batch = rsids[i : i + batch_size]
            # Obtener datos de Myvariant para los campos especificados:
            results = mv.getvariants(batch, fields=fields)
            # Iterar para cada variante del lote:
            for record in results:
                # Recuperar rsid consultado:
                rsid = record.get("query")
                # Si no hay rsid, saltarse el resto de pasos del bucle:
                if not rsid:
                    continue
                # Buscar las filas del dataframe "df" cuyo valor en la columna "dbSNP" sea igual al rsid consultado:
                matching_rows = df[df["dbSNP"] == rsid]
                # Rellenar diccionario de anotaciones con los valores de las filas obtenidas:
                self.fill_annotations(df, matching_rows, record)

        # Procesar datos a partir del hgvs_id
        # Seleccionar índices de las variantes sin rsid:
        indexes_without_rsid = df[df["dbSNP"].isna()].index
        # Obtener hgvs_ids de las variantes sin rsid:
        hgvs_ids = (
            df.loc[indexes_without_rsid, "hgvs_id"] # Obtener valores de la columna "hgvs_id" de las variantes sin rsid a partir de sus indices. 
            .drop_duplicates() # Eliminar elementos duplicados.
            .tolist() # Añadir hgvs_ids a una lista.
        )
        # Procesar hgvs_ids en lotes de 100.
        for i in range(0, len(hgvs_ids), batch_size):
            # Procesar hgvs_ids por lotes:
            batch = hgvs_ids[i : i + batch_size]
            # Obtener datos de Myvariant para los campos especificados:
            results = mv.getvariants(batch, fields=fields)
            # Iterar para cada variante del lote:
            for record in results:
                # Recuperar hgvs_id consultado:
                hgvs_id = record.get("query")
                # Si no hay hgvs_id, saltarse el resto de pasos del bucle:
                if not hgvs_id or record.get("notfound", False):
                    continue
                # Buscar las filas del dataframe "df" cuyo valor en la columna "hgvs_id" sea igual al hgvs_id consultado:
                matching_rows = df[df["hgvs_id"] == hgvs_id]
                # Rellenar diccionario de anotaciones con los valores de las filas obtenidas:
                self.fill_annotations(df, matching_rows, record)
        # Devolver dafaframe con las columnas CADD, gnomAD_AF, REVEL y PhyloP añadidas:
        return df

    def save_df(self, input_file, csv_file, opt):
        # Leer dataframe del archivo csv:
        df = pd.read_csv(input_file, sep=",",low_memory=False)
        # Obtener CADD, gnomAD_AF, REVEL y PhyloP y añadirlos al dataframe:
        features_df = self.add_features_myvariant(df)
        # Guardar dataframe resultante en archivo csv:
        features_df.to_csv(csv_file, index=False)
        