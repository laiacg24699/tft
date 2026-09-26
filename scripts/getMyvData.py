import pandas as pd
import myvariant
import re
import ast

class MyvData:
    ################################################
    # Funciones MLModelCreator
    # NORMALIZAR RSID
    def normalize_rsid(self, value):
        if pd.isna(value) or value is None:
            return None
        value = str(value).strip()
        if value in ["", "-", "na", ".", "NA", "-1", "nan"]:
            return None
        if value.startswith("rs"):
            return value
        try:
            return f"rs{int(float(value))}"
        except (ValueError, TypeError):
            return None
    # EXTRAER PHYLOP DEL REGISTRO DE MYVARIANT
    def extract_phylop_from_record(self, record):
        if not record or record.get("notfound", False):
            return None
        phylop_val = None
        # Buscar PhyloP en dbNSFP
        dbnsfp = record.get("dbnsfp")
        dbnsfp_dict = dbnsfp[0] if isinstance(dbnsfp, list) and len(dbnsfp) > 0 else (dbnsfp if isinstance(dbnsfp, dict) else {})
        
        if dbnsfp_dict:
            for p_key in [
                "phylop100way_vertebrate",
                "phyloP100way_vertebrate",
                "phylop30way_mammalian",
                "phyloP30way_mammalian",
                "phylop17way_primate",
                "phyloP17way_primate"
            ]:
                p_val = dbnsfp_dict.get(p_key)
                if p_val is not None:
                    if isinstance(p_val, dict):
                        phylop_val = p_val.get("score")
                    elif isinstance(p_val, list) and len(p_val) > 0:
                        phylop_val = p_val[0].get("score") if isinstance(p_val[0], dict) else p_val[0]
                    elif isinstance(p_val, (int, float, str)):
                        phylop_val = p_val
                    
                    if phylop_val is not None:
                        break

        # Si no se encuentra en dbNSFP, buscar en CADD
        if phylop_val is None:
            cadd = record.get("cadd")
            cadd_dict = cadd[0] if isinstance(cadd, list) and len(cadd) > 0 else (cadd if isinstance(cadd, dict) else {})
            if cadd_dict:
                cadd_phylop = cadd_dict.get("phylop", {})
                if isinstance(cadd_phylop, dict):
                    phylop_val = cadd_phylop.get("vertebrate") or cadd_phylop.get("mammalian") or cadd_phylop.get("primate")

        # Intentar convertir a float
        if phylop_val is not None:
            try:
                return float(phylop_val)
            except (ValueError, TypeError):
                return None

        return None
    # OBTENER VALOR DE PHYLOP POR RSID
    def get_Phylop_rsid(self, df, batch_size, mv, fields):
        indexes_with_rsid = df[df["dbSNP_norm"].notna()].index
        rsids = df.loc[indexes_with_rsid, "dbSNP_norm"].drop_duplicates().tolist()
        for i in range(0, len(rsids), batch_size):
            batch = rsids[i : i + batch_size]
            try:
                results = mv.getvariants(batch, fields=fields)
                for record in results:
                    query_rsid = record.get("query")
                    if not query_rsid:
                        continue
                    score = self.extract_phylop_from_record(record)
                    if score is not None:
                        matching_rows = df[df["dbSNP_norm"] == query_rsid].index
                        df.loc[matching_rows, "PhyloP"] = score
            except Exception as e:
                continue
        return df
    # OBTENER VALOR DE PHYLOP POR HGVS
    def get_Phylop_hgvs(self, df, batch_size, mv, fields):
        pending_indexes = df[df["PhyloP"].isna()].index
        hgvs_ids = df.loc[pending_indexes, "hgvs_id"].drop_duplicates().tolist()
        for i in range(0, len(hgvs_ids), batch_size):
            batch = hgvs_ids[i : i + batch_size]
            try:
                results = mv.getvariants(batch, fields=fields)
                for record in results:
                    query_hgvs = record.get("query")
                    if not query_hgvs or record.get("notfound", False):
                        continue
                    score = self.extract_phylop_from_record(record)
                    if score is not None:
                        matching_rows = df[df["hgvs_id"] == query_hgvs].index
                        df.loc[matching_rows, "PhyloP"] = score
            except Exception as e:
                continue
        return df
    # OBTENER VALOR DE PHYLOP
    def get_Phylop(self, input_file, batch_size=100):
        # Cargar datos del csv
        df = pd.read_csv(input_file, low_memory=False)
        # Crear columnas auxiliares para las consultas
        if "RS# (dbSNP)" in df.columns:
            df["dbSNP_norm"] = df["RS# (dbSNP)"].apply(self.normalize_rsid)
        else:
            df["dbSNP_norm"] = None
        df["hgvs_id"] = (
            "chr"
            + df["Chromosome"].astype(str).str.replace("chr", "").str.strip()
            + ":g."
            + df["PositionVCF"].astype(str).str.strip()
            + df["ReferenceAlleleVCF"].astype(str).str.strip().str.upper()
            + ">"
            + df["AlternateAlleleVCF"].astype(str).str.strip().str.upper()
        )
        df["PhyloP"] = pd.NA
        mv = myvariant.MyVariantInfo()
        fields = "cadd.phylop,dbnsfp"
        # Procesar por lotes primero las variantes que tienen rsID
        df = self.get_Phylop_rsid(df, batch_size, mv, fields)
        # Processar por HGVS las variantes que NO tienen rsID o que han quedado sin PhyloP
        df = self.get_Phylop_hgvs(df, batch_size, mv, fields)
        # Limpiar columnas temporales
        df.drop(columns=["dbSNP_norm", "hgvs_id"], inplace=True, errors="ignore")
        # Devolver dataframe resultante
        return df
    # EXTRAER REF Y ALT DEL IDENTIFICADOR HGVS
    def extract_ref_alt(self, variant_id):
        # Si el variant id no es un string, devolver valores nulos.
        if not isinstance(variant_id, str):
            return None, None
        # Buscar estructura :g.POSICIONREF>ALT
        match = re.search(
            r":g\.\d+([ACGT]+)>([ACGT]+)",
            variant_id
        )
        # Si se encuentra la estructura, devolver REF y ALT.
        if match:
            return match.group(1), match.group(2)
        return None, None
    # EXTRAER EL VALOR DE REVEL DE LA CONSULTA A MYVARIANT
    def extract_revel(self, record, ref, alt):
        # Extraer variant id.
        variant_id = record.get("_id", "")
        # Extraer REF y ALT del variant id.
        record_ref, record_alt = self.extract_ref_alt(variant_id)
        # Comprobar los alelos cuando están disponibles.
        if ref and alt:
            if (
                record_ref is not None
                and record_alt is not None
                and (
                    record_ref != str(ref).upper()
                    or record_alt != str(alt).upper()
                )
            ):
                return None
        dbnsfp = record.get("dbnsfp")
        # MyVariant puede devolver dbnsfp como lista.
        if isinstance(dbnsfp, list) and len(dbnsfp) > 0:
            dbnsfp_dict = dbnsfp[0]
        # O como diccionario.
        elif isinstance(dbnsfp, dict):
            dbnsfp_dict = dbnsfp
        else:
            dbnsfp_dict = {}
        if not dbnsfp_dict:
            return None
        revel = dbnsfp_dict.get("revel")
        if isinstance(revel, dict):
            return revel.get("score")
        elif isinstance(revel, (list, tuple, str, int, float)):
            return revel
        return None
    # PROCESAR REVEL
    def arrange_revel(self, revel):
        # Si REVEL es una lista o una tupla.
        if isinstance(revel, (list, tuple)):
            valid_floats = []
            for x in revel:
                if x is None:
                    continue
                try:
                    value = float(x)
                    if pd.notna(value):
                        valid_floats.append(value)
                except (ValueError, TypeError):
                    continue
            return max(valid_floats) if valid_floats else None
        # Si REVEL es una cadena.
        elif isinstance(revel, str):
            try:
                parsed = ast.literal_eval(revel.strip())
                if isinstance(parsed, (list, tuple, set)):
                    valid_floats = []
                    for x in parsed:
                        if x is None:
                            continue
                        try:
                            value = float(x)
                            if pd.notna(value):
                                valid_floats.append(value)
                        except (ValueError, TypeError):
                            continue
                    return max(valid_floats) if valid_floats else None
                else:
                    return float(parsed)
            except (ValueError, SyntaxError, TypeError):
                return None
        # Si REVEL es un número.
        elif isinstance(revel, (float, int)):
            if pd.notna(revel):
                return float(revel)
        return None
    # RELLENAR REVEL EN LA FILA ADECUADA DEL DATAFRAME RESULTANTE
    def fill_revel(self, df, matching_rows, record):
        for idx in matching_rows.index:
            # Recuperar alelo de referencia.
            ref = df.at[idx, "ReferenceAlleleVCF"]
            # Recuperar alelo alternativo.
            alt = df.at[idx, "AlternateAlleleVCF"]
            # Obtener anotaciones.
            revel = self.extract_revel(
                record,
                ref,
                alt
            )
            # Rellenar REVEL.
            if revel is not None:
                revel_value = self.arrange_revel(revel)
                if revel_value is not None:
                    df.at[idx, "REVEL"] = revel_value
    # OBTENER DATAFRAME CON EL REVEL RELLENADO
    def get_Revel(self, df):
        # Preparar identificadores
        df["dbSNP"] = df["RS# (dbSNP)"].apply(self.normalize_rsid)
        df["hgvs_id"] = (
            "chr"
            + df["Chromosome"].astype(str)
            + ":g."
            + df["PositionVCF"].astype(str)
            + df["ReferenceAlleleVCF"].astype(str)
            + ">"
            + df["AlternateAlleleVCF"].astype(str)
        )
        # Crear columnas
        df["REVEL"] = pd.NA
        # Cliente MyVariant
        mv = myvariant.MyVariantInfo()
        fields = "dbnsfp"
        batch_size = 100
        # Procesar rsids
        indexes_with_rsid = df[df["dbSNP"].notna()].index
        rsids = (df.loc[indexes_with_rsid,"dbSNP"].drop_duplicates().tolist())
        for i in range(0, len(rsids), batch_size):
            batch = rsids[i:i + batch_size]
            results = mv.getvariants(batch,fields=fields)
            for record in results:
                rsid = record.get("query")
                if not rsid:
                    continue
                matching_rows = df[df["dbSNP"] == rsid]
                self.fill_revel(df,matching_rows,record)
        # Procesar HGVS
        indexes_without_rsid = df[df["dbSNP"].isna()].index
        hgvs_ids = (df.loc[indexes_without_rsid,"hgvs_id"].drop_duplicates().tolist())
        for i in range(0, len(hgvs_ids), batch_size):
            batch = hgvs_ids[i:i + batch_size]
            results = mv.getvariants(batch,fields=fields)
            for record in results:
                hgvs_id = record.get("query")
                if (not hgvs_id or record.get("notfound", False)):
                    continue
                matching_rows = df[df["hgvs_id"] == hgvs_id]
                self.fill_revel(df,matching_rows,record)
        return df
    ################################################
    # Funcion principal
    def getRevelPhylop(self, features_vep_file, features_myv_file):
        df = self.get_Phylop(features_vep_file, batch_size=100)
        df = self.get_Revel(df)
        df.to_csv(features_myv_file, index=False, na_rep="NA")
    ################################################
    # Funciones vcf_reader
    # OBTENER DATOS DE LAS ANOTACIONES EN LOS DATOS DE CONSULTA DE MYVARIANT
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
    # RELLENAR DATAFRAME CON LOS DATOS DE LAS 4 ANOTACIONES
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
    # OBTENER INFORMACIÓN DE LAS 4 ANOTACIONES DESDE MYVARIANT
    def getAllInfoMyVariant(self, df):
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
    ################################################
    # Funcion principal
    def getAllAnnotations(self, input_file, csv_file):
        # Leer dataframe del archivo csv:
        df = pd.read_csv(input_file, sep=",",low_memory=False)
        # Obtener CADD, gnomAD_AF, REVEL y PhyloP y añadirlos al dataframe:
        features_df = self.getAllInfoMyVariant(df)
        # Guardar dataframe resultante en archivo csv:
        features_df.to_csv(csv_file, index=False)

