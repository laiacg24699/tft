import pandas as pd
from cyvcf2 import VCF

class vcf:
    #################################################
    # Funciones:
    # GUARDAR ARCHIVO CSV CON DATOS DEL VCF:
    def get_vcf_df(self, vcf_file, vcf_df_file, myv, annotations_file):
        # Leer archivo vcf y almacenar su contenido en la variable "vcf"
        vcf = VCF(vcf_file)
        # Crear lista para ir almacenando las variantes del archivo vcf:
        variants = []
        # Iterar para cada variante en el contenido del archivo vcf:
        for variant in vcf:
            # Si la variante es de nucleótido simple, es decir, tiene un solo nucleótido en el alelo de referencia y en el alelo alternatio (hay que murar longitud de variant.REF, variant.ALT, y la longitud del primer elemento de variant.ALT, puesto que este campo puede ser una lista con varios elementos):
            if (len(str(variant.REF)) == 1 and len(variant.ALT) == 1 and len(variant.ALT[0]) == 1):
                # Añadir diccionario con los campos de "variant" a la lista "variants" (para el campo "variant.ALT", como puede ser una lista o un string de los cuales ya se ha comprobado que tengan una longitud de 1 elemento o 1 caracter respectivamente, se elige el primer elemento o caracter de ese campo):
                variants.append({
                    "Chromosome": variant.CHROM,
                    "RS# (dbSNP)": variant.ID,
                    "ReferenceAlleleVCF": variant.REF,
                    "AlternateAlleleVCF": variant.ALT[0],
                    "PositionVCF": variant.POS,
                    "ReferenceAllele": variant.REF,
                })
            vcf_df = pd.DataFrame(variants)
            #Añadir columna "RS# (dbSNP)" vacía (se usará para consultar datos de Myvariant si vcf contiene los "RS# (dbSNP)"): 
            # vcf_df["RS# (dbSNP)"] = pd.NA
            vcf_df.to_csv(vcf_df_file, index=False, na_rep="NA")
        myv.getAllAnnotations(vcf_df_file, annotations_file)
    #################################################
    #Función principal:
    def vcf_classifier(self, vcf_file, vcf_df_file, myv, annotations_file, ml, model_file, model, classified_file, txt_classified_file):
        self.get_vcf_df(vcf_file, vcf_df_file, myv, annotations_file)
        myv.getAllAnnotations(vcf_df_file, annotations_file)
        classified_df = pd.read_csv(annotations_file, sep=",",low_memory=False)
        variants = classified_df[["CADD","gnomAD_AF","REVEL", "PhyloP" ]].apply(pd.to_numeric,errors="coerce")
        classified_df.drop(columns=["RS# (dbSNP)"], inplace=True)
        classified_df.drop(columns=["dbSNP"], inplace=True)
        classified_df.drop(columns=["Chromosome"], inplace=True)
        classified_df.drop(columns=["ReferenceAlleleVCF"], inplace=True)
        classified_df.drop(columns=["AlternateAlleleVCF"], inplace=True)
        classified_df.drop(columns=["PositionVCF"], inplace=True)
        
        classified_df["ClinicalSignificance"] = ml.classify_variants(variants, model_file, model)
        classified_df.to_csv(classified_file , sep="\t", index=None) #Archivo csv
        classified_df.to_string(txt_classified_file, index=False) #Archivo txt
