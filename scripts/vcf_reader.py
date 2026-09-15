import pandas as pd
from cyvcf2 import VCF

class vcf:
    def read(self, vcf_file):
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
        # Devolver dataframe creado a partir de la lista de diccionarios "variants".
        return pd.DataFrame(variants)