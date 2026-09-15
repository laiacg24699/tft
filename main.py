import pandas as pd
from scripts.filter import filter
from scripts.features import features
from scripts.machine_learning import machine_learning
from scripts.vcf_reader import vcf

#########################################################
# Nombres de los archios:
#  Nombre del archivo a leer
input_file = "data/raw/variant_summary.txt.gz"
# Nombre del archivo csv resultante
filter_file = "data/processed/filtered.csv"
features_file = "data/processed/features.csv"
# Nombre del archivo del modelo
rf_model_file = "models/acmg_classifier_rf.joblib"
gb_model_file = "models/acmg_classifier_gb.joblib"
#  Nombre del archivo vcf
# vcf_file = "data/raw/chr7_C150.vcf"
vcf_file = "data/raw/variants_chr1.vcf"
#  Nombre del archivo a leer
vcf_input_file = "data/raw/vcf_variant_summary.txt.gz"
#########################################################
#Para el vcf:
# Nombre del archivo csv resultante
rf_classified_file = "data/processed/rf_classified.csv"
gb_classified_file = "data/processed/gb_classified.csv"
#Nombre del archivo txt resultant
rf_txt_classified_file = "data/processed/rf_classified.txt"
gb_txt_classified_file = "data/processed/gb_classified.txt"

# Creando objetos a partir de las clases de los scripts:
ft = filter()
feat = features()
ml = machine_learning()
vcf_reader = vcf()
########################################################
# #GENERAR MODELOS:
# # Generar archio filterin.csv a partir de archio .txt.gz
#  ft.save_df(input_file, filter_file)
# Generar archivo features.csv añadiendo datos consultados en MyVariant a filtering.csv:
# feat.save_df(filter_file, features_file, 'myvariant')
# #Generar archivos de modelo:
# ml.save_model(features_file, 'random_forest', rf_model_file)
# ml.save_model(features_file, 'gradient_boosting', gb_model_file)
########################################################
#CLASIFICAR DATOS:
#Leer archio vcf:
vcf_df = vcf_reader.read(vcf_file)
#Añadir columna "RS# (dbSNP)" vacía (se usará para consultar datos de Myvariant si vcf contiene los "RS# (dbSNP)"): 
vcf_df["RS# (dbSNP)"] = pd.NA
#Añadir columnas de datos a vcf_df:
classified_df = feat.add_features_myvariant(vcf_df)
#Seleccionar columnas del dataframe releantes para hacer la clasificación:
variants = classified_df[["CADD", "REVEL", "PhyloP", "gnomAD_AF"]].apply(pd.to_numeric,errors="coerce")
#Clasificar variables con el modelo indicado:
classified_df["ClinicalSignificance"] = ml.classify_variants(variants.values.tolist(), 'random_forest')
# classified_df["ClinicalSignificance"] = ml.classify_variants(variants.values.tolist(), 'gradient_boosting')
#Eliminar columnas inecesarias para el archio de resultados:
classified_df.drop(columns=["RS# (dbSNP)"], inplace=True)
classified_df.drop(columns=["dbSNP"], inplace=True)
classified_df.drop(columns=["Chromosome"], inplace=True)
classified_df.drop(columns=["ReferenceAlleleVCF"], inplace=True)
classified_df.drop(columns=["AlternateAlleleVCF"], inplace=True)
classified_df.drop(columns=["PositionVCF"], inplace=True)
#Guardar resultados:
#classified_df.to_csv(rf_classified_file , sep="\t", index=None) #Archio csv
classified_df.to_string(rf_txt_classified_file, index=False) #Archivo txt