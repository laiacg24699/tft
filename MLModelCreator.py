from scripts.filter import filter
from scripts.getEnsemblData import EnsemblData
from scripts.getMyvData import MyvData
from scripts.machine_learning import machine_learning

#########################################################
# Creación de los modelos de machine learning:
# Nombre del archivo a leer:
input_file = "data/raw/variant_summary.txt.gz"
# Nombre archivos csv intermedios:
filter_file = "data/processed/filtered.csv"
features_ensembl_file = "data/processed/features_ensembl.csv"
features_myv_file = "data/processed/features_myv.csv"
# Nombre del archivo del modelo:
rf_model_file = "models/acmg_classifier_rf.joblib"
gb_model_file = "models/acmg_classifier_gb.joblib"
# Archivos de salida de la validación de los modelos:
rf_model_info_file = "results/rf_val.txt"
gb_model_info_file = "results/gb_val.txt"
#########################################################
# Creando objetos a partir de las clases de los scripts:
ft = filter()
ensembl = EnsemblData()
myv = MyvData()
ml = machine_learning()
########################################################
# GENERAR MODELOS:
# Generar archivo filter.csv a partir de archivo .txt.gz
ft.save_df(input_file, filter_file)
# Generar archivo features_ensembl.csv añadiendo cadd y gnomead_af consultados en Ensembl REST API al archivo filtered.csv
ensembl.getCaddGnomead(filter_file, features_ensembl_file)
# Generar archivo features_ensembl.csv añadiendo revel y phylop consultados en MyVariant al archivo features_vep.csv
myv.getRevelPhylop(features_ensembl_file, features_myv_file)
# Generar archivos de modelo:
ml.save_model(features_myv_file, 'random_forest', rf_model_file, rf_model_info_file)
ml.save_model(features_myv_file, 'gradient_boosting', gb_model_file, gb_model_info_file)
