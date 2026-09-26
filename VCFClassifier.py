from scripts.getMyvData import MyvData
from scripts.machine_learning import machine_learning
from scripts.vcf_reader import vcf

#########################################################
# Clasificación con los modelos de machine learning:
# Nombre del archivo vcf:
vcf_file = "data/raw/Pfeiffer.vcf"
vcf_df_file = "data/processed/Pfeiffer_read.csv"
annotations_file = "data/processed/Pfeiffer_annotated.csv"
# Nombre de los archivos de los modelos:
rf_model_file = "models/acmg_classifier_rf.joblib"
gb_model_file = "models/acmg_classifier_gb.joblib"
# Nombre del archivo csv resultante:
rf_classified_file = "results/rf_classified_Pfeiffer.csv"
gb_classified_file = "results/gb_classified_Pfeiffer.csv"
#Nombre del archivo txt resultante:
rf_txt_classified_file = "results/rf_classified_Pfeiffer.txt"
gb_txt_classified_file = "results/gb_classified_Pfeiffer.txt"
#########################################################
# Creando objetos a partir de las clases de los scripts:
myv = MyvData()
ml = machine_learning()
vcf_reader = vcf()
#CLASIFICAR DATOS:
vcf_reader.vcf_classifier(vcf_file, vcf_df_file, myv, annotations_file, ml, rf_model_file, 'random_forest', rf_classified_file, rf_txt_classified_file)
vcf_reader.vcf_classifier(vcf_file, vcf_df_file, myv, annotations_file, ml, gb_model_file, 'gradient_boosting', gb_classified_file, gb_txt_classified_file)