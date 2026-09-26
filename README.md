
# Introducción

Este repositorio contiene el código desarrollado para el Trabajo de Fin de Máster (TFM) del Máster en Bioinformática de la Universidad Internacional de Valencia.

El código desarrollado en Python permite crear modelos de Machine Learning para la clasificación de variantes genéticas según las cinco categorías establecidas por las recomendaciones del American College of Medical Genetics and Genomics (ACMG): Pathogenic, Likely pathogenic, Uncertain significance, Likely benign y Benign.

El proyecto utiliza un conjunto de datos de variantes genéticas procedentes de ClinVar, correspondientes a agosto de 2026. El conjunto de datos se filtra mediante el script filter.py y, posteriormente, se incorporan características adicionales mediante los scripts getEnsemblData.py y getMyvData.py. La información obtenida se utiliza como entrada para el entrenamiento de los modelos de Machine Learning, proceso que se ejecuta mediante el script MLModelCreator.py.

En este proyecto se han desarrollado dos modelos de clasificación basados en diferentes algoritmos de Machine Learning: Random Forest y Gradient Boosting.

Una vez entrenados, los modelos pueden utilizarse para clasificar variantes genéticas a partir de un archivo en formato VCF mediante el script VCFClassifier.py.

# Archivos utilizados

El repositorio se organiza en diferentes carpetas según la función de los archivos dentro del proyecto:

.
├── MLModelCreator.py
├── README.md
├── VCFClassifier.py
├── data/
│   ├── raw/
│   │   ├── variant_summary.txt.gz
│   │   └── Pfeiffer.vcf
│   └── processed/
│       ├── filtered.csv
│       ├── features_ensembl.csv
│       ├── features_myv.csv
│       ├── Pfeiffer_annotated.csv
│       └── Pfeiffer_read.csv
├── models/
│   ├── acmg_classifier_gb.joblib
│   └── acmg_classifier_rf.joblib
├── results/
│   ├── gb_classified.csv
│   ├── gb_classified.txt
│   ├── gb_val.txt
│   ├── rf_classified.csv
│   ├── rf_classified.txt
│   └── rf_val.txt
├── scripts/
│   ├── filter.py
│   ├── getEnsemblData.py
│   ├── getMyvData.py
│   ├── machine_learning.py
│   └── vcf_reader.py
└── requirements.txt

## Datos

La carpeta data/ contiene los archivos utilizados durante las diferentes etapas del procesamiento de las variantes.

raw/ contiene los datos de entrada utilizados en el proyecto. Incluye el archivo de variantes de ClinVar (variant_summary.txt.gz) y el archivo VCF utilizado para las prueba de clasificación,  obtenido de https://github.com/davetang/learning_vcf_file/tree/main/eg. 

processed/ contiene los archivos generados durante el procesamiento de los datos. Estos incluyen el conjunto de variantes filtrado (filtered.csv), los archivos con las características obtenidas mediante Ensembl y MyVariant.info (features_ensembl.csv y features_myv.csv), y los archivos intermedios generados durante la lectura y anotación de los archivos VCF.

## Modelos

La carpeta models/ contiene los modelos de Machine Learning entrenados y almacenados en formato .joblib:

acmg_classifier_rf.joblib: modelo basado en Random Forest.
acmg_classifier_gb.joblib: modelo basado en Gradient Boosting.
Resultados

## Results

La carpeta results/ contiene los resultados obtenidos durante la validación y clasificación de variantes.

Los archivos rf_* corresponden a los resultados obtenidos mediante el modelo Random Forest, mientras que los archivos gb_* corresponden a los obtenidos mediante Gradient Boosting. Los archivos *_val.txt contienen los resultados de la validación de los modelos.

## Scripts

La carpeta scripts/ contiene los módulos Python utilizados para las diferentes etapas del proyecto:

- filter.py: realiza el filtrado del conjunto de variantes de ClinVar.
- getEnsemblData.py: obtiene información adicional de las variantes mediante Ensembl.
- getMyvData.py: obtiene información adicional mediante MyVariant.info.
- machine_learning.py: contiene las funciones relacionadas con la creación y entrenamiento de los modelos de Machine Learning.
- vcf_reader.py: contiene las funciones utilizadas para la lectura y procesamiento de archivos VCF.

## Scripts principales
- MLModelCreator.py: ejecuta el proceso de creación y entrenamiento de los modelos de clasificación.
- VCFClassifier.py: utiliza los modelos entrenados para clasificar las variantes contenidas en un archivo VCF.

## Otros archivos
- requirements.txt: contiene las dependencias de Python necesarias para ejecutar el proyecto.
- README.md: contiene la documentación y las instrucciones de uso del repositorio.

