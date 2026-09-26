import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer

import joblib
from sklearn.metrics import (
    accuracy_score,
    classification_report,
    confusion_matrix
)

class machine_learning:
    # CONFIGURACIÓN
    # Columnas con las variables de entrada del modelo:
    feature_columns = [
        "CADD",
        "gnomAD_AF",
        "REVEL",
        "PhyloP"
    ]
    #################################################
    # Funciones:
    # GUARDAR DATOS DEL MODELO
    def save_model_info(self, accuracy, x_report, c_matrix, model_info_file):
        # Limpiar contenido del archivo txt:
        open(model_info_file, 'w').close()
        #Ir guardando nueva información al archio txt:
        with open(model_info_file, 'a') as f:
            f.write('Classification report:\n')
            f.write(x_report)
            f.write('\nConfusion matrix:\n')
            f.write(np.array_str(c_matrix))
    # CREAR MODELO DE MACHINE LEARNING
    def create_model(self, machinelearning_df, classifier, model_info_file):
        # Seleccionar variables de entrada que se usaran para entrenar el modelo:
        X = machinelearning_df[self.feature_columns]
        # Seleccionar ariable objetivo que deberá predecir el modelo:
        y = machinelearning_df["ClinicalSignificance"]
        # Dividir aleatoriamente variantes de la muestra en variantes de entrenamiento (80%), que se usaran para entrenar el modelo, y variantes de prueba (20%), que se usaran para probar el modelo:
        # Para la diisión aleatoria se usa una semilla random_state con un alor de 42. Esto permitirá obtener siempre las mismas muestras aleatorias en cada grupo.
        X_train, X_test, y_train, y_test = train_test_split(
            X,
            y,
            test_size=0.2,
            random_state=42,
            stratify=y
        )
        # Si se selecciona el clasificador "random_forest", se genera un Random Forest con 100 estimadores:
        if classifier == 'random_forest':
            model = RandomForestClassifier(
                n_estimators=100,
                random_state=42
            )
        # Si no se selecciona el classificador "random_forest", entendemos que se quiere usar "gradient_boosting".
        # Se crea un Gradient Boosting con pipeline que rellenará primero los valores con valor NaN usando la mediana del resto de valores de la columna correspondiente (CAAD, gnomAD_AF, etc.) y luego generará el Gradient Boosting con 100 estimadores.
        else:
            model = Pipeline([
                ("imputer", SimpleImputer(strategy="median")),
                ("classifier", GradientBoostingClassifier(
                    n_estimators=100,
                    random_state=42
                ))
            ])
        # Entrenar el modelo con los datos de entrenamiento:
        model.fit(X_train, y_train)
        # Predecir valores objetivo a partir de las columnas de variable de entrada de los datos de prueba:
        y_pred = model.predict(X_test)
        self.save_model_info(accuracy_score(y_test, y_pred), classification_report(y_test, y_pred, digits=4), confusion_matrix(y_test, y_pred), model_info_file)
        # Generar las métricas:
        # Presición: para cada clase, el porcentaje de elementos predecidos en los resultados que realmente pertenecian a esa clase.
        # Recall: para todas las variantes de cada clase de los datos de entrada, el porcentaje de cuantas ha identificado correctamente.
        # F1-score: F1=2×precision×recall​/(precision+recall)
        # Devolver modelo entrenado:
        return model
    #################################################
    # Funciones principales
    def save_model(self, input_file, classifier, model_file, model_info_file):
        # Leer dataframe de un archio csv:
        df = pd.read_csv(input_file, sep=",",low_memory=False)
        # Entrenar y generar el modelo:
        model = self.create_model(df, classifier, model_info_file)
        # Guardar modelo como archio joblib.
        joblib.dump(
            model,
            model_file
        )

    def classify_variants(self, new_features, model_file, classifier):
        # Si se selecciona el clasificador "random_forest", se carga el modelo de random forest generado anteriormente:
        if classifier == 'random_forest':
            model = joblib.load(model_file)
        # Si no se selecciona el clasificador "random_forest", se carga el modelo de gradient boosting generado anteriormente:
        else:
            model = joblib.load(model_file)

        # Usar modelo seleccionado para obtener el "ClinicalSignificance" de cada variante de la lista "new_features". Cada elemento de esta lista corresponde a una variante, y es una lista que contiene cada uno de los valores de las columnas de entrada del modelo:
        return model.predict(new_features)