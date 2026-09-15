import pandas as pd
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
    # Columnas con las variables de entrada del modelo:
    feature_columns = [
        "CADD",
        "gnomAD_AF",
        "REVEL",
        "PhyloP"
    ]
    def create_model(self, machinelearning_df, classifier):
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
        # Comparar valores objetios predichos con los valores objetivo de los datos de prueba para obtener la precisión global del modelo:
        accuracy = accuracy_score(y_test, y_pred)
        print("Accuracy:", accuracy)
        # Generar las métricas:
        # Presición: para cada clase, el porcentaje de elementos predecidos en los resultados que realmente pertenecian a esa clase.
        # Recall: para todas las variantes de cada clase de los datos de entrada, el porcentaje de cuantas ha identificado correctamente.
        # F1-score: F1=2×precision×recall​/(precision+recall)
        # Support: para cada clase, el número real de variantes que pertenecían a esa clase. 
        print("\nClassification report:")
        print(classification_report(y_test, y_pred))
        # Generar matriz de confusión (indica clases que confunde el modelo): 
        print("\nConfusion matrix:")
        print(confusion_matrix(y_test, y_pred))

        #Report obtenido:
        # accuracy: 0.5606

        # Classification report:
        #                         precision    recall  f1-score   support

        #                 Benign       0.91      0.83      0.87      2000
        #         Likely benign       0.77      0.47      0.59      2000
        #     Likely pathogenic       0.32      0.77      0.46      2000
        #             Pathogenic       0.67      0.38      0.49      2000
        # Uncertain significance       0.64      0.35      0.45      2000

        #             accuracy                           0.56     10000
        #             macro avg       0.66      0.56      0.57     10000
        #         weighted avg       0.66      0.56      0.57     10000


        # Confusion matrix:
        # [[1668  135  114    4   79]
        # [ 135  948  803    7  107]
        # [   7   54 1531  295  113]
        # [   2   29 1117  764   88]
        # [  17   65 1147   76  695]]


        # Devolver modelo entrenado:
        return model

    def save_model(self, input_file, classifier, model_file):
        # Leer dataframe de un archio csv:
        df = pd.read_csv(input_file, sep=",",low_memory=False)
        # Entrenar y generar el modelo:
        model = self.create_model(df, classifier)
        # Guardar modelo como archio joblib.
        joblib.dump(
            model,
            model_file
        )

    def classify_variants(self, new_features, classifier):
        # Si se selecciona el clasificador "random_forest", se carga el modelo de random forest generado anteriormente:
        if classifier == 'random_forest':
            model = joblib.load(self.rf_model_file)
        # Si no se selecciona el clasificador "random_forest", se carga el modelo de gradient boosting generado anteriormente:
        else:
            model = joblib.load(self.gb_model_file)
        # Usar modelo seleccionado para obtener el "ClinicalSignificance" de cada variante de la lista "new_features". Cada elemento de esta lista corresponde a una variante, y es una lista que contiene cada uno de los valores de las columnas de entrada del modelo:
        return model.predict(new_features)