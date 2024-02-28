def predict_function(X, file_name="C:/Users/nakamukai/Desktop/RO_prediction_of_PNPA/"): 
    import shap
    import numpy as np
    import pandas as pd
    import joblib

    import os 
    os.chdir("{}/src".format(file_name))
    
    loaded_model = joblib.load("{}/model/model_mix_count_PredRet_Nematostella_pca_200.pkl".format(file_name))

    return loaded_model.predict(X, X)


def shap_calculation(df_pair):
    
    import numpy as np
    import shap
    import pandas as pd

    # Create "explainer" to calculate SHAP values. 
    explainer = shap.Explainer(predict_function, df_pair)
    
    shap_values = explainer.shap_values(df_pair)
    
    shap_values_class0 = shap_values[0] 
    shap_values_class1 = shap_values[1] 


    shap_values_mean_class0 = np.mean(shap_values_class0, axis=1)
    shap_values_mean_class1 = np.mean(shap_values_class1, axis=1)

    shap_1 = np.array([shap_values_mean_class0, shap_values_mean_class1])
    
    return(shap_1)


import shap
import numpy as np
import pandas as pd
import joblib

file_name="C:/Users/nakamukai/Desktop/RO_prediction_of_PNPA/"

df_test = pd.read_csv("{}Fingerprints/PredRet_Nematostella/PredRet_Nematostella_pca_200_test.csv".format(file_name))

df_count = df_test.iloc[:, 3:]

list_1 = [n for n in range(0, 22)]
list_2 = [n for n in range(22, 31)]

shap_n = []
k = 0
for i in list_1:
    for j in list_2:
        df_pair = df_count.iloc[[i, j], :]
        df_pair.values
        shap_1 = shap_calculation(df_pair)
        shap_n.append(shap_1)
        k += 1
        print(k)

import matplotlib.pyplot as plt

shap_values_2 = np.concatenate(shap_n, axis=0)

k = 0
df_pairs = []
for i in list_1:
    for j in list_2:
        df_pair = df_count.iloc[[i, j], :]
        df_pair.values
        df_pairs.append(df_pair)
        df_pairs_2 = np.concatenate(df_pairs, axis=0)
df_pairs_2.shape

feature_names = [f"PC{i+1}" for i in range(200)]
df_pairs_2 = pd.DataFrame(df_pairs_2, columns=feature_names)

print(shap.summary_plot(shap_values_2, df_pairs_2, plot_type="bar"))
print(shap.summary_plot(shap_values_2, df_pairs_2))