import pandas as pd
import numpy as np
import joblib
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split


TARGET_SYSTEMS = [
    "Eawag_XBridgeC18", "FEM_long","RIKEN", "UFZ_Phenomenex",
    "Nematostella_1", "Nematostella_2", "Nematostella_3", "Nematostella_4", "Nematostella_5","Reference_Compounds"
]
MAX_PAIRS_PER_SYSTEM = None  


rts = pd.read_csv("C:/Users/nakamukai/Desktop/RO_prediction_of_PNPA/rts.csv")
fps = pd.read_csv("C:/Users/nakamukai/Desktop/RO_prediction_of_PNPA/Fingerprints/PredRet_Nematostella/PredRet_Nematostella_pca_200.csv")


rts = rts[rts['system'].isin(TARGET_SYSTEMS)]
data = pd.merge(rts, fps, on='inchi', how='inner')


fp_cols = [col for col in data.columns if col.startswith("PC")]
def generate_pairwise_data(df, feature_cols, max_pairs=None):
    pairs = []
    labels = []

    n = len(df)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            rt_i = df.iloc[i]["rt"]
            rt_j = df.iloc[j]["rt"]
            if rt_i == rt_j:
                continue
            vec_i = df.iloc[i][feature_cols].values
            vec_j = df.iloc[j][feature_cols].values

            combined = np.concatenate([vec_i, vec_j]).astype(np.float32)
            label = 1 if rt_i > rt_j else 0
            pairs.append(combined)
            labels.append(label)

    if max_pairs is not None and len(labels) > max_pairs:
        idx = np.random.choice(len(labels), max_pairs, replace=False)
        return np.array(pairs, dtype=np.float32)[idx], np.array(labels)
    return np.array(pairs, dtype=np.float32), np.array(labels)

X_all, y_all = generate_pairwise_data(data, fp_cols, max_pairs=None)

final_clf = RandomForestClassifier(n_estimators=100, random_state=42)
final_clf.fit(X_all, y_all)

joblib.dump(final_clf, "C:/Users/nakamukai/Desktop/CCS_Prediction/rf_model_PredRet_Nematostella.pkl")
