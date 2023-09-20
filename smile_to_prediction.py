# Function to convert SMILES to MACCS fingerprints
def smile_to_maccs(smile, file_name="C:/Users/nakamukai/Desktop/"):
  
  import pyper

  r = pyper.R(RCMD= "C:/Program Files/R/R-4.2.1/bin/R.exe")

  r("library(parse)")
  r("library(rcdk)")
  r.assign("smile", smile)
  
  # Load a file containing MACCS fingerprint descriptors
  maccs_list_file = "{}maccs_list_166.csv".format(file_name)
  r.assign("maccs_list_file", maccs_list_file)
  
  r('''
  # Create a vector of descriptors
  DF <- read.csv(maccs_list_file)
  DF_vector <- c(DF$maccs)
  
  # Convert to MACCS fingerprints
  mol <- parse.smiles(smile)[[1]]
  fps_count_1 <- get.fingerprint(mol, type='substructure', fp.mode='count', substructure.pattern=DF_vector)
  Nums <- c()
  for (j in 1:166) {
    x <- floor(fps_count_1@features[[j]]@count)
    Nums <- c(Nums, c(x))
  }

  DF_finger <- data.frame(x = Nums)
  DF_finger.t <- t(DF_finger)
  
  # Rename the columns of the data frame
  new_colnames <- paste("maccs_", 1:166, sep="") 
  
  ''')
  
  fps = r.get("DF_finger.t")
  return(fps)


# A function to convert SMILES to bit data of counting fingerprints
def smile_to_fps(smile, type_1, file_name="C:/Users/nakamukai/Desktop/"):
  
  import pyper

  r = pyper.R(RCMD= "C:/Program Files/R/R-4.2.1/bin/R.exe")

  r("library(parse)")
  r("library(rcdk)")
  r.assign("smile", smile)
  r.assign("type_1", type_1)
  r.assign("file_name", file_name)
  
  r('''
  # Convert to fingerprints for each type
  mol <- parse.smiles(smile)[[1]]
  fps_count <- get.fingerprint(mol, type= type_1, fp.mode='count')
  Nums <- c()
  columns <- c()
  
  # Extract bit data
  for (j in 1:length(fps_count@features)) {
  x <- floor(fps_count@features[[j]]@count)
  feature <- fps_count@features[[j]]@feature
  Nums <- c(Nums, c(x))
  columns <- c(columns, c(feature))
  }
  
  # Create a data frame with the existing vector Nums
  df <- data.frame(matrix(Nums, nrow = 1))

  # Set column names
  colnames(df) <- c(columns)
  
  ''')
  fps = r.get("df")
  return(fps)


# A function to convert MACCS counting fingerprint's bit data derived from R into a bit sequence
def maccs_from_R(smile, file_name="C:/Users/nakamukai/Desktop/"):
    
    import pandas as pd
    # Read bit data derived from R
    fps = smile_to_maccs(smile)

    column_names = ["maccs_" + str(i) for i in range(1, 167)] 

    row_names = ["{}".format(smile)]

    df = pd.DataFrame(fps, columns=column_names, index=row_names)

    # Read the CSV file containing column information for MACCS fingerprints
    maccs_columns = pd.read_csv("{}Fingerprints/PredRet_Nematostella/maccs_columns.csv".format(file_name))

    # Identify columns in the fingerprint that do not exist in the DataFrame
    new_columns = [col for col in maccs_columns.columns if col not in df.columns]

    # Add new columns to the fingerprint DataFrame and initialize their values to 0
    for new_col in new_columns:
        df[new_col] = 0

    # Reorder the columns of the fingerprint to match the order in maccs_columns
    fps_maccs = df[maccs_columns.columns]

    return(fps_maccs)



def convert_by_PCA(smile, file_name): 
    
    import numpy as np
    import pandas as pd
    from sklearn.preprocessing import MinMaxScaler
    from sklearn.decomposition import PCA
    
    # Read the training data
    df = pd.read_csv("{}Fingerprints/PredRet_Nematostella/PredRet_Nematostella_90_percent.csv".format(file_name), index_col=0)

    # explanatory variables(X)
    header = df.columns
    col_name = list(header)
    X = np.array(df[col_name])

    scaler = MinMaxScaler()

    scaler.fit(X)
    X_standardization = scaler.transform(X)
    
    # Get the maximum values of each column in `df`
    max_values = df.max()
    
    # Obtain and concatenate the bit sequences of the target compound
    fps_maccs = maccs_from_R(smile, file_name)
    fps_extended = fps_from_R(smile, "extended", file_name)
    fps_circular = fps_from_R(smile, "circular", file_name)
    fps_graph = fps_from_R(smile, "graph", file_name)
    df_concat = pd.concat([fps_maccs, fps_extended, fps_circular, fps_graph], axis = 1)
    
    # Replace values ​​exceeding max_value in each column of df_concat with the values ​​in the corresponding column of max_value
    list_target = [i for i in df_concat.iloc[0,:]]
    
    for col in range(len(list_target)):
        max_value = max_values[col]
        if list_target[col] > max_value: 
            list_target[col] = max_value
        else:
            pass

    df_target = pd.DataFrame(list_target, index=df_concat.columns, columns=df_concat.index.values)
    df_target = df_target.T
    X_standardization_target = scaler.transform(np.array(df_target))

    #PCA
    pca = PCA()
    pca.fit(X_standardization)
    np_pca = pca.transform(X_standardization)
    
    # Apply PCA to target fingerprints
    feature_target = pca.transform(X_standardization_target)

    index_name_target = df_target.index.values

    df_pca_target = pd.DataFrame(feature_target,  columns=["PC{}".format(x + 1) for x in range(len(df.columns))], index=index_name_target)
    
    # make every element 0 or higher than 0
    df_pca_target = df_pca_target - np.amin(np_pca)
    df_pca_target[df_pca_target < 0]  = 0
    
    return(df_pca_target)


# Function to perform retention order prediction from SMILE
def RO_prediction(smile, file_name="C:/Users/nakamukai/Desktop/"):
    
    import os 
    import pandas as pd
    import joblib
    import numpy as np

    os.chdir("{}src".format(file_name))
    print(os.getcwd())

    # Load reference compounds and prediction model
    variables = 200
    df_test = convert_by_PCA(smile, file_name)
    df_reference = pd.read_csv("{}Fingerprints/PredRet_Nematostella/PredRet_Nematostella_pca_all_test.csv".format(file_name), index_col=0)
    model_file = "{}model/model_no_peptides_200.pkl".format(file_name)
    loaded_model = joblib.load(model_file)

    df_target_variables = df_test.iloc[:, :variables]
    df_reference_variables = df_reference.iloc[:22, 2:variables + 2]

    df_reference_target = pd.concat([df_reference_variables, df_target_variables])
    
   
    

    # Create output result labels
    TRU = np.array([[ 0.,  1.], [-1.,  0.]])
    FAL = np.array([[ 0., -1.], [ 1.,  0.]])

    # Perform retention order prediction
    for i in range(len(df_reference_variables)):
        df_pair = df_reference_target.iloc[[i, 22], :]
        
        test = loaded_model.predict(df_pair, df_pair)
        if np.allclose(TRU, test):
            print(df_reference.iat[i, 0] + ": 1")
        elif np.allclose(FAL, test):
            print(df_reference.iat[i, 0] + ": 0")
        else:
            print(df_reference.iat[i, 0] + ": 2") 



