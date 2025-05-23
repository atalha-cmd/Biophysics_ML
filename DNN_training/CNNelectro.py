#Read in train/test features
import pickle
from matplotlib import pyplot as plt
import numpy as np
from numpy import savetxt
from sklearn.preprocessing import StandardScaler
import pandas as pd
import os 
import re
from tensorflow.keras.preprocessing.image import ImageDataGenerator
from sklearn.model_selection import train_test_split

import tensorflow.keras as keras
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Reshape
from tensorflow.keras.layers import Dense, Dropout, Activation, Flatten, BatchNormalization
from tensorflow.keras.layers import Conv1D, Conv2D, MaxPooling2D, AveragePooling1D
from tensorflow.keras.optimizers import Adam
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score, mean_absolute_percentage_error 
from keras import regularizers
#from scikeras.wrappers import KerasClassifier
from tensorflow.keras.regularizers import l2

from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from numpy import savetxt
import os 
import re
import tensorflow as tf
import tensorflow.keras as keras
from tensorflow.keras.layers import Dense, Dropout, Activation, Flatten
from tensorflow.keras.layers import Conv1D, Conv2D, MaxPooling2D, AveragePooling1D
from tensorflow.keras.optimizers import Adam
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from keras.initializers import GlorotUniform
from keras import regularizers
from keras import layers
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score, mean_absolute_percentage_error 
from sklearn.preprocessing import MinMaxScaler
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA  # to apply PCA
from sklearn.model_selection import KFold
import sys
from scipy.stats import pearsonr
import sys
from tensorflow.keras.layers import Input

# Read in  features (generated by running run_all.py)
def model(p, L, k):

    #  Get the labels
    comps_df = pd.read_csv('comp_17k_list.txt', sep='\s+', header = None, names = ['PDB_IDs'])
    print(comps_df)
    New_Set_labels = pd.read_csv('17k_CE_labels.txt', sep=',', header = None, names = ['PDB_IDs', 'CE'])
    New_Set_labels = New_Set_labels.drop_duplicates(subset=['PDB_IDs']).reset_index(drop=True)
    print(New_Set_labels)

    y_df = pd.merge(comps_df, New_Set_labels, on='PDB_IDs', how='inner')

    print(y_df)
    print(y_df.shape)
    print("number of elements in labels:",len(y_df))

    # Now load in electrostatic features

    X_electrostatic_df =  pd.read_csv('X_electrostatic_p' + str(p) + '_L' + str(L) + '.csv', index_col=0)

    # Line everything up based on protein IDs
    merged_df = pd.merge(X_electrostatic_df, y_df, on='PDB_IDs', how = 'inner') 
    merged_df = merged_df.drop('PDB_IDs', axis=1)
    merged_df.dropna(subset=['CE'], inplace=True)
    mask = merged_df['CE'].isin([np.inf, -np.inf])
    merged_df = merged_df[~mask]

    # Get rid of rows where the entire feature is equal to zero
    columns_to_check = merged_df.iloc[:, 1:] 
    are_rows_zero = (columns_to_check == 0).all(axis = 1)  # Check if all values in these columns are 0
    merged_df = merged_df[~are_rows_zero]
    print("size of merged_df after cleaning: ", merged_df.shape)


    X = merged_df.drop(['CE'], axis=1)
    y = merged_df['CE']
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
        
    scaler_X = StandardScaler()
    scaler_y = StandardScaler()
    X_train_scaled = scaler_X.fit_transform(X_train)
    X_test_scaled = scaler_X.transform(X_test)
    y_train_scaled = scaler_y.fit_transform(np.array(y_train).reshape(-1, 1)).flatten()
    y_test_scaled = scaler_y.transform(np.array(y_test).reshape(-1, 1)).flatten()

    print("y:", y)
    # Resetting indices
    X_train.reset_index(drop=True, inplace=True)
    y_train.reset_index(drop=True, inplace=True)

    print("X_train.shape", X_train.shape)
    print("X_test.shape", X_test.shape)
    print("y_train.shape", y_train.shape)
    print("y_test.shape", y_test.shape)
       
       
    kf = KFold(n_splits=k, shuffle=True, random_state=42)

    histories = []
    y_true = []
    y_pred = []
    evaluation_metrics = []

    counter = 1
    for train_index, val_index in kf.split(X_train):
        X_train_fold, X_val = X_train.iloc[train_index], X_train.iloc[val_index]
        y_train_fold, y_val = y_train.iloc[train_index], y_train.iloc[val_index]

        print(X_train_fold.shape)
        print(X_val.shape)

        print(y_train_fold.shape)
        print(y_val.shape)

        scaler_X_fold = StandardScaler()
        scaler_y_fold = StandardScaler()
        X_train_fold_scaled = scaler_X_fold.fit_transform(X_train_fold)
        X_val_scaled = scaler_X_fold.transform(X_val)
        y_train_fold_scaled = scaler_y_fold.fit_transform(np.array(y_train_fold).reshape(-1, 1)).flatten()
        y_val_scaled = scaler_y_fold.transform(np.array(y_val).reshape(-1, 1)).flatten()

        # Sequential Model
        model = Sequential()
        model.add(Input(shape=(X.shape[1],)))
        model.add(layers.Dense(128, activation="relu",kernel_initializer='he_uniform', kernel_regularizer=regularizers.l2(0.01)))
        model.add(Dropout(0.15))
        model.add(layers.Dense(64))
        model.add(BatchNormalization())
        model.add(tf.keras.layers.Activation('relu'))
        model.add(Dropout(0.15))
        model.add(layers.Dense(32))
        model.add(BatchNormalization())
        model.add(tf.keras.layers.Activation('relu'))
        model.add(Dropout(0.15))
        model.add(layers.Dense(16))
        model.add(BatchNormalization())
        model.add(tf.keras.layers.Activation('relu'))
        model.add(Dropout(0.15))
        model.add(layers.Dense(8))
        model.add(BatchNormalization())
        model.add(tf.keras.layers.Activation('relu'))
        model.add(Dropout(0.15))
        model.add(layers.Dense(4))
        model.add(BatchNormalization())
        model.add(tf.keras.layers.Activation('relu'))
        model.add(Dropout(0.15))
        model.add(layers.Dense(1))
        
        adam = Adam(learning_rate = 0.0001)

        model.compile(loss = 'mean_squared_error', optimizer = adam, metrics=['mean_squared_error'])   
                
        history = model.fit(X_train_fold_scaled, y_train_fold_scaled, batch_size = 32, epochs=500, 
                            validation_data = (X_val_scaled, y_val_scaled), verbose=1)
                
        histories.append(history)

        y_val_pred_scaled = model.predict(X_val_scaled)
        y_val_pred = scaler_y_fold.inverse_transform(y_val_pred_scaled)
        
        # Append unscaled y_val and y_val_pred to y_true and y_pred lists
        print("Filling y_true with y_val", y_val.shape)
        print("Filling y_pred with y_val_pred", y_val_pred.shape)

        y_true.append(y_val)
        y_pred.append(y_val_pred)

        print("Filling y_true with y_val", len(y_true))
        print("Filling y_pred with y_val_pred", len(y_pred))

        # Calculate evaluation metrics
        mse = mean_squared_error(y_val, y_val_pred)
        mape = mean_absolute_percentage_error(y_val, y_val_pred)
        r2 = r2_score(y_val, y_val_pred)
        cc = np.corrcoef(np.array(y_val).squeeze(), np.array(y_val_pred).squeeze())[0, 1]
        
        print("Evaluation Metrics for Validation Data for k = ", counter)
        print(f"Mean Squared Error (MSE): {mse}")
        print(f"Mean Absolute Error (MAE): {mape}")
        print(f"R-squared (R2) Score: {r2}")
        print(f"Correlation Coefficient: {cc}")

        evaluation_metrics.append((mse, mape, r2, cc))
        counter = counter + 1
        print("evaluation metrics: ", evaluation_metrics)

    # Calculate mean evaluation metrics across all folds
    mean_metrics = np.mean(evaluation_metrics, axis=0)
    print("mean evaluation metrics: ", mean_metrics)

    # Evaluate on the test set
    X_test_scaled = scaler_X.transform(X_test)
    y_pred_scaled = model.predict(X_test_scaled)
    y_pred_inv = scaler_y.inverse_transform(y_pred_scaled)
    
    test_mse_scaled = mean_squared_error(y_test_scaled, y_pred_scaled)
    test_rmse_scaled = np.sqrt(test_mse_scaled)
    test_mape = mean_absolute_percentage_error(y_test, y_pred_inv)
    test_cc = np.corrcoef(np.array(y_test).squeeze(), np.array(y_pred_inv).squeeze())[0, 1]
    test_r2 = r2_score(y_test, y_pred_inv)
    

    return histories, y_test, y_pred_inv, mean_metrics, test_mse_scaled, test_rmse_scaled, test_mape, test_cc, test_r2


def plot_mean_loss(histories):
    mean_train_loss = np.mean([history.history['loss'] for history in histories], axis=0)
    mean_val_loss = np.mean([history.history['val_loss'] for history in histories], axis=0)

    plt.figure(figsize=(9,6))
    plt.plot(mean_train_loss, label='Mean Training Loss')
    plt.plot(mean_val_loss, label='Mean Validation Loss')
    plt.xlabel('Epochs', fontsize = 20)
    plt.ylabel('Loss', fontsize = 20)
    # plt.title('Model loss across ' + str(k) + ' folds for $p$ = ' + str(p) + ' and $L$ = ' + str(L)+"\nUsing Electrostatic Features", fontsize = 20)
    plt.legend(fontsize=18)  
    plt.xticks(fontsize=14)  
    plt.yticks(fontsize=14) 
    plt.tight_layout()
    plt.savefig('loss_electro_p' + str(p) + '_L' + str(L) + '.png')

def plot_scatter(y_true, y_pred):
    y_true = np.array(y_true)
    y_pred = np.array(y_pred) 
    plt.figure(figsize=(9,6))
    plt.scatter(y_true, y_pred, marker='o', facecolors='none', edgecolors='b')
    plt.plot([min(y_true.min(), y_pred.min()), max(y_true.max(), y_pred.max())], [min(y_true.min(), y_pred.min()), max(y_true.max(), y_pred.max())], color='red')
    plt.xlabel('Reference Values',fontsize=20)
    plt.ylabel('Predicted Values',fontsize = 20)
    # plt.title('True vs Predicted Values for $p$ = ' + str(p) + ' and $L$ = ' + str(L) + " \nUsing Electrostatic Features",fontsize = 20)
    plt.xticks(fontsize=14)  
    plt.yticks(fontsize=14)
    plt.tight_layout()
    plt.savefig('scatter_electro_p' + str(p) + '_L' + str(L) + '.png')
  
if __name__ == "__main__":
    # Extract command-line arguments
    p = int(sys.argv[1])
    L = int(sys.argv[2])
    k = int(sys.argv[3])

    histories, y_true, y_pred, mean_metrics, test_mse_scaled, test_rmse_scaled, test_mape, test_cc, test_r2 = model(p, L, k)
    
    plot_mean_loss(histories)
 
    plot_scatter(y_true, y_pred)

    #  Save output data   
    pd.DataFrame(y_true).to_csv('y_true_electro_p' + str(p) + '_L' + str(L) + '.csv')
    pd.DataFrame(y_pred).to_csv('y_pred_electro_p' + str(p) + '_L' + str(L) + '.csv')

    
    metrics_data = {
         "Metric": [
              "p",
              "L",
              "Mean Squared Error (MSE) Scaled",
              "Root Mean Squared Error (RMSE) Scaled",
              "Mean Absolute Percentage Error (MAPE)",
              "Correlation coefficient",
              "R-squared (R2) Score"],

        "Value": [
             p,
             L,
             test_mse_scaled,
             test_rmse_scaled,
             test_mape,
             test_cc,
             test_r2]}


    # Create the DataFrame
    df = pd.DataFrame(metrics_data)

    # Transpose the DataFrame
    df_metrics = df.set_index("Metric").T

    print(df_metrics)


    # Save the DataFrame to a text file
    df_metrics.to_csv("evaluation_metrics_p" + str(p) + "_L" + str(L) + "_CE.txt", sep='\t', index=False)




