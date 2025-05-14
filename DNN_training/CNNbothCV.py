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
from tensorflow.keras.layers import Reshape, concatenate
import tensorflow.keras as keras
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Reshape
from tensorflow.keras.layers import Dense, Dropout, Activation, Flatten, BatchNormalization
from tensorflow.keras.layers import Conv1D, Conv2D, MaxPooling2D, AveragePooling1D
from tensorflow.keras.optimizers import Adam
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score, mean_absolute_percentage_error 
from keras import regularizers
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


def model(p, L, k):

    pickle_file = open("TDA_17k_features.txt", "rb")
    features = []
    while True:
        try:
            features.append(pickle.load(pickle_file))
        except EOFError:
            break
    pickle_file.close()

    X_topological = np.asarray(features[0])

    # Now get the labels
    comps_df = pd.read_csv('comp_17k_list.txt', sep='\s+', header = None, names = ['PDB_IDs'])
    print(comps_df)
    New_Set_labels = pd.read_csv('17k_CE_labels.txt', sep=',', header = None, names = ['PDB_IDs', 'CE'])
    New_Set_labels = New_Set_labels.drop_duplicates(subset=['PDB_IDs']).reset_index(drop=True)
    print(New_Set_labels)

    y_df = pd.merge(comps_df, New_Set_labels, on='PDB_IDs', how='inner')
    
    print(y_df.shape)
    print("number of elements in labels:",len(y_df))

    X_topological = np.asarray(features[0])
    print("Shape X: ", X_topological.shape)

    # Now load in electrostatic features
    X_electrostatic_df =  pd.read_csv('X_electrostatic_p' + str(p) + '_L' + str(L) + '.csv', index_col=0)

    y_df = y_df.drop(['PDB_IDs'], axis = 1)
    X_electrostatic_df = X_electrostatic_df.drop(['PDB_IDs'], axis = 1)

    X_topological_reshaped = X_topological.reshape(X_topological.shape[0], -1) # Flatten the image features
    X_combined = np.concatenate((X_topological_reshaped, np.array(X_electrostatic_df)), axis=1)

    print(X_topological.shape)
    print(X_electrostatic_df.shape)

    histories = []
    y_true = []
    y_pred = []
    evaluation_metrics = []
    # Split the combined data into training and testing sets
    X_train_combined, X_test_combined, y_train, y_test = train_test_split(X_combined, y_df, test_size=0.2, random_state=42)
    
    kf = KFold(n_splits=k, shuffle=True, random_state=42)
    counter = 1
    
    for train_index, val_index in kf.split(X_train_combined):

        X_train_fold, X_val = X_train_combined[train_index], X_train_combined[val_index]
        y_train_fold, y_val = y_train.iloc[train_index], y_train.iloc[val_index]

        # Get topological features from fold data
        X_train_images_fold = X_train_fold[:, :X_topological_reshaped.shape[1]].reshape(-1, 100, 12, 1) 
        X_val_images = X_val[:, :X_topological_reshaped.shape[1]].reshape(-1, 100, 12, 1) 
        X_test_images = X_test_combined[:, :X_topological_reshaped.shape[1]].reshape(-1, 100, 12, 1)
 
        
        # Center train images
        datagen = ImageDataGenerator(featurewise_center=True, featurewise_std_normalization=True)
        # calculate the mean on the training dataset
        datagen.fit(X_train_images_fold)
        # demonstrate effect on entire training dataset
        iterator = datagen.flow(X_train_images_fold, batch_size=len(X_train_images_fold), shuffle=False)
        # get a batch
        X_train_centered = iterator.__next__()
        # pixel stats in the batch
        print(X_train_centered.shape, X_train_centered.mean(), X_train_centered.std())

        # Center val images
        # create generator that centers pixel values
        datagen = ImageDataGenerator(featurewise_center=True, featurewise_std_normalization=True)
        # calculate the mean on the training dataset
        datagen.fit(X_val_images)
        # demonstrate effect on entire training dataset
        iterator = datagen.flow(X_val_images, batch_size=len(X_val_images), shuffle=False)
        # get a batch
        X_val_centered = iterator.__next__()
        # pixel stats in the batch
        print(X_val_centered.shape, X_val_centered.mean(), X_val_centered.std())

        # Center test images
        # create generator that centers pixel values
        datagen = ImageDataGenerator(featurewise_center=True, featurewise_std_normalization=True)
        # calculate the mean on the training dataset
        datagen.fit(X_test_images)
        # demonstrate effect on entire training dataset
        iterator = datagen.flow(X_test_images, batch_size=len(X_test_images), shuffle=False)
        # get a batch
        X_test_centered = iterator.__next__()
        # pixel stats in the batch
        print(X_test_centered.shape, X_test_centered.mean(), X_test_centered.std())

        # Reshape 
        X_train_centered = X_train_centered.reshape((X_train_centered.shape[0], X_train_centered.shape[1], X_train_centered.shape[2]))
        X_test_centered = X_test_centered.reshape((X_test_centered.shape[0], X_test_centered.shape[1], X_test_centered.shape[2]))
        X_val_centered = X_val_centered.reshape((X_val_centered.shape[0], X_val_centered.shape[1], X_val_centered.shape[2]))

        scaler_X = StandardScaler()
        X_train_ef_scaled = scaler_X.fit_transform(X_train_fold[:, X_topological_reshaped.shape[1]:])
        X_test_ef_scaled = scaler_X.transform(X_test_combined[:, X_topological_reshaped.shape[1]:])
        X_val_ef_scaled = scaler_X.transform(X_val[:, X_topological_reshaped.shape[1]:])

        # SCALE y
        scaler_y = StandardScaler()
        y_train_scaled = scaler_y.fit_transform(np.array(y_train_fold).reshape(-1, 1)).flatten()
        y_test_scaled = scaler_y.transform(np.array(y_test).reshape(-1, 1)).flatten()
        y_val_scaled = scaler_y.transform(np.array(y_val).reshape(-1, 1)).flatten()

        inputCNN = keras.Input(shape = X_train_centered.shape[1:])
        inputEF = keras.Input(shape = (X_train_ef_scaled.shape[1],)) 

        # CNN model
        x = Conv1D(filters=128, kernel_size=3, activation='tanh', kernel_regularizer=l2(0.02))(inputCNN)
        x = AveragePooling1D(pool_size=2)(x)
        x = Dropout(0.3)(x)
        x = BatchNormalization()(x)
        x = Conv1D(filters=64, kernel_size=3, activation='relu', kernel_regularizer=l2(0.02))(x)
        x = AveragePooling1D(pool_size=2)(x)
        x = Dropout(0.3)(x)
        x = BatchNormalization()(x)
        x = Conv1D(filters=64, kernel_size=3, activation='tanh', kernel_regularizer=l2(0.02))(x)
        x = AveragePooling1D(pool_size=2)(x)
        x = Dropout(0.3)(x)
        x = Flatten()(x)
        x = Dense(64, activation='relu')(x)
        modelCNN = keras.Model(inputs=inputCNN, outputs=x)

        # EF model
        y = Dense(128, activation='relu', kernel_initializer='he_uniform', kernel_regularizer=l2(0.01))(inputEF)
        y = Dropout(0.15)(y)
        y = BatchNormalization()(y)
        y = Dense(64, activation='relu', kernel_regularizer=l2(0.02))(y)
        modelEF = keras.Model(inputs=inputEF, outputs=y)

        # Merged model
        mergedOutput = concatenate([modelCNN.output, modelEF.output])  
        z = Dense(32, activation='relu', kernel_initializer='he_normal', kernel_regularizer=l2(0.01))(mergedOutput)
        z = Dropout(0.15)(z)
        z = BatchNormalization()(z)
        z = Dense(16, activation='relu', kernel_regularizer=l2(0.02))(z)
        z = Dropout(0.15)(z)
        z = BatchNormalization()(z)
        z = Dense(8, activation='relu', kernel_regularizer=l2(0.02))(z)
        finalOutput = Dense(1)(z)
        model_merged = keras.Model(inputs=[modelCNN.input, modelEF.input], outputs=finalOutput)

        adam = Adam(learning_rate = 0.0001)
        model_merged.compile(loss = 'mean_squared_error',
                        optimizer = adam,
                        metrics=['mean_squared_error'])

        history = model_merged.fit([X_train_centered, X_train_ef_scaled], y_train_scaled, batch_size = 16, epochs=5,
                        validation_data = ([X_val_centered, X_val_ef_scaled], y_val_scaled), verbose=1)

        histories.append(history)

        y_val_pred_scaled = model_merged.predict([X_val_centered, X_val_ef_scaled])
        y_val_pred = scaler_y.inverse_transform(y_val_pred_scaled)
        
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
        print(f"Mean Absolute Percentage Error (MAE): {mape}")
        print(f"R-squared (R2) Score: {r2}")
        print(f"Correlation Coefficient: {cc}")

        evaluation_metrics.append((mse, mape, r2, cc))
        counter = counter + 1
    
    print("evaluation metrics: ", evaluation_metrics)    
    # Calculate mean evaluation metrics across all folds
    mean_metrics = np.mean(evaluation_metrics, axis=0)

    # Evaluate on the test set
    y_pred_scaled = model_merged.predict([X_test_centered, X_test_ef_scaled])
    y_pred_inv = scaler_y.inverse_transform(y_pred_scaled)
    y_test = np.array(y_test)

    test_mse = mean_squared_error(y_test, y_pred_inv)
    test_mse_scaled = mean_squared_error(y_test_scaled, y_pred_scaled)
    test_rmse = np.sqrt(test_mse)
    test_rmse_scaled = np.sqrt(test_mse_scaled)
    test_mae = mean_absolute_error(y_test, y_pred_inv)
    test_mape = mean_absolute_percentage_error(y_test, y_pred_inv)
    test_cc = np.corrcoef(np.array(y_test).squeeze(), np.array(y_pred_inv).squeeze())[0, 1]
    test_r2 = r2_score(y_test, y_pred_inv)
    pearson_corr, _ = pearsonr(np.array(y_test).squeeze(), np.array(y_pred_inv).squeeze())
    print("Pearson correlation coefficient for test set (R_p):", pearson_corr)

    return histories, y_test, y_pred_inv, mean_metrics, test_mse, test_mse_scaled, test_rmse, test_rmse_scaled, test_mae, test_mape, test_cc, test_r2

def plot_mean_loss(histories):
    mean_train_loss = np.mean([history.history['loss'] for history in histories], axis=0)
    mean_val_loss = np.mean([history.history['val_loss'] for history in histories], axis=0)

    plt.figure(figsize=(9,6))
    plt.plot(mean_train_loss, label='Mean Training Loss')
    plt.plot(mean_val_loss, label='Mean Validation Loss')
    plt.xlabel('Epochs', fontsize = 20)
    plt.ylabel('Loss', fontsize = 20)
    # plt.title('Model loss across ' + str(k) + ' folds for $p$ = ' + str(p) + ' and $L$ = ' + str(L)+"\nUsing Electrostatic & Topological Features", fontsize = 20)
    plt.legend(fontsize=18)  
    plt.xticks(fontsize=14)  
    plt.yticks(fontsize=14)      
    plt.tight_layout()
    plt.savefig('loss_both_p' + str(p) + '_L' + str(L) + '.png')

def plot_scatter(y_true, y_pred):
    y_true = np.array(y_true)
    y_pred = np.array(y_pred) 
    plt.figure(figsize=(9,6))
    plt.scatter(y_true, y_pred, marker='o', facecolors='none', edgecolors='b')
    plt.plot([min(y_true.min(), y_pred.min()), max(y_true.max(), y_pred.max())], [min(y_true.min(), y_pred.min()), max(y_true.max(), y_pred.max())], color='red')
    plt.xlabel('Reference Values', fontsize=20)
    plt.ylabel('Predicted Values',fontsize = 20)
    # plt.title('True vs Predicted Values for $p$ = ' + str(p) + ' and $L$ = ' + str(L)  + "\nUsing Electrostatic & Topological Features", fontsize = 20)
    plt.xticks(fontsize=14)  
    plt.yticks(fontsize=14)
    plt.tight_layout()
    plt.savefig('scatter_both_p' + str(p) + '_L' + str(L) + '.png')

  
if __name__ == "__main__":
    # Extract command-line arguments
    p = int(sys.argv[1])
    L = int(sys.argv[2])
    k = int(sys.argv[3])

    histories, y_true, y_pred, mean_metrics,  test_mse, test_mse_scaled, test_rmse, test_rmse_scaled, test_mae, test_mape, test_cc, test_r2 = model(p, L, k)
    
    plot_mean_loss(histories)
 
    plot_scatter(y_true, y_pred)

    #  Save output data   
    pd.DataFrame(y_true).to_csv('y_true_both_p' + str(p) + '_L' + str(L) + '.csv')
    pd.DataFrame(y_pred).to_csv('y_pred_both_p' + str(p) + '_L' + str(L) + '.csv')

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

    # Create a DataFrame
    df = pd.DataFrame(metrics_data)

    # Transpose the DataFrame
    df_metrics = df.set_index("Metric").T

    print(df_metrics)

    # Save the DataFrame to a text file
    df_metrics.to_csv("evaluation_metrics_both_p" + str(p) + "_L" + str(L) + "_CE.txt", sep='\t', index=False)




