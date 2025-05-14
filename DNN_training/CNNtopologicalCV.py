#Read in train/test features
import numpy as np
import pandas as pd
import tensorflow as tf
import tensorflow.keras as keras
import os 
import re
import pickle
import matplotlib.pyplot as plt
import sys
from scipy.stats import pearsonr
from matplotlib import pyplot as plt
from numpy import savetxt
from keras import regularizers
from keras.initializers import GlorotUniform
from keras import regularizers
from keras import layers
from tensorflow.keras.models import save_model, load_model
from tensorflow.keras.layers import Reshape, concatenate
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Reshape
from tensorflow.keras.layers import Dense, Dropout, Activation, Flatten, BatchNormalization
from tensorflow.keras.layers import Conv1D, Conv2D, MaxPooling2D, AveragePooling1D
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.regularizers import l2
from tensorflow.keras.layers import Dense, Dropout, Activation, Flatten
from tensorflow.keras.layers import Conv1D, Conv2D, MaxPooling2D, AveragePooling1D
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.preprocessing.image import ImageDataGenerator
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score, mean_absolute_percentage_error 
from sklearn.preprocessing import MinMaxScaler
from sklearn.decomposition import PCA  # to apply PCA
from sklearn.model_selection import KFold



def model(k):
    
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

    print("Shape X: ", X_topological.shape)

    y_df = y_df.drop(['PDB_IDs'], axis = 1)


    histories = []
    y_true = []
    y_pred = []
    evaluation_metrics = []

    # Split the combined data into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(X_topological, y_df, test_size=0.2, random_state=42)
    
    kf = KFold(n_splits=k, shuffle=True, random_state=42)

    best_mse = float('inf')  # Initialize the best validation MSE as infinity
    # best_model_filename = "best_model.h5"     # Variable to store the filename of the best model

    counter = 1
    for train_index, val_index in kf.split(X_train):
        X_train_fold, X_val = X_train[train_index], X_train[val_index]
        y_train_fold, y_val = y_train.iloc[train_index], y_train.iloc[val_index]

        #reshape X_train and X_test for data generator
        width, height, channels = X_train_fold.shape[1], X_train_fold.shape[2], 1
        X_train_fold = X_train_fold.reshape((X_train_fold.shape[0], width, height, channels))
        width, height, channels = X_test.shape[1], X_test.shape[2], 1
        X_test = X_test.reshape((X_test.shape[0], width, height, channels))
        width, height, channels = X_val.shape[1], X_val.shape[2], 1
        X_val = X_val.reshape((X_val.shape[0], width, height, channels))

        # # report pixel means and standard deviations
        
        # Fit the datagen on the training set and use the same generator for test and validation sets
        # Center training data
        # create generator that centers pixel values
        datagen = ImageDataGenerator(featurewise_center=True, featurewise_std_normalization=True)
        datagen.fit(X_train_fold)

        # Transform training data
        train_iterator = datagen.flow(X_train_fold, batch_size=len(X_train_fold), shuffle=False)
        X_train_centered = train_iterator.__next__()

        # Transform test data using training set statistics
        test_iterator = datagen.flow(X_test, batch_size=len(X_test), shuffle=False)
        X_test_centered = test_iterator.__next__()

        # Transform validation data using training set statistics
        val_iterator = datagen.flow(X_val, batch_size=len(X_val), shuffle=False)
        X_val_centered = val_iterator.__next__()

        # Print stats after transformation
        print(X_train_centered.shape, X_train_centered.mean(), X_train_centered.std())
        print(X_test_centered.shape, X_test_centered.mean(), X_test_centered.std())
        print(X_val_centered.shape, X_val_centered.mean(), X_val_centered.std())


        # Reshape both
        X_train_centered = X_train_centered.reshape((X_train_centered.shape[0], X_train.shape[1], X_train.shape[2]))
        X_test_centered = X_test_centered.reshape((X_test_centered.shape[0], X_test.shape[1], X_test.shape[2]))
        X_val_centered = X_val_centered.reshape((X_val_centered.shape[0], X_val.shape[1], X_val.shape[2]))

        # SCALE y
        scaler_y = StandardScaler()
        y_train_scaled = scaler_y.fit_transform(np.array(y_train_fold).reshape(-1, 1)).flatten()
        y_test_scaled = scaler_y.transform(np.array(y_test).reshape(-1, 1)).flatten()
        y_val_scaled = scaler_y.transform(np.array(y_val).reshape(-1, 1)).flatten()

        inputCNN = keras.Input(shape = X_train_centered.shape[1:])

        # CNN model
        x = Conv1D(filters=128, kernel_size=3, activation='tanh', kernel_regularizer=l2(0.02))(inputCNN)
        x = AveragePooling1D(pool_size=2)(x)
        x = Dropout(0.3)(x)
        x = BatchNormalization()(x)
        x = Conv1D(filters=256, kernel_size=3, activation='relu', kernel_regularizer=l2(0.02))(x)
        x = AveragePooling1D(pool_size=2)(x)
        x = Dropout(0.3)(x)
        x = BatchNormalization()(x)
        x = Conv1D(filters=128, kernel_size=3, activation='tanh', kernel_regularizer=l2(0.02))(x)
        x = AveragePooling1D(pool_size=2)(x)
        x = Dropout(0.3)(x)
        x = Flatten()(x)
        x = Dense(64, activation='relu')(x)
        x = Dropout(0.15)(x)
        x = Dense(1)(x)

        modelCNN = keras.Model(inputs=inputCNN, outputs=x)

        adam = Adam(learning_rate = 0.0001)
        modelCNN.compile(loss = 'mean_squared_error', optimizer = adam, metrics=['mean_squared_error'])

        history = modelCNN.fit(X_train_centered, y_train_scaled, batch_size = 16, epochs=500, 
                        validation_data = (X_val_centered, y_val_scaled), verbose=0)

        histories.append(history)

        
        y_val_pred_scaled = modelCNN.predict(X_val_centered)
        y_val_pred = scaler_y.inverse_transform(y_val_pred_scaled)
        
        # Append unscaled y_val and y_val_pred to y_true and y_pred lists
        print("Filling y_true with y_val", y_val.shape)
        print("Filling y_pred with y_val_pred", y_val_pred.shape)

        y_true.append(y_val)
        y_pred.append(y_val_pred)

        # Calculate evaluation metrics
        mse = mean_squared_error(y_val_scaled, y_val_pred_scaled)
        mape = mean_absolute_percentage_error(y_val, y_val_pred)
        r2 = r2_score(y_val, y_val_pred)
        cc = np.corrcoef(np.array(y_val).squeeze(), np.array(y_val_pred).squeeze())[0, 1]
        
        print("Evaluation Metrics for Validation Data for k = ", counter)
        print(f"Mean Squared Error (MSE) scaled: {mse}")
        print(f"Mean Absolute Percentage Error (MAPE): {mape}")
        print(f"R-squared (R2) Score: {r2}")
        print(f"Correlation Coefficient: {cc}")

        # Save the model if it has the best MSE so far
        if mse < best_mse:
            best_mse = mse
            filename = f"best_model_fold_{counter}.keras"
            modelCNN.save(filename)
            print(f"\n New best model saved with MSE: {best_mse} at fold {counter} \n")

        evaluation_metrics.append((mse, mape, r2, cc))
        print("evaluation metrics: \n", evaluation_metrics)
        counter = counter + 1 

    # Calculate mean evaluation metrics across all folds
    mean_metrics = np.mean(evaluation_metrics, axis=0)
    print("\n mean evaluation metrics: ", mean_metrics)


    print(f"\n Best Model MSE: {best_mse}, saved as {filename}\n ")

    # Load the best model
    best_model = load_model(filename)
    best_model.compile(
        optimizer=Adam(learning_rate=0.0001),  # Same optimizer as during training
        loss='mean_squared_error',            # Same loss function
        metrics=['mse', 'mape'])     # Same metrics
    

    # y_pred_scaled = modelCNN.predict(X_test_centered)
    y_pred_scaled = best_model.predict(X_test_centered)
    y_pred_inv = scaler_y.inverse_transform(y_pred_scaled)
    y_test = np.array(y_test)

    # Evaluate on the test set
    test_mse_scaled = mean_squared_error(y_test_scaled, y_pred_scaled)
    test_rmse_scaled = np.sqrt(test_mse_scaled)
    test_mape = mean_absolute_percentage_error(y_test, y_pred_inv)
    test_cc = np.corrcoef(np.array(y_test).squeeze(), np.array(y_pred_inv).squeeze())[0, 1]
    test_r2 = r2_score(y_test, y_pred_inv)

    print("y_test:", y_test[0:10])
    print("y_pred_inv:", y_pred_inv[0:10])

    return histories, y_test, y_pred_inv, mean_metrics, test_mse_scaled, test_rmse_scaled, test_mape, test_cc, test_r2


def plot_mean_loss(histories):
    mean_train_loss = np.mean([history.history['loss'] for history in histories], axis=0)
    mean_val_loss = np.mean([history.history['val_loss'] for history in histories], axis=0)

    plt.figure(figsize=(9,6))
    plt.plot(mean_train_loss, label='Mean Training Loss')
    plt.plot(mean_val_loss, label='Mean Validation Loss')
    plt.xlabel('Epochs', fontsize = 20)
    plt.ylabel('Loss', fontsize = 20)
    # plt.title('Model loss across ' + str(k) + ' folds \nUsing Topological Features', fontsize = 20)
    plt.legend(fontsize=18)  
    plt.xticks(fontsize=14)  
    plt.yticks(fontsize=14)  
    plt.tight_layout()
    plt.savefig('plots/loss_topological.png')

def plot_scatter(y_true, y_pred):
    y_true = np.array(y_true)
    print("y_true in scatter", y_true[0:10])
    y_pred = np.array(y_pred) 
    print("y_pred in scatter", y_pred[0:10])

    plt.figure(figsize=(9,6))
    plt.scatter(y_true, y_pred, marker='o', facecolors='none', edgecolors='b')
    plt.plot([min(y_true.min(), y_pred.min()), max(y_true.max(), y_pred.max())], [min(y_true.min(), y_pred.min()), max(y_true.max(), y_pred.max())], color='red')
    plt.xlabel('Reference Values', fontsize=20)
    plt.ylabel('Predicted Values', fontsize = 20)
    # plt.title('True vs Predicted Values \nUsing Topological Features', fontsize = 20)
    plt.xticks(fontsize=14)  
    plt.yticks(fontsize=14)
    plt.tight_layout()
    plt.savefig('plots/scatter_topological.png')
  
if __name__ == "__main__":
    # Extract command-line arguments
    k = int(sys.argv[1])

    histories, y_true, y_pred, mean_metrics,  test_mse_scaled,  test_rmse_scaled,  test_mape, test_cc, test_r2 = model(k)
    
    plot_mean_loss(histories)
 
    plot_scatter(y_true, y_pred)
    # Save output data
    pd.DataFrame(y_true).to_csv('outputs/y_true_topological.csv')
    pd.DataFrame(y_pred).to_csv('outputs/y_pred_topological.csv')

    metrics_data = {
         "Metric": [
              "Mean Squared Error (MSE) scaled",
              "Root Mean Squared Error (MSE) scaled",
              "Mean Absolute Percentage Error (MAPE)",
              "Correlation coefficient",
              "R-squared (R2) Score"],
        "Value": [
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
    df_metrics.to_csv("evals/evaluation_metrics_topological.txt", sep='\t', index=False)




