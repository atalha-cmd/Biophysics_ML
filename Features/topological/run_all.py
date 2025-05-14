# THIS CODE IS USED TO GENERATE TOPOLOGICAL FEATURES

import time
from Get_structure import Get_structure
from Run_alpha import Run_alpha
from Run_alpha_hydro import Run_alpha_hydro
from PrepareData import PrepareData
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import os 
import re
import pickle


# Read in labels from MIBPB_Sol_PQRs folder

csv_filename = "comp_test_protein.txt"
# csv_filename = "comp_test_file.txt"

# Read data from the CSV file into a DataFrame
df = pd.read_csv(csv_filename, header = None, names = ['PDB_IDs'])
print(df)
print(df['PDB_IDs'].shape)
# Now generate 1D image-like input for each complex in the New Set
# Only the New set has .pdb files to generate the topological features

# list of output files
new_set_features = []
new_total_time = 0
counter = 0


for comp in df['PDB_IDs']:
  start_time = time.time()
  print(counter + 1)
  print(comp)

  dir = "/users/atalha/Biophysics_ML/Data/test_protein/" + comp
  # dir = "/users/atalha/HPC_test/test_data/" + comp

  os.chdir(dir) # change to directory

  pdb_file = [x for x in os.listdir(dir) if x.endswith('.pdb')][0]


  Get_structure(pdb_file,'complex.npz')
  Run_alpha('complex.npz', 'protein.PH')
  Run_alpha_hydro('complex.npz', 'protein_C-C.PH')
  PrepareData('protein.PH', 'protein_C-C.PH', 'complex_digit.npy')
  
  # Final Output
  feature = np.load('complex_digit.npy') 
  print("feature size", feature.shape)
  new_set_features.append(feature)

  # Computing Time
  end_time = time.time()
  total_time = (end_time - start_time)
  print("total time: ", total_time, " seconds")
  new_total_time = new_total_time + total_time
  counter = counter + 1

print("test_total_time:", new_total_time)

print("final shape:", len(new_set_features ))

# now save 
dir = "/users/atalha/Biophysics_ML/Features/topological"
os.chdir(dir) # change to directory
with open("TDA_test_features.txt", "wb") as fp:   #Pickling
    pickle.dump(new_set_features, fp)

