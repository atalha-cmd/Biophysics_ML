
# Biophysics Machine Learning

Welcome to the **Biophysics Machine Learning** repository! This repository contains scripts and tools for generating electrostatic and topological features, as well as training Convolutional Neural Networks (CNNs) to analyze biophysical properties using these features.

## Overview

This repository includes the following key functionalities:
- **Electrostatic Feature Generation**: A script to compute electrostatic features for biomolecules.
- **Topological Feature Generation**: A script to compute topological features.
- **CNN Model Training**: Several scripts to train CNN models using electrostatic and topological features.

## Prerequisites

Before running the scripts, ensure you have the following:
1. A **Python 3** environment.
2. A **virtual environment** set up (if not already done).
3. Required Python packages installed (can be set up in the virtual environment).

### Scripts Overview

This repository includes several scripts designed for generating features and training machine learning models for biophysics research. Below is a summary of each script's purpose and how to use them.

### 1. **Electrostatic Features Generation (`get_electrostatic_features.py`)**

This script calculates electrostatic features for biomolecules, which are essential for understanding their interaction with the environment.

**Usage:**
```bash
python get_electrostatic_features.py 0 1  # Generate electrostatic features (Arguments: p=0, L=1)
```
### 2. **Topological Features Generation (`run_all.py`)** 

This script computes topological features, which describe the shape and connectivity of the biomolecules.

**Usage:**
```bash
python run_all.py  # Generate topological features
```
### 3. **CNN Model Training** 

This script runs multiple CNN models using both electrostatic and topological features. It includes cross-validation for performance evaluation.

**Train CNN using both electrostatic and topological features:**
```bash
python CNNbothCV.py 0 1 5  # Run CNN with both features and cross-validation
```
**Train CNN using only electrostatic features:**
```bash
python CNNelectro.py 0 1 5  # Run CNN with only electrostatic features and cross-validation
```
**Train CNN using only topological features:**
```bash
python CNNtopologicalCV.py 5  # Run CNN with only topological features and cross-validation
```

