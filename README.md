
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
source /users/atalha/.venv/bin/activate  # Activate the virtual environment
python get_electrostatic_features.py 0 1  # Generate electrostatic features (Arguments: p=0, L=1)

