"""

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  ML-STIM: Machine Learning for SubThalamic nucleus Intraoperative Mapping  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Python script for processing and classifying MicroElectrode Recordings (MERs). 
This script loads MERs, processes it (applies filters and removes artifacts), 
extracts features, and uses a pre-trained MLP model to make predictions. The 
results are then tabulated and saved to CSV files.

Author(s): Fabrizio SCISCENTI (fabrizio.sciscenti@polito.it)
           PolitoBIOMed Lab and BIOLAB, Politecnico di Torino, Turin, Italy

           Marco GHISLIERI (marco.ghislieri@polito.it)
           PolitoBIOMed Lab and BIOLAB, Politecnico di Torino, Turin, Italy 

Last Update: 15-05-2025

Functions:
----------
load_data(filepath, metapath)
    Loads MERs data .npz and metadata .csv from specified file paths.
load_model(model_path)
    Loads a pre-trained MLP model from the specified model path.
tabulate_results(predictions, dts, meta)
    Tabulates the predictions and processing times.
process_recording(model, recording, fsamp, b, a)
    Processes a single MER through filtering, artifact removal, feature 
    extraction,     and prediction using the MLP.
main(model, raw_data, lens, fsamp, b, a, partition_count, num_partitions)
    Main processing function that processes recordings in batches and 
    aggregates results.

Files:
------
    testing.py
    lib.py

"""

# import necessary libraries
# --------------------------
import os
import numpy as np
import pandas as pd
import torch
import time
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
from scipy.signal import resample
import lib
import warnings
warnings.filterwarnings("ignore")

# Check for GPU availability
# --------------------------

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
if device.type == 'cuda':
    print("Running on GPU!")
else:
    print("GPU not detected. Running on CPU.")

# Load data
# ---------
def load_data(filepath, metapath):
    with np.load(filepath) as npfh:
        data = npfh['data']
        meta = pd.read_csv(metapath, sep=';')
        lens = meta['length'].to_numpy()
        clss = meta['class'].to_numpy()

    print(f"Data loaded with shape: {data.shape}")
    return data, clss, lens, meta

# Load trained model
# ------------------
from trained_model.MLP_architecture import MLP_STIM

def load_model(model_path):
    model = MLP_STIM(9, 1).to(device)
    mapped_state_dict = torch.load(os.path.join(model_path,'MLP_parameters.pth'), 
                                weights_only=True, map_location=device)
    model.load_state_dict(mapped_state_dict)
    model.eval()
    print(f"Model loaded from {model_path}")
    return model


# Tabulate results
# ----------------
def tabulate_results(predictions, dts, meta):
    all_predictions = []
    mean_predictions = []
    all_meta = []
    for i, preds in enumerate(predictions):
        mean_predictions.append(np.mean(preds))
        for pred in preds:
            all_predictions.append(pred[0])
            all_meta.append(meta.iloc[i])
    preds = pd.DataFrame(all_meta)
    preds['probability'] = all_predictions

    preds_per_rec = pd.DataFrame(meta)
    preds_per_rec['probability'] = np.array(mean_predictions)
    preds_per_rec['time'] = np.array(dts,dtype=np.float32)

    return preds, preds_per_rec

#%% Main processing functions
# -------------------------
def process_recording(model, recording, fsamp, b, a):
    if fsamp != 24000:
        num_samples = int(recording.shape[0] * 24000 / fsamp)
        artifact_free_data = resample(recording, num_samples)
        fsamp = 24000
    start_time = time.time()
    filtered_data = lib.filter_data(recording, b, a)
    artifact_free_data, _ = lib.remove_artifacts(filtered_data, fsamp)
    if artifact_free_data.shape[0] < fsamp:
        predictions = np.array([])
        dt = time.time() - start_time
    else:
        features = lib.extract_features(artifact_free_data, fsamp)
        features = torch.tensor(features, dtype=torch.float32).to(device)
        with torch.no_grad():
            predictions = model(features)
            predictions = torch.sigmoid(predictions)
        dt = time.time() - start_time
        predictions = predictions.cpu().numpy()
    return predictions, dt

def main(model, raw_data, lens, fsamp, b, a):
    batch_size = 1
    results = []
    num_batches = (raw_data.shape[0] + batch_size - 1) // batch_size
    for i in range(0, raw_data.shape[0], batch_size):
        batch_count = i // batch_size + 1
        print(f"batch (x/{num_batches}): {batch_count}", end='\r')
        batch_results = Parallel(n_jobs=4, timeout=3600)(
            delayed(process_recording)(model, raw_data[j, :lens[j]], fsamp, b, a)
            for j in range(i, min(i + batch_size, raw_data.shape[0])))
        results.extend(batch_results)
    predictions, dts = zip(*results)
    return predictions, dts

#%% 

# Load data and model
# -------------------
model_path = os.path.join(os.getcwd(),'trained_model')
model = load_model(model_path)

filepath = 'data.npz'
metapath = 'metadata.csv'
print('Loading data...')
raw_data, clss, lens, meta = load_data(filepath, metapath)

# Initialize filter coefficients
# ------------------------------
fsamp = 24000 # adjust sampling frequency if needed
b, a = lib.initialize_filter_coefficients(fsamp)

# Initialize empty dataframes that will store the results
# ---------------------------------------------------------
all_preds = pd.DataFrame()              
all_preds_per_rec = pd.DataFrame()

# Process recordings in batches
# -----------------------------
num_recordings = 200 
num_partitions = raw_data.shape[0] // num_recordings + 1

for start_idx in range(0, raw_data.shape[0], num_recordings):
    partition_count = start_idx // num_recordings + 1
    end_idx = min(start_idx + num_recordings, raw_data.shape[0])
    sub_raw_data = raw_data[start_idx:end_idx]
    sub_lens = lens[start_idx:end_idx]
    sub_meta = meta.iloc[start_idx:end_idx]
    print("")
    print(f"Partition {partition_count}/{num_partitions}")
    predictions, dts = main(model, sub_raw_data, sub_lens, fsamp, b, a)
    preds, preds_per_rec = tabulate_results(predictions, dts, sub_meta)
    all_preds = pd.concat([all_preds, preds], ignore_index=True)
    all_preds_per_rec = pd.concat([all_preds_per_rec, preds_per_rec], ignore_index=True)

print()
print('Done!')

# Save the combined results
# -------------------------
all_preds.to_csv('preds.csv', index=False, sep=';')
all_preds_per_rec.to_csv('preds_per_rec.csv', index=False, sep=';')
