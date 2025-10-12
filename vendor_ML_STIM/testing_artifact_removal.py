"""

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  ML-STIM: Machine Learning for SubThalamic nucleus Intraoperative Mapping  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Python script for processing and classifying MicroElectrode Recordings (MERs). 
This script loads MERs, processes it (applies filters and removes artifacts), 
compares annotations from experts with ML-STIM artifact removal algorithm, 
and evaluates performance (accuracy, sensitivity, specificity, precision, F1).

Author(s): Fabrizio SCISCENTI (fabrizio.sciscenti@polito.it)
           PolitoBIOMed Lab and BIOLAB, Politecnico di Torino, Turin, Italy

           Marco GHISLIERI (marco.ghislieri@polito.it)
           PolitoBIOMed Lab and BIOLAB, Politecnico di Torino, Turin, Italy 

Last Update: 15-05-2025
Extended:  Pavla Mašková (15-09-2025)

Files:
------
    testing_artifact_removal.py
    convert_data_to_mlstim_format.py
    ML_STIM_lib.py
"""

import os
import numpy as np
import pandas as pd
import torch
import time
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
from scipy.signal import resample
from sklearn.metrics import accuracy_score, recall_score, precision_score, f1_score, confusion_matrix
import warnings
import sys
sys.path.append(os.path.join(os.getcwd(), 'vendor_ML_STIM'))
import ML_STIM_lib as lib
warnings.filterwarnings("ignore")

# --------------------------
# Device setup
# --------------------------
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print("Running on GPU!" if device.type == 'cuda' else "GPU not detected. Running on CPU.")

# --------------------------
# Data loading
# --------------------------
def load_data(filepath, metapath, artpath):
    with np.load(filepath) as npfh:
        data = npfh['data']
    meta = pd.read_csv(metapath, sep=';')
    lens = meta['length'].to_numpy()
    clss = meta['class'].to_numpy()
    artifacts = np.load(artpath)
    print(f"Data loaded: signals {data.shape}, artifacts {artifacts.shape}")
    return data, clss, lens, meta, artifacts

# --------------------------
# Utility functions
# --------------------------
def decode_artifacts(annotations, order=None):
    """
    Decode artifact annotations into a boolean indicator.
    
    Parameters
    ----------
    annotations : array-like
        Numeric annotations (0 = clean, >0 = some artifact, or sums of powers of 2).
    order : int or None, optional
        - If int: check specific artifact type 
          (order=1 -> mask=1, order=2 -> mask=2, order=3 -> mask=4, etc).
        - If None: return True if any artifact is present (annotation > 0).
    
    Returns
    -------
    binary : ndarray
        Boolean array of shape (len(annotations),).
    """
    annotations = np.asarray(annotations, dtype=int)

    if order is None:
        return annotations > 0
    else:
        mask = 2**(order-1)
        return (annotations & mask) > 0

def adjust_annotations_to_halfsec(annotations_1s):
    """
    Expand 1-second expert annotations into 0.5-second resolution.
    Removes NaNs before expansion.
    """
    # Drop NaNs at the tail (padding)
    valid = annotations_1s[~np.isnan(annotations_1s)]
    # Convert to binary (any artifact type > 0 → order = None, otherwise checks for specified artifact type by order)
    valid = decode_artifacts(valid, order=1)
    # Expand each 1s label into two 0.5s labels
    return np.repeat(valid, 2)

def align_annotations(expert_halfsec, mlstim_mask):
    """Truncate to common length for comparison."""
    L = min(len(expert_halfsec), len(mlstim_mask))
    return expert_halfsec[:L], mlstim_mask[:L]

def compute_metrics(y_true, y_pred):
    acc = accuracy_score(y_true, y_pred)
    sens = recall_score(y_true, y_pred, zero_division=0)  # sensitivity
    spec = recall_score(y_true, y_pred, pos_label=0, zero_division=0)  # specificity
    prec = precision_score(y_true, y_pred, zero_division=0)
    f1 = f1_score(y_true, y_pred, zero_division=0)
    cm = confusion_matrix(y_true, y_pred)
    return {
        "accuracy": acc,
        "sensitivity": sens,
        "specificity": spec,
        "precision": prec,
        "f1": f1,
        "confusion_matrix": cm.tolist()
    }

# --------------------------
# Artifact processing (ML-STIM)
# --------------------------
def process_recording_artifacts(recording, fsamp, b, a):
    if fsamp != 24000:
        num_samples = int(recording.shape[0] * 24000 / fsamp)
        recording = resample(recording, num_samples)
        fsamp = 24000

    filtered_data = lib.filter_data(recording, b, a)
    _, artifact_mask = lib.remove_artifacts(filtered_data, fsamp)

    # Convert sample-level mask → 0.5s windows (100 ms windows, 50% overlap)
    win_len = int(100e-3 * fsamp)
    step = win_len // 2
    winMask = [int(artifact_mask[start:start+win_len].any())
               for start in range(0, len(artifact_mask), step)]
    winMask = np.array(winMask, dtype=int)

    return winMask

# --------------------------
# Batch processing
# --------------------------
def process_all_recordings(raw_data, lens, artifacts, fsamp, meta,
                           excel_path="artifact_stats.xlsx"):
    metrics_list = []
    all_winMasks = []
    all_annotations = [] 
    
    b, a = lib.initialize_filter_coefficients(fsamp)

    for idx in range(100): #raw_data.shape[0]):
        print(f"Processing recording {idx+1}/{raw_data.shape[0]}", end='\r')
        recording = raw_data[idx, :lens[idx]]
        mlstim_mask = process_recording_artifacts(recording, fsamp, b, a)

        # Expert annotations: already in 1s windows
        expert_vec = artifacts[idx, :lens[idx]]
        if np.all(np.isnan(expert_vec)):
            continue  # skip if no expert annotations

        # Expand from 1s → 0.5s windows
        expert_halfsec = adjust_annotations_to_halfsec(expert_vec)

        # Align with ML-STIM
        y_true, y_pred = align_annotations(expert_halfsec, mlstim_mask)

        # Save per-window annotations
        df_ann = pd.DataFrame({
            "recording": idx,
            "window_idx": np.arange(len(y_true)),
            "true_class": y_true,
            "pred_class": y_pred
        })
        all_annotations.append(df_ann)

        # Metrics
        metrics = compute_metrics(y_true, y_pred)
        metrics["recording"] = idx
        metrics_list.append(metrics)

        all_winMasks.append(mlstim_mask)

    # Combine into DataFrames
    df_metrics = pd.DataFrame(metrics_list)
    df_annotations = pd.concat(all_annotations, ignore_index=True)

    # Summary (mean ± std)
    metrics_only = df_metrics.drop(columns=["confusion_matrix", "recording"])
    summary_mean = metrics_only.mean()
    summary_std = metrics_only.std()
    summary_df = pd.DataFrame({
        "metric": summary_mean.index,
        "mean": summary_mean.values,
        "std": summary_std.values
    })

    # Save everything to one Excel file with sheets
    with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
        df_metrics.to_excel(writer, sheet_name="metrics", index=False)
        df_annotations.to_excel(writer, sheet_name="annotations", index=False)
        summary_df.to_excel(writer, sheet_name="summary", index=False)

    print(f"\nSaved all results to {excel_path}")
    print("\nOverall Performance Summary (mean ± std):")
    for metric, m, s in zip(summary_mean.index, summary_mean.values, summary_std.values):
        print(f"{metric:12s}: {m:.3f} ± {s:.3f}")

    return df_metrics, all_winMasks, summary_df, df_annotations


# --------------------------
# Main
# --------------------------
if __name__ == "__main__":
    filepath = 'data_external_microrecordingCZSK/converted_data/data.npz'
    metapath = 'data_external_microrecordingCZSK/converted_data/metadata.csv'
    artpath = 'data_external_microrecordingCZSK/converted_data/artifacts.npy'

    print('Loading data...')
    raw_data, clss, lens, meta, artifacts = load_data(filepath, metapath, artpath)
    fsamp = 24000

    df_metrics, all_winMasks, summary_df, df_annotations = process_all_recordings(
        raw_data, lens, artifacts, fsamp, meta,
        excel_path="artifact_metrics_ARTIF.xlsx"
    )

