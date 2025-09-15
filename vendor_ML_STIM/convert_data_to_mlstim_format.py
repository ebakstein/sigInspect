"""
This script converts the data from the original format to the ML-STIM format.
It saves transformed data, metadata, and artifacts in the output directory.

Original data:
- Metadata: .mat file with struct 'infoCZSK' (paths, patient name, sampling freq, #samples, #channels, artifacts, …).
- Data:     per-patient .mat file with struct 'd' (2D array: [channels x samples]).

Desired ML-STIM format (described in README.md):
- data.npz  : NumPy array with all signals (zero-padded).
- metadata.csv : Pandas DataFrame with metadata.
- artifacts.npy : NumPy array with all artifacts (NaN-padded).

Edited: 11.09.2025 Pavla Mašková
"""

import os
import numpy as np
import pandas as pd
import scipy.io
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

def load_metadata_info(metadata_path):
    """
    Load the main metadata file (microrecordingCZSK.mat)
    Returns:
        metadata (dict): patient info, paths, sampling, etc.
        artifacts (list): list of artifact matrices (n_channels x n_samples) per patient
    """
    print(f"Loading metadata from: {metadata_path}")
    mat_data = scipy.io.loadmat(metadata_path)
    
    # Extract the infoCZSK struct
    info_struct = mat_data['infoCZSK']
    
    # Get the field names
    field_names = info_struct.dtype.names
    print(f"Available fields: {field_names}")
    print(f"Struct shape: {info_struct.shape}")
    
    # Initialize metadata dict
    metadata = {field: [] for field in ['patient', 'path', 'samplingFreq', 'Nsamples', 'Nchannels']}
    artifacts = []

    # Extract information per patient
    for i in range(info_struct.shape[1]):  # shape[1] = number of patients
        entry = info_struct[0, i]
        
        # Normal metadata 
        for field in metadata.keys():
            try:
                field_value = entry[field]
                if field_value.size > 0:
                    if field in ['path', 'patient']:
                        metadata[field].append(str(field_value[0]))
                    else:
                        metadata[field].append(float(field_value))
                else:
                    metadata[field].append('' if field in ['path', 'patient'] else 0)
            except:
                metadata[field].append('' if field in ['path', 'patient'] else 0)
        
        # Artifacts field
        if "artifacts" in field_names:
            try:
                art_matrix = entry["artifacts"]
                if art_matrix.size > 0:
                    artifacts.append(np.array(art_matrix, dtype=np.float32))
                else:
                    artifacts.append(None)
            except Exception as e:
                print(f"Error reading artifacts for patient {metadata['patient'][-1]}: {e}")
                artifacts.append(None)
        else:
            artifacts.append(None)

    return metadata, artifacts

def process_patient_data(patient_file, art_matrix, 
                         all_recordings, all_artifacts, all_metadata,
                         patient_name, samp_freq, n_samp, n_chan):
    """
    Process all recordings and artifacts from a single patient file.
    """
    print(f"Processing patient data in file: {patient_file}", end='\r')

    try:
        # Load the .mat file
        mat_data = scipy.io.loadmat(str(patient_file))
        if 'd' not in mat_data:
            print(f"No 'd' object found in {patient_file}")
            return

        # Get the data object 'd'
        data_obj = mat_data['d']

        # Verify dimensions match expected values
        if data_obj.shape[0] != n_chan:
            print(f"Warning: Expected {n_chan} channels, found {data_obj.shape[0]}")
        if data_obj.shape[1] != n_samp:
            print(f"Warning: Expected {n_samp} samples, found {data_obj.shape[1]}")

        # Convert artifacts matrix if available
        if art_matrix is not None:
            art_matrix = np.array(art_matrix, dtype=np.float32)
        else:
            art_matrix = None

        # Loop through channels
        for ch in range(data_obj.shape[0]):
            # Signal
            recording = np.nan_to_num(
                data_obj[ch, :], nan=0.0, posinf=0.0, neginf=0.0
            ).astype(np.float32)
            all_recordings.append(recording)

            # Artifact
            if art_matrix is not None and ch < art_matrix.shape[0]:
                artifact_vec = art_matrix[ch, :].astype(np.float32)
                all_artifacts.append(artifact_vec)
            else:
                all_artifacts.append(None)

            # Metadata
            all_metadata.append({
                'patient': patient_name,
                'side': 'UNKNOWN',
                'electrode': f'Channel{ch+1}',
                'depth': 0,
                'length': len(recording),
                'class': 0,
                'source_file': patient_file,
                'channel': ch + 1,
                'sampling_freq': samp_freq,
                'expected_samples': n_samp,
                'expected_channels': n_chan
            })

    except Exception as e:
        print(f"Error processing {patient_file}: {e}")


def convert_to_mlstim_format(metadata_path, base_path, output_dir="converted_data"):
    """
    Main function to convert data to ML-STIM format
    """
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Load main metadata + artifacts
    metadata_info, artifacts_info = load_metadata_info(metadata_path)
    
    # Extract patient paths and other info
    patient_paths = metadata_info.get('path', [])
    patient_names = metadata_info.get('patient', [])
    sampling_freqs = metadata_info.get('samplingFreq', [])
    n_samples = metadata_info.get('Nsamples', [])
    n_channels = metadata_info.get('Nchannels', [])
    
    print(f"Found {len(patient_paths)} patients")
    
    # Initialize collections
    all_recordings, all_artifacts, all_metadata = [], [], []
    
    # Process each patient entry
    for i, (patient_path, patient_name, samp_freq, n_samp, n_chan, art_matrix) in enumerate(zip(patient_paths, patient_names, sampling_freqs, n_samples, n_channels, artifacts_info)):
        print(f"Processing patient {i+1}/{len(patient_paths)}: {patient_name}", end='\r')
        
        # Skip patients whose names start with "Kos_" to avoid memory issues - sampling frequency is float
        if str(patient_name).startswith('Kos_'):
            print(f"Skipping {patient_name} (starts with 'Kos_')", end='\r')
            continue
        
        # Construct full path
        full_path = os.path.join(base_path, patient_path)
        
        if os.path.exists(full_path):
            process_patient_data(full_path, art_matrix, all_recordings, all_artifacts, all_metadata, patient_name, samp_freq, n_samp, n_chan)
        else:
            print(f"Patient path not found: {full_path}")
    
    if not all_recordings:
        print("No recordings found!")
        return
    
    # Find the maximum length of data for zero-padding
    max_length_data = max(len(recording) for recording in all_recordings)
    print(f"\nMaximum recording length: {max_length_data}")
    
    # Zero-pad all recordings to the same length
    data_matrix = np.array([
        np.pad(r, (0, max_length_data - len(r)), constant_values=0)
        if len(r) < max_length_data else r for r in all_recordings
    ], dtype=np.float32)

    # Find the maximum length of data for zero-padding
    max_length_artifacts = max(len(artifact) for artifact in all_artifacts)
    print(f"Maximum recording length: {max_length_artifacts}")

    # Nan-pad all artifacts
    artifacts_matrix = np.array([
        np.pad(a, (0, max_length_artifacts - len(a)), constant_values=np.nan)
        if (a is not None and len(a) < max_length_artifacts)
        else (np.full(max_length_artifacts, np.nan, dtype=np.float32) if a is None else a[:max_length_artifacts])
        for a in all_artifacts
    ], dtype=np.float32)

    # Save
    np.savez(os.path.join(output_dir, "data.npz"), data=data_matrix)
    pd.DataFrame(all_metadata).to_csv(
        os.path.join(output_dir, "metadata.csv"), index=False, sep=';'
    )
    np.save(os.path.join(output_dir, "artifacts.npy"), artifacts_matrix)

    print("\nConversion Summary:")
    print(f"  Data shape      : {data_matrix.shape}")
    print(f"  Artifacts shape : {artifacts_matrix.shape}")
    print(f"  Metadata entries: {len(all_metadata)}")
    print(f"  Output dir      : {output_dir}")

    return data_matrix, pd.DataFrame(all_metadata), artifacts_matrix

if __name__ == "__main__":
    # Specify the path to the data and metadata file
    base_path = "data_external_microrecordingCZSK"
    metadata_file = os.path.join(base_path, "microrecordingCZSK.mat") 
    # Specify the output directory
    output_dir = os.path.join(base_path, "converted_data")
    # Convert the data to ML-STIM format
    data, metadata, artifacts = convert_to_mlstim_format(metadata_file, base_path, output_dir)
