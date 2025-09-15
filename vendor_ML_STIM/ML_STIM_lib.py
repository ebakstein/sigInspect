"""
This module provides functions for signal processing, including filtering, 
artifact removal, and feature extraction from MicroElectrode Recordings (MERs).

Author(s): Fabrizio SCISCENTI (fabrizio.sciscenti@polito.it)
           PolitoBIOMed Lab and BIOLAB, Politecnico di Torino, Turin, Italy

           Marco GHISLIERI (marco.ghislieri@polito.it)
           PolitoBIOMed Lab and BIOLAB, Politecnico di Torino, Turin, Italy 

Last Update: 02-04-2025

Functions:
----------
    initialize_filter_coefficients(fsamp)
        Initializes filter coefficients for highpass, lowpass, and notch filters.
    filter_data(data, b, a)
        Applies a series of filters to the input data.
    remove_artifact(data, fsamp)
        Removes artifacts from the input data based on variance thresholding.
    extract_features(data, fsamp)
        Extracts features from data segments

"""

from scipy import signal
import numpy as np
from scipy.signal import lfilter, welch
from scipy.stats import skew, kurtosis

def initialize_filter_coefficients(fsamp):
    """
    Initializes filter coefficients for highpass, lowpass, and notch filters.
    Parameters:
        fsamp (float): Sampling frequency of the signal.
    Returns:
        tuple: Two dictionaries containing the filter coefficients (b, a).
    """
    b = {}
    a = {}
    nyquist = 0.5 * fsamp

    f_band = (250, 5000)
    b['band'], a['band'] = signal.butter(6, [f_band[0] / nyquist, f_band[1] / nyquist], btype='bandpass')

    notch_freqs = np.arange(200, 5000, 50)
    for freq in notch_freqs:
        b[f'notch_{freq}Hz'], a[f'notch_{freq}Hz'] = signal.iirnotch(freq, 50, fsamp)
    
    return b, a

def filter_data(data, b, a):
    """
    Applies a series of filters to the input data.
    Parameters:
        data (numpy.ndarray): Input MER.
        b (dict): Dictionary of numerator filter coefficients.
        a (dict): Dictionary of denominator filter coefficients.
    Returns:
        numpy.ndarray: Filtered MER.
    """
    filtered_data = data
    for key in b:
        filtered_data = lfilter(b[key], a[key], filtered_data).astype(np.float32)
    return filtered_data

def remove_artifacts(data, fsamp):
    """
    Removes artifacts from the input data based on variance thresholding.
    Parameters:
        data (numpy.ndarray): Input MER.
        fsamp (float): Sampling frequency of the signal.
    Returns:
        tuple: Cleaned MER and artifact mask.
    """
    L = data.shape[0]
    params = {'w': int(100e-3 * fsamp), 'threshold': 1.33}
    artMask = np.zeros(L, dtype=bool)
    Wbis = 0
    t = 0
    start = 0
    stop = params['w']
    pvar = 1e14
    first_it = True    
    s = params['w'] // 2
    w = params['w']
    while start < L:
        epoch = data[start:stop]
        nvar = np.var(epoch)
        if nvar >= params['threshold'] * pvar:
            W = [1,1]
        else:
            W = [0,0]
            pvar = nvar
        if first_it:
            Wbis = W[0]
            first_it = False

        if stop - start < s:
            s = stop - start

        cls = np.mean([W[0], Wbis])
        artMask[t:t + s] = cls >= 0.5

        Wbis = W[1]

        t += s
        start += s
        stop = start + w
        if stop > L:
            stop = L
            w = stop - start

    clearSignal = data[~artMask]
    return clearSignal, artMask

def extract_features(data, fsamp):
    """
    Extracts features from 1-second long data segments.
    Parameters:
        data (numpy.ndarray): Input MER.
        fsamp (float): Sampling frequency of the signal.
    Returns:
        numpy.ndarray: Extracted and normalized features.
    """

    segment_length = fsamp # 1 second
    all_features = []

    # normalization intervals computed from the training set
    max_vals = [8.707813,0.016264,0.073682,0.133387,0.16357,2769,0.435164,0.936993,2.680869]
    min_vals = [1.357211,0.002564,0.019996,-0.250695,0.041139,1644.88,0.153334,0.833393,-0.096321]

    overlap = 0.5
    step = int(segment_length * (1 - overlap))

    for start in range(0, data.shape[0] - segment_length + 1, step):
        segment = data[start:start + segment_length]

        # PSD calculation
        freqs, Pxx = welch(segment, fsamp, nperseg=fsamp // 10)
        P_tot = np.sum(Pxx)
        PSD = Pxx / P_tot
        freqsA, PxxA = welch(np.abs(segment), fsamp, nperseg=fsamp // 10)
        P_totA = np.sum(PxxA)
        PSDA = PxxA / P_totA
       
        # Feature extraction
        avgAbsDiff = np.mean(np.abs(segment - np.mean(segment)))
        pr_8_13Hz = np.sum(PSDA[(freqsA >= 8) & (freqsA <= 13)])
        pr_30_70Hz = np.sum(PSDA[(freqsA >= 30) & (freqsA <= 70)])
        skewness = skew(segment)
        PSDratio = np.sum(Pxx[(freqs >= 3000) & (freqs <= 4000)]) / np.sum(Pxx[(freqs >= 1600) & (freqs <= 2200)])
        zc_count = ((segment[:-1] * segment[1:]) < 0).sum()
        pr_1_2kHz = np.sum(PSD[(freqs >= 1000) & (freqs <= 2000)])
        PSDindex = 0
        for j in range(22):
            PSDindex += np.sum(PSD[(freqs >= 255 + j * 100) & (freqs <= 245 + (j + 1) * 100)])   
        kurt = kurtosis(segment)

        features = np.array([avgAbsDiff, pr_8_13Hz, pr_30_70Hz, skewness, PSDratio, zc_count, pr_1_2kHz, PSDindex, kurt])

        all_features.append(features)

    all_features = np.stack(all_features)

    for i in range(all_features.shape[1]):
        all_features[:, i] = (all_features[:, i] - min_vals[i]) / (max_vals[i] - min_vals[i])
    all_features = np.clip(all_features, 0, 1)
    return all_features