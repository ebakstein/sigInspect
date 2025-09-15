# ML-STIM: Machine Learning for SubThalamic nucleus Intraoperative Mapping

<p align="center">

<img  src="https://github.com/Biolab-PoliTO/ML-STIM/blob/main/docs/ML-STIM_cover.jpg" style="width:100%; height:auto;"/></p>

Deep Brain Stimulation (DBS) of the SubThalamic Nucleus (STN) is an effective electroceutical therapy for treating motor symptoms in patients with Parkinson’s Disease (PD). Accurate placement of the stimulating electrode within the STN is essential for achieving optimal therapeutic outcomes. To this end, MicroElectrode Recordings (MERs) are acquired during surgery to provide intraoperative visual and auditory confirmation of the electrode position. MERs are traditionally analyzed by trained operators. However, this approach is time-consuming and subject to variability. <br>

This work introduces ```ML-STIM```, a machine learning-based pipeline for real-time classification of MERs to identify the STN during DBS procedures. ```ML-STIM``` is designed to ensure high classification accuracy and real-time applicability, incorporating interpretable machine learning techniques to ensure compatibility with clinical practices.

## What ```ML-STIM``` algorithm does:
1.	Loads `numpy` arrays (`.npz`) storing MERs as rows and the relative `.csv` metafile
2.	Applies `ML-STIM` pipeline to each MER:
	-	Filters and removes artifacts from raw MERs
	-	Extract temporal and spectral features from cleared signals
	-	Classifies MERs as inside-STN or outside-STN 
3. 	Exports results in ```.csv``` format.

## Files description:
The following files are provided within this GitHub repository:
- `testing.py`: main script to test `ML-STIM` on the input dataset
- `lib.py`: is a collection of functions called within `testing.py` to perform pre-processing (filtering + artifact removal) and feature extraction.
- `trained_model`: the folder includes the trained model:
	- `MLP_architecture.py` defines a MultiLayer Perceptron (MLP) with predefined architecture
	- `MLP_parameters.pth`containes the trained parameters for the model
 - `data.npz`: example data file containing MERs from a representative subject's hemisphere
 - `metadata.csv`: infos about MERs in `data.npz`.
</p>

## How to prepare your data:
To use this analysis framework, your data must be structured in ```.npz``` data file and a ```.csv``` metafile.
- **Data** : your data file must contain MERs as rows of a NxM matrix where N is the number of recordings, M is the length of the longest recording.
Recordings shorter than M must be zero-padded to length M. The actual recording length must be reported in the metadata file described in the following.
- **Metadata** : Your metadata file should contain a table with N rows and variables (columns):
	- *patient*: patient id (e.g., `P7`)
	- *side*: hemisphere (`LEFT` or `RIGHT`)
	- *electrode*: recording electrode (e.g., `Electrode1`)
	- *depth*: Estimated Distance from Target (EDT) expressed in *μm*
	- *length*: signal length (in samples) before zero-padding
	- *class*: label (`0` for **outside the STN**, `1` for **inside the STN**)

```
	patient	side	electrode	depth	length	class
0	P07	RIGHT	Electrode1	-8000	240000	0
1	P07	RIGHT	Electrode1	-1000	240000	1
2	P07	RIGHT	Electrode1	 2000	227648	0
```
For a representative example of the expected input format, refer to the ```metadata.csv``` and ```data.npz``` file.

</p>

## A simple workflow
A simplified workflow for MER processing and classification looks as follows.

1. **Data loading:**

```r
import numpy as np
import pandas as pd

# Path definition
filepath = "path/to/data.npz"
metapath = "path/to/metadata.csv"

with np.load(filepath) as npfh:
	raw_data = npfh['data']		# Load data matrix
meta = pd.read_csv(metapath)		# Load metadata
```
2. **Signal processing:**

```r
# Import library
import lib

# process signal
fsamp = 24000		# Sampling frequency (Hz)
b, a = lib.initialize_filter_coefficients(fsamp)

recording = raw_data[0, :meta['length'][0].to_numpy()]	# Select a raw signal
filtered_rec = lib.filter_data(recording, b, a)		# Apply filters (band-pass + notch)
artifact_free_data, art_mask = lib.remove_artifact(filtered_data, fsamp)	# Remove artifacts
```
Here's an example of artifact segmentation:
<img  src="https://github.com/Biolab-PoliTO/ML-STIM/blob/main/docs/artifact_segmentation.png" style="width:100%; height:auto;"/> </p>

3. **Feature extraction:**
```r
features = lib.extract_segment_features(artifact_free_data, fsamp)	# from 1-second segments
```

4. **Classification:**
```r
import torch
# Import architecture
from trained_model.MLP_architecture import MLP_STIM

# Define model Path
model_path = 'path/to/trained_model/folder'

# Import model
model = MLP_STIM.to(device)						# Initialize empty model
params = torch.load(os.path.join(model_path,'MLP_parameters.pth'), 	# Import weight and biases
                    weights_only=True, 
		    map_location=torch.device("cuda" if torch.cuda.is_available() else "cpu"))
model.load_state_dict(params)						# Fill MLP architecture

# Classify recordings
prediction = model(features)

# Apply the sigmoid to get the probability of being inside the STN.
prediction = torch.sigmoid(prediction)

# Apply thresholding to get the class label.
class = (prediction >= .51).astype(int)
```

## References
[1] Sciscenti, F., Agostini, V., Rizzi, L., Lanotte, M., & Ghislieri, M. (2025). ML-STIM: Machine Learning for SubThalamic nucleus Intraoperative Mapping (1.0.0) [Paper]. Journal of Neural Engineering DOI: https://doi.org/10.1088/1741-2552/adf579.

[2] Sciscenti, F., Agostini, V., Rizzi, L., Lanotte, M., & Ghislieri, M. (2025). ML-STIM: Machine Learning for SubThalamic nucleus Intraoperative Mapping (1.0.0) [Dataset]. Zenodo. DOI: https://doi.org/10.5281/zenodo.14894226.

##  How to contribute to ```ML-STIM```
Contributions are the heart of the open-source community, making it a fantastic space for learning, inspiration, and innovation. While we've done our best, our code may contain inaccuracies or might not fully meet your needs. If you come across any issues—or have ideas for improvements—we encourage you to contribute! Follow the instructions below to suggest edits or enhancements. Every contribution is **greatly appreciated**!

Bugs are tracked as **GitHub issues**. Whenever you report an issue, please make sure to:
1.	Use a concise and descriptive title
2.	Report your Python version
3.	Report whether the code ran successfully on the test data available within the repository.


## Contacts
**Fabrizio Sciscenti**, Ph.D. Candidate - [BIOLAB@Polito](https://biolab.polito.it/people/fabrizio-sciscenti/) <br>
[@FSciscenti](https://x.com/FSciscenti) - fabrizio.sciscenti@polito.it

**Marco Ghislieri**, Ph.D. - [BIOLAB@Polito](https://biolab.polito.it/people/marco-ghislieri/) <br>
[@MarcoGhislieri](https://twitter.com/MarcoGhislieri) - marco.ghislieri@polito.it
