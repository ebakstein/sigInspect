function [loadedSignals] = loadSignalsFromMat(matFilePath, saveData, savePath)
% Loads signals and annotations from .mat files with metadata.
%
% Inputs:
%   matFilePath - path to the metadata .mat file containing 'infoCZSK' structure
%   saveData    - (optional) boolean, save the result (default = false)
%   savePath    - (optional) path for saving output (default = '../loadedSignalsCZSK.mat')
%
% Output:
%   loadedSignals - structure with loaded and processed signals

% Default parameters
if nargin < 2, saveData = false; end
[parentDir, ~, ~] = fileparts(matFilePath);
if nargin < 3
    [grandParentDir, ~, ~] = fileparts(parentDir);
    savePath = fullfile(grandParentDir, 'loadedSignalsCZSK.mat');
end

fprintf('--- Loading signals and annotations ---\n');
fprintf('Metadata file: %s\n', matFilePath);
if saveData
    fprintf('Output will be saved to: %s\n', savePath);
else
    fprintf('Output will not be saved.\n');
end

% Load metadata
matData = load(matFilePath);
if ~isfield(matData, 'infoCZSK')
    error('The provided MAT file does not contain the "infoCZSK" structure.');
end
infoCZSK = matData.infoCZSK;
numSignals = numel(infoCZSK);
fprintf('Found %d signals in metadata.\n', numSignals);

% Initialize output structure
loadedSignals = struct();

progressMarks = round([0.25, 0.5, 0.75, 1] * numSignals);

for i = 1:numSignals
    meta = infoCZSK(i);
    sigId = meta.sigId;
    signalPath = fullfile(parentDir, meta.path);
    fprintf('Processing %d/%d: %s\n', i, numSignals, sigId);

    if ~isfile(signalPath)
        warning('Missing file: %s', signalPath);
        continue;
    end

    % Load actual signal .mat file
    sigData = load(signalPath);
    dataFields = fieldnames(sigData);
    if isempty(dataFields)
        warning('Empty .mat file for signal: %s', sigId);
        continue;
    end
    signalMatrix = double(sigData.(dataFields{1}));
    artifactsMatrix = double(meta.artifacts);

    % Check sampling frequency
    if ~isfield(meta, 'samplingFreq') || isempty(meta.samplingFreq)
        warning('Missing sampling frequency for %s, setting to 24000 Hz.', sigId);
        fs = 24000;
    else
        fs = meta.samplingFreq;
        if abs(fs - round(fs)) > 1e-6
            fprintf('Non-integer sampling frequency detected (%.4f). Rounding to 24000 Hz.\n', fs);
            fs = 24000;
        end
    end

    % Compute samples per 1-second window
    samplesPerWin = round(fs); % since fs is Hz
    totalSamples = size(signalMatrix, 2);
    fullWindows = round(double(totalSamples) / double(samplesPerWin));
    keepSamples = fullWindows * samplesPerWin;

    % Trim incomplete last window if there is less then half samples in the
    % window
    if keepSamples < totalSamples
        signalMatrix = signalMatrix(:, 1:keepSamples);
        artifactsMatrix = artifactsMatrix(:, 1:fullWindows);
        percTrimmed = 100.0 * (double(totalSamples-keepSamples) / double(totalSamples));
        fprintf('Trimmed last incomplete window for %s (%.6f%% of data removed).\n', ...
            sigId, percTrimmed);
        meta.Nsamples = keepSamples;

    end

    % --- Store results ---
    for ch = 1:size(signalMatrix, 1)
        chName = matlab.lang.makeValidName(sprintf('sig_%s_ch%d', sigId, ch));
        loadedSignals.(chName).data = signalMatrix(ch, :);
        loadedSignals.(chName).artif = artifactsMatrix(ch, :);
        loadedSignals.(chName).sigId = sigId;
        loadedSignals.(chName).channelNumber = ch;
        loadedSignals.(chName).samplingFreq = fs;

        % Copy metadata fields if present
        extraFields = {'Nsamples', 'Nchannels', 'info', 'center', ...
                       'artifactAuthors', 'patient'};
        for f = 1:numel(extraFields)
            fld = extraFields{f};
            if isfield(meta, fld)
                loadedSignals.(chName).(fld) = meta.(fld);
            end
        end

        % Placeholder for classification results
        loadedSignals.(chName).classif = struct();
    end

    % Progress updates
    if ismember(i, progressMarks)
        fprintf('Progress: %.0f%%\n', 100 * i / numSignals);
    end
end

% Save output if requested
if saveData
    fprintf('Saving processed data to %s...\n', savePath);
    save(savePath, 'loadedSignals', '-v7.3');
end
fprintf('--- Done. %d signals processed successfully. ---\n', numSignals);

end
