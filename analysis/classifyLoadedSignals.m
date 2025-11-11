function loadedSignals = classifyLoadedSignals(loadedSignals, method, varargin)
% Runs artifact classification on all signals in loadedSignals
% and adds results into each signal structure.
%
% Inputs:
%   loadedSignals - structure from loadSignalsFromMat
%   method        - classification method, e.g. 'psd', 'tree', 'svm'
%   varargin      - additional classifier parameters
%
% Output:
%   loadedSignals - same struct with new field .classif.(method)

signalNames = fieldnames(loadedSignals);

fprintf('Starting classification of %d signals using method "%s"...\n', numel(signalNames), method);

for i = 1:numel(signalNames)
    sigName = signalNames{i};
    sigData = loadedSignals.(sigName);

    if ~isfield(sigData, 'data') || isempty(sigData.data)
        warning('Skipping %s — no signal data.', sigName);
        continue;
    end

    fs = sigData.samplingFreq;
    signal = sigData.data;
    signalId = sigData.sigId;

    % Ensure signal, fs is double for classification
    if ~isa(signal,'double'), signal = double(signal); end
    if ~isa(fs,'double'), fs = double(fs); end

    try
        annot = sigInspectClassify(signal, signalId, fs, method, varargin{:});
        loadedSignals.(sigName).classif.(method).annot = annot;
        fprintf('✓ Classified %s\n', sigName);
    catch ME
        warning('❌ Failed classification for %s: %s', sigName, ME.message);
    end
end

fprintf('Classification finished for all signals.\n');
end