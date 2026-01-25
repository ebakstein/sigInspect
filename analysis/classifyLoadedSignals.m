function [loadedSignals, classifErrors] = classifyLoadedSignals(loadedSignals, method, varargin)
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

classifErrors = struct();
signalNames = fieldnames(loadedSignals);

fprintf('Starting classification of %d signals using method "%s"...\n', numel(signalNames), method);

for i = 1:numel(signalNames)
    sigName = signalNames{i};
    sigData = loadedSignals.(sigName);

    if ~isfield(sigData, 'data') || isempty(sigData.data)
        msg = 'No signal data.';
        warning('Skipping %s — %s', sigName, msg);
        classifErrors.(method).(sigName) = msg;
    end


    % % Skip if already computed
    % if isfield(sigData, 'classif') && ...
    %    isfield(sigData.classif, method) && ...
    %    isfield(sigData.classif.(method), 'annot')
    % 
    %     fprintf('Annotation known for %s, skipping\n', sigName);
    %     continue;
    % end

    fs = sigData.samplingFreq;
    signal = sigData.data;
    signalId = sigData.sigId;

    % Ensure signal, fs is double for classification
    if ~isa(signal,'double'), signal = double(signal); end
    if ~isa(fs,'double'), fs = double(fs); end

    try
        [annot, featVals, scores] = sigInspectClassify(signal, signalId, fs, method, varargin{:});
        loadedSignals.(sigName).classif.(method).annot = annot;
        loadedSignals.(sigName).classif.(method).featVals = featVals;
        % Only save scores if they are not empty (primarily for SVM)
        if ~isempty(scores)
            loadedSignals.(sigName).classif.(method).scores = scores;
        end
        fprintf('✓ Classified %s\n', sigName);
    catch ME
        msg = ME.message;
        warning('❌ Failed classification for %s: %s', sigName, msg);
        classifErrors.(method).(sigName) = msg;
    end
end

fprintf('Classification finished for all signals.\n');
end