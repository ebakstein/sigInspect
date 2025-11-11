clc; clear;

sigInspectAddpath;

%% Parameters
% Define paths
dataFolder = 'data_external_microrecordingCZSK/';
% Metadata saved in mat file
matFile = 'microrecordingCZSK.mat';
matFilePath = fullfile(dataFolder, matFile);
loadedSignalsPath = fullfile(dataFolder, 'loadedSignalsCZSK.mat');

%% LOAD SIGNALS
disp('--- Loading signals and annotations ---');
% Load or create loadedSignals
if isfile(loadedSignalsPath)
    fprintf('Loaded signals file exists. Loading data...\n');
    load(loadedSignalsPath, 'loadedSignals');
else
    fprintf('Loaded signals file not found. Running `loadSignalsFromMat`...\n');
    loadedSignals = loadSignalsFromMat(matFilePath, true, loadedSignalsPath);
end

%% RUN CLASSIFICATION
disp('--- Starting classification ---');

% PSD method
disp('> Running PSD classification...');
loadedSignals = classifyLoadedSignals(loadedSignals, 'psd');

% COV method
disp('> Running COV classification...');
% parameters: [threshold, winLength, aggregPerc]
covParams = {1.2, 0.25, 0.25};
loadedSignals = classifyLoadedSignals(loadedSignals, 'cov', covParams);

% SVM method
disp('> Running SVM classification...');
loadedSignals = classifyLoadedSignals(loadedSignals, 'svm');


%% SAVE RESULTS
classifiedSignalsPath = fullfile(dataFolder, 'classifiedSignalsCZSK.mat');
disp(['Saving classified signals to ', classifiedSignalsPath]);
save(classifiedSignalsPath, 'loadedSignals', '-v7.3');

disp('--- Classification completed successfully ---');

%%

results = evaluateClassification(loadedSignals, 'psd')


function results = evaluateClassification(loadedSignals, method)
    signalNames = fieldnames(loadedSignals);
    results = struct();
    
    for i = 1:numel(signalNames)
        sigName = signalNames{i};
        sigData = loadedSignals.(sigName);
    
        if ~isfield(sigData.classif, method) || isempty(sigData.artif)
            continue;
        end
    
        ref = logical(sigData.artif);  % reference annotations
        pred = logical(sigData.classif.(method).annot);  % predicted
        
        % ensure same length
        L = min(numel(ref), numel(pred));
        ref = ref(1:L);
        pred = pred(1:L);
    
        tp = sum(ref & pred);
        tn = sum(~ref & ~pred);
        fp = sum(~ref & pred);
        fn = sum(ref & ~pred);
    
        acc = (tp + tn) / (tp + tn + fp + fn);
        sens = tp / (tp + fn);
        spec = tn / (tn + fp);
    
        results.(sigName) = struct('Accuracy', acc, 'Sensitivity', sens, 'Specificity', spec);
    end
end


%%
% Define the file
cacheFile = 'featuresCache.mat';

% Load all variables into a struct
S = load(cacheFile);

% Get field names
fields = fieldnames(S);

% Pattern: x_ followed by one or more digits
pattern = '^x_\d+$';

% Find matching field names
toRemove = fields(~cellfun('isempty', regexp(fields, pattern)));

% Remove those fields
for i = 1:numel(toRemove)
    S = rmfield(S, toRemove{i});
end

% Save back (overwrite existing file)
save(cacheFile, '-struct', 'S');

fprintf('Removed %d fields matching "%s" and saved back to %s.\n', numel(toRemove), pattern, cacheFile);
