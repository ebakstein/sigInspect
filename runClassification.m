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

%% PREPROCESSING - trimming or removing NaN
ratio = 0.5;
%trimmedSignals = trimIncompleteWindows(loadedSignals, ratio);

%[loadedSignalsClean, stats] = removeSignalsWithNaN(loadedSignals, 'window');

%% Covnerting annotations
signals = loadedSignals; %trimmedSignals;
signals = convertLoadedSignalsAnnotations(signals);

%% RUN CLASSIFICATION
disp('--- Starting classification ---');

% PSD method
disp('> Running PSD classification...');
psdThr = 0.01;
[signals, classifErrorsPSD] = classifyLoadedSignals(signals, 'psd', psdThr);

%%
% COV method
disp('> Running COV classification...');
% parameters: [threshold, winLength, aggregPerc]
covParams = {1.2, 0.25, 0.25};
[signals, classifErrorsCOV] = classifyLoadedSignals(signals, 'cov', covParams);

% SVM method
disp('> Running SVM classification...');
[signals, classifErrorsSVM] = classifyLoadedSignals(signals, 'svm');

% SAVE RESULTS
classifiedSignalsPath = fullfile(dataFolder, 'classifiedSignalsCZSKNaNRemovalWindow.mat');
disp(['Saving classified signals to ', classifiedSignalsPath]);
save(classifiedSignalsPath, 'signals', '-v7.3');

disp('--- Classification completed successfully ---');

%%
% Classification
% Optional loading of existing data
% classifiedSignalsPath = fullfile(dataFolder, 'classifiedSignalsCZSKNaNRemoval.mat');
% load(classifiedSignalsPath, 'signals');

% Evaluate classification
[results, signals] = evaluateClassificationMethods(signals);

disp('--- Evaluation completed successfully ---');

%%
function [resultsTable, loadedSignals] = evaluateClassificationMethods(loadedSignals)
    % Configuration
    binaryMethods   = {'psd','cov'};
    svmMethod       = 'svm';
    svmLabels       = {'POW','BASE','FREQ'};
    centerKeys      = {'sig_Bra','sig_Brn','sig_Kos','sig_Olo'};
    centerNames     = {'Bra','Brn','Kos','Olo'};
    overallCenter   = 'All';
    
    % Build method list
    methodNames = [binaryMethods, {svmMethod}];
    signalNames = fieldnames(loadedSignals);
    
    % Helper: detect which center a signal belongs to
    sigCenter = repmat({overallCenter}, numel(signalNames), 1);
    for i = 1:numel(signalNames)
        for c = 1:numel(centerKeys)
            if startsWith(signalNames{i}, centerKeys{c})
                sigCenter{i} = centerNames{c};
                break;
            end
        end
    end
    
    % Accumulator: a dynamic list of rows for the final table
    rows = [];
    
    % Row builder
    function addRow(method, center, label, counts)
        m = countsToMetrics(counts);
        newRow = table( ...
            string(method), string(center), string(label), ...
            counts.TP, counts.FP, counts.TN, counts.FN, ...
            m.accuracy, m.sensitivity, m.specificity, ...
            m.precision, m.recall, m.f1, m.balancedAcc, ...
            'VariableNames', {'Method','Center','Label','TP','FP','TN','FN', ...
                              'Accuracy','Sensitivity','Specificity','Precision', ...
                              'Recall','F1','BalancedAcc'});
        rows = [rows; newRow]; %#ok<AGROW>
    end
    
    % Process each signal
    allCounts = struct();   % aggregated per method/center/label
    
    function acc(method, center, label, TP,FP,TN,FN)
        key = sprintf('%s_%s_%s', method, center, label);
        if ~isfield(allCounts, key)
            allCounts.(key) = struct('TP',0,'FP',0,'TN',0,'FN',0);
        end
        c = allCounts.(key);
        c.TP = c.TP + TP;
        c.FP = c.FP + FP;
        c.TN = c.TN + TN;
        c.FN = c.FN + FN;
        allCounts.(key) = c;
    end
    
    % Loop through signals
    for i = 1:numel(signalNames)
        try
            name = signalNames{i};
            sig  = loadedSignals.(name);
            center = sigCenter{i};
        
            if ~isfield(sig,'classif'), continue; end
        
            yTrueBin   = logical(sig.artifBin);
            yTrueMulti = logical(sig.artifMulti);
        
            % ---------- binary methods ----------
            for m = 1:numel(binaryMethods)
                method = binaryMethods{m};
                if ~isfield(sig.classif,method), continue; end
                pred = logical(sig.classif.(method).annot);
        
                [TP,FP,TN,FN] = getCounts(yTrueBin, pred);
        
                % store into loadedSignals
                loadedSignals.(name).classif.(method).counts = struct('TP',TP,'FP',FP,'TN',TN,'FN',FN);
        
                % accumulate
                labelBinary = 'bin';
                acc(method, center, labelBinary, TP,FP,TN,FN);
                acc(method, overallCenter, labelBinary, TP,FP,TN,FN);
            end
        
            % ---------- SVM multilabel ----------
            if isfield(sig.classif, svmMethod)
                % Skip signals with unknown artifact type (artifMulti row 5)
                if any(yTrueMulti(5,:) == 1)
                    fprintf('Skipping %s (contains unknown artifact in row 5)\n', name);
                    continue;
                end
    
                pred = sig.classif.(svmMethod).annot;  % 1 x N x 3
                for l = 1:numel(svmLabels)
                    yPred = squeeze(pred(1,:,l));
                    yTrue = yTrueMulti(l,:);
                    [TP,FP,TN,FN] = getCounts(yTrue, yPred);
        
                    loadedSignals.(name).classif.(svmMethod).counts(l) = struct('TP',TP,'FP',FP,'TN',TN,'FN',FN);
        
                    label = svmLabels{l};
                    acc(svmMethod, center, label, TP,FP,TN,FN);
                    acc(svmMethod, overallCenter, label, TP,FP,TN,FN);
                end
            end
        catch ME
            msg = ME.message;
            warning('❌ Failed classification for %s: %s', name, msg);
        end
    end

    % Convert accumulated counts into output table
    keys = fieldnames(allCounts);
    for k = 1:numel(keys)
        key = keys{k};
        parts = split(string(key), "_");
        method = parts{1};
        center = parts{2};
        label  = parts{3};
        addRow(method, center, label, allCounts.(key));
    end
    
    resultsTable = rows;
end

% =========================================================================
% Helper Functions
% =========================================================================
function [TP, FP, TN, FN] = getCounts(yTrue, yPred)
    yTrue = logical(yTrue(:));
    yPred = logical(yPred(:));
    
    TP = sum(yTrue &  yPred);
    TN = sum(~yTrue & ~yPred);
    FP = sum(~yTrue &  yPred);
    FN = sum( yTrue & ~yPred);
end

function metrics = countsToMetrics(c)
    TP=c.TP; FP=c.FP; TN=c.TN; FN=c.FN;
    
    metrics.accuracy    = (TP+TN) / max(1, TP+TN+FP+FN);
    metrics.sensitivity = TP / max(1, TP+FN);
    metrics.specificity = TN / max(1, TN+FP);
    metrics.precision   = TP / max(1, TP+FP);
    metrics.recall      = metrics.sensitivity;
    metrics.f1          = 2*TP / max(1, 2*TP+FP+FN);
    metrics.balancedAcc = (metrics.sensitivity + metrics.specificity)/2;
end

