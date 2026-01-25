%% SIGNAL PROCESSING & CLASSIFICATION PIPELINE
% Script for loading Prague + external (CZ/SK) data, classifaction using
% PSD, COV and SVM method using sigInspectClassify function (the same as in
% sigInspect) and for evaluation
% Author: Pavla Mašková
% Date: 11 25

clc; clear; close all;

%%
sigInspectAddpath;

%% CENTRALIZED CONFIGURATION
% Change settings here 

cfg = struct();
% --- Paths ---
cfg.paths.saveDir    = 'data_merged/';
cfg.paths.merged     = fullfile(cfg.paths.saveDir, 'loadedSignals.mat'); % unprocessed loaded signals from Prague plus external data
cfg.paths.classified = fullfile(cfg.paths.saveDir, 'classifiedSignals.mat'); % unprocessed classified signals using all methods
cfg.paths.final = fullfile(cfg.paths.saveDir, 'finalSignals.mat'); % final classified and proccessed data 

% --- External Data (CZ/SK) ---
cfg.ext.folder    = 'data_external_microrecordingCZSK/';
cfg.ext.matFile   = fullfile(cfg.ext.folder, 'microrecordingCZSK.mat');
cfg.ext.cacheFile = fullfile(cfg.ext.folder, 'loadedSignalsCZSK.mat');

% --- Training Data (Prague) ---
cfg.train.folder    = 'data/';
cfg.train.csv       = '_metadataMER2020.csv';
cfg.train.sigFolder = 'signals/';
cfg.train.cache     = fullfile(cfg.train.folder, 'loadedSignalsPrague.mat');
cfg.train.fs        = 24000;

% --- Parameters & Preprocessing ---
cfg.split.seed   = 3; 
cfg.split.ratios = struct('train', 0.65, 'val', 0.2, 'test', 0.15);
cfg.proc.trimRatio        = NaN;      % Set to value (e.g. 0.8) to trim
cfg.proc.NaNRemovalMethod = 'signal'; % Options: 'signal', 'window', or '' (how to handle signals containing NaN)

% --- Classification Hyperparameters ---
cfg.params.psdThr    = 0.012;
cfg.params.covParams = {1.2, 0.25, 0.25}; % [threshold, winLength, aggregPerc]

%%
cfg.params.psdThr    = 0.01;
fprintf('>>> Running PSD (%f)...\n', cfg.params.psdThr);
[signals, ~] = classifyLoadedSignals(signals, 'psd', cfg.params.psdThr);


%% MAIN EXECUTION LOGIC
% At first loading data -> creating mat file from Prague and from CZSK
% folder separately -> then merging and classifying unprocessed ->
% classification made later - user can then experiment with the trimming
% signals or removing NaN using different methods...

% STEP 1: LOAD DATA
if isfile(cfg.paths.classified)
    % A. Priority 1: Load already classified data
    fprintf('>>> SUCCESS: Found classified data - not processed. Loading: %s\n', cfg.paths.classified);
    load(cfg.paths.classified, 'signals');

elseif isfile(cfg.paths.merged)
    % B. Priority 2: Load merged signals and start classification
    fprintf('>>> Merged file found. Proceeding to classification...\n');
    load(cfg.paths.merged, 'signals');
    signals = runPipelineClassification(signals, cfg);

else
    % C. Priority 3: No processed files found. Load from scratch, merge, then classify.
    fprintf('>>> No processed files found. Starting full load/merge/classify sequence...\n');
    signals = loadAndMergeSources(cfg);
    save(cfg.paths.merged, 'signals', '-v7.3');
    
    signals = runPipelineClassification(signals, cfg);
end

%%
% STEP 2: POSTPROCESSING (Trimming / NaN Handling / converting annotation / aggragation of multilabel SVM)
signals = applyPostprocessing(signals, cfg);

%%
% STEP 3: EVALUATION
fprintf('--- Finalizing Results ---\n');
[resultsTable, signals] = evaluateClassificationMethods(signals);

%%
% STEP 4: STATISTICS & VISUALS
plotMetricsByCenter(signals, 'psd');
statsTable = generateParticipantStatistics(signals);

disp('--- Pipeline complete. Results table summary: ---');
disp(head(statsTable));

%%
% STEP 5: SAVE FINAL MAT
save(cfg.paths.final, 'signals', '-v7.3');


%% COMPARISON NEW VS OLD COMPUTATION
% Load OLD fatures using script runSVMClassification
% X = featureSetUndersampled.X;        % Cell array
% Y = featureSetUndersampled.Y;        % Cell array
% signalIds = featureSetUndersampled.signalIds; % Cell array
% featNames = featureSetUndersampled.featNames; % Cell array

% Filter out only Prague signals
origSignals = signals;
% [signals, stats] = filterSignalsByCenter(signals, 'Prague');
% disp(stats)

% Identify Missing Signals
structSignals = fieldnames(signals); % From your first dataset
cellSignals = signalIds;             % From your second dataset

% Find signals present in Struct but missing in Cell
missingInCell = setdiff(structSignals, cellSignals);
% Find common signals to compare
commonSignals = intersect(structSignals, cellSignals);

fprintf('Analysis Results:\n');
fprintf('-----------------\n');
fprintf('Signals in Struct: %d\n', numel(structSignals));
fprintf('Signals in Cell:   %d\n', numel(cellSignals));
fprintf('Signals missing from second dataset: %d\n', numel(missingInCell));
fprintf('Signals common: %d\n', numel(commonSignals));

% Analysis Results - Prague dataset:
% -----------------
% Signals in Struct: 5342
% Signals in Cell:   3565
% Signals missing from second dataset: 1777
% Signals common: 3565

% Remove signals
matchingIndices = find(ismember(structSignals, commonSignals));
filteredSignalNames = structSignals(matchingIndices);
removeNames = setdiff(structSignals, filteredSignalNames);
filteredSignals = rmfield(signals, removeNames);

% Join back the nonPrague signals
isOtherCenter = ~cellfun(@isempty, regexp(structSignals, '^sig_[^0-9]'));
otherNames = structSignals(isOtherCenter);
otherCentersStruct = rmfield(origSignals, structSignals(~isOtherCenter));
finalSignals = cell2struct([struct2cell(filteredSignals); struct2cell(otherCentersStruct)], ...
                           [fieldnames(filteredSignals); fieldnames(otherCentersStruct)], 1);

fprintf('Added %d non-Prague signals based on naming convention.\n', numel(otherNames));
signals = finalSignals;
clear finalSignals filteredSignals otherCentersStruct 
% Then run STEP 2 postprocessing,...

%% HELPER FUNCTIONS 

function signals = loadAndMergeSources(cfg)
    % Handles the loading of External and Prague datasets
    fprintf('Loading External Source...\n');
    if isfile(cfg.ext.cacheFile)
        tmp = load(cfg.ext.cacheFile);
        signalsExt = tmp.loadedSignals;
    else
        signalsExt = loadSignalsFromMat(cfg.ext.matFile, true, cfg.ext.cacheFile);
    end
    
    fprintf('Loading Prague Source...\n');
    if isfile(cfg.train.cache)
        tmp = load(cfg.train.cache);
        rawPrague = tmp.loadedSignals;
    else
        rawPrague = loadSignalsPrague(cfg.train.folder, cfg.train.csv, cfg.train.sigFolder, true, cfg.train.cache);
    end
    
    signalsPrague = upgradePragueStruct(rawPrague, cfg.train.fs);
    signalsPrague = applyPragueLabels(signalsPrague, cfg.split);
    signals = mergeSignalStructs(signalsExt, signalsPrague);
end

function signals = runPipelineClassification(signals, cfg)
    % Takes some time, especially svm...
    % Sequentially runs PSD, COV, and SVM
    fprintf('>>> Converting annotations...\n');
    signals = convertLoadedSignalsAnnotations(signals);
    
    fprintf('>>> Running PSD (%f)...\n', cfg.params.psdThr);
    [signals, ~] = classifyLoadedSignals(signals, 'psd', cfg.params.psdThr);
    
    fprintf('>>> Running COV...\n');
    [signals, ~] = classifyLoadedSignals(signals, 'cov', cfg.params.covParams);
    
    fprintf('>>> Running SVM...\n');
    [signals, ~] = classifyLoadedSignals(signals, 'svm');
    
    % Save checkpoint
    save(cfg.paths.classified, 'signals', '-v7.3');
end

function signals = applyPostprocessing(signals, cfg)
% Trim signals, remove Nan and unknown artifact type windows/signals,
% convert annotations to binary format, compute aggregation of SVM results
    % Trimming and NaN removal
    if ~isnan(cfg.proc.trimRatio)
        fprintf('> Trimming incomplete windows...\n');
        signals = trimIncompleteWindows(signals, cfg.proc.trimRatio);
    end
    if ~isempty(cfg.proc.NaNRemovalMethod)
        fprintf('> Removing NaNs and unknown artifact types (Method: %s)...\n', cfg.proc.NaNRemovalMethod);
        winSec = 1;
        [signals, stats] = cleanSignals(signals, cfg.proc.NaNRemovalMethod, winSec);
        disp(stats)
    end

    fprintf('Applying Prague Train/Val/Test split (Seed: %d)...\n', cfg.split.seed);
    signals = applyPragueLabels(signals, cfg.split);
    signals = convertLoadedSignalsAnnotations(signals);
    signals = computeSVMUnion(signals);
end

function signals = computeSVMUnion(signals)
    % Logical aggregation of SVM scores
    sigNames = fieldnames(signals);
    for i = 1:numel(sigNames)
        name = sigNames{i};
        if isfield(signals.(name), 'classif') && isfield(signals.(name).classif, 'svm')
            svmRes = signals.(name).classif.svm;
            if isfield(svmRes, 'annot') && ~isempty(svmRes.annot)
                signals.(name).classif.svmunion.annot = max(svmRes.annot, [], 3);
                signals.(name).classif.svmunion.scores = max(svmRes.scores, [], 3);
            end
        end
    end
end

function signalsOut = upgradePragueStruct(rawStruct, fs)
    % Only adds missing fields. Does NOT split data.
    signalsOut = struct();
    names = fieldnames(rawStruct);
    
    for i = 1:numel(names)
        fName = names{i};
        sig = rawStruct.(fName);
        
        % Add missing fields
        if ~isfield(sig, 'samplingFreq'), sig.samplingFreq = fs; end
        if ~isfield(sig, 'channelNumber'), sig.channelNumber = 1; end
        if ~isfield(sig, 'classif'), sig.classif = struct(); end
        
        % Default center (will be overwritten by split function later)
        sig.center = 'Prague'; 
        
        signalsOut.(fName) = sig;
    end
end

function signals = applyPragueLabels(signals, splitConfig)
    % Takes ANY signals struct, finds "Prague" entries, and re-assigns 
    % Train/Val/Test based on the seed.
    
    sigNames = fieldnames(signals);
    praguePatientIds = {};
    
    % 1. Collect Patients
    for i = 1:numel(sigNames)
        name = sigNames{i};
        % Check if this signal belongs to Prague (generic or already split)
        c = string(signals.(name).center);
        if contains(c, 'Prague', 'IgnoreCase', true)
            praguePatientIds{end+1} = getPatientId(name); %#ok<AGROW>
        end
    end
    
    uniquePatients = unique(praguePatientIds);
    numPatients = numel(uniquePatients);

    % 2. Calculate Split
    rng(splitConfig.seed);
    shuffledPatients = uniquePatients(randperm(numPatients));

    nTrain = round(splitConfig.ratios.train * numPatients);
    nVal   = round(splitConfig.ratios.val * numPatients);

    trainSet = shuffledPatients(1:nTrain);
    valSet   = shuffledPatients(nTrain+1 : nTrain+nVal);
    testSet  = shuffledPatients(nTrain+nVal+1 : end);

    % 3. Apply Labels
    for i = 1:numel(sigNames)
        name = sigNames{i};
        c = string(signals.(name).center);

        if contains(c, 'Prague', 'IgnoreCase', true)
            pId = getPatientId(name);

            if ismember(pId, trainSet)
                signals.(name).center = 'PragueTrain';
            elseif ismember(pId, valSet)
                signals.(name).center = 'PragueVal';
            elseif ismember(pId, testSet)
                signals.(name).center = 'PragueTest';
            else
                signals.(name).center = 'PragueUnknown';
            end
        end
    end
end

function merged = mergeSignalStructs(s1, s2)
    % Merges two structures. If field names collide, it warns and keeps s1.
    merged = s1;
    f2 = fieldnames(s2);
    
    for i = 1:numel(f2)
        fn = f2{i};
        if isfield(merged, fn)
            warning('Duplicate signal ID found during merge: %s. Keeping original.', fn);
        else
            merged.(fn) = s2.(fn);
        end
    end
end

function statsTable = generateParticipantStatistics(signals)
    % Iterates through signals and aggregates stats by Participant ID
    
    sigNames = fieldnames(signals);
    
    % Map to store aggregated data: Key = ParticipantID
    % Value = Struct(Center, NumSignals, TotalWins, CleanWins, ArtifWins)
    pMap = containers.Map();
    
    fprintf('Processing %d signals...\n', numel(sigNames));
    
    for i = 1:numel(sigNames)
        name = sigNames{i};
        sig  = signals.(name);
        
        % 1. Extract Participant ID based on naming convention
        pID = getPatientId(name);
        
        % 2. Identify Center
        if isfield(sig, 'center') && ~isempty(sig.center)
            centerName = string(sig.center);
        end
        
        % 3. Count Windows (Ground Truth)
        if isfield(sig, 'artif') && ~isempty(sig.artif)
            % Ground truth vector (0 = clean, >0 = artifact)
            annot = sig.artif(:); 
            nTotal = numel(annot);
            nArt   = sum(annot ~= 0);
            nClean = sum(annot == 0);
        else
            warning('Signal %s has no artifact annotations.', name);
            nTotal = 0; nArt = 0; nClean = 0;
        end
        
        % 4. Aggregate Results
        if isKey(pMap, pID)
            curr = pMap(pID);
            curr.NumSignals      = curr.NumSignals + 1;
            curr.TotalWindows    = curr.TotalWindows + nTotal;
            curr.CleanWindows    = curr.CleanWindows + nClean;
            curr.ArtifactWindows = curr.ArtifactWindows + nArt;
            
            % Update center if it was previously unknown/generic
            if curr.Center == "Prague" && centerName ~= "Prague"
                curr.Center = centerName;
            end
            
            pMap(pID) = curr;
        else
            % Initialize new participant entry
            s = struct();
            s.Center          = centerName;
            s.NumSignals      = 1;
            s.TotalWindows    = nTotal;
            s.CleanWindows    = nClean;
            s.ArtifactWindows = nArt;
            pMap(pID) = s;
        end
    end
    
    % 5. Convert Map to Table
    allIDs = keys(pMap);
    numP   = numel(allIDs);
    
    % Pre-allocate arrays for table
    t_ID        = string(allIDs(:));
    t_Center    = strings(numP, 1);
    t_NumSig    = zeros(numP, 1);
    t_TotalWin  = zeros(numP, 1);
    t_CleanWin  = zeros(numP, 1);
    t_ArtWin    = zeros(numP, 1);
    
    for k = 1:numP
        data = pMap(allIDs{k});
        t_Center(k)   = data.Center;
        t_NumSig(k)   = data.NumSignals;
        t_TotalWin(k) = data.TotalWindows;
        t_CleanWin(k) = data.CleanWindows;
        t_ArtWin(k)   = data.ArtifactWindows;
    end
    
    statsTable = table(t_ID, t_Center, t_NumSig, t_TotalWin, t_CleanWin, t_ArtWin, ...
        'VariableNames', {'ParticipantID', 'Center', 'NumSignals', ...
                          'TotalWindows', 'CleanWindows', 'ArtifactWindows'});
                      
    % Sort by Center then ID
    statsTable = sortrows(statsTable, {'Center', 'ParticipantID'});
end


function [loadedSignals, stats] = filterSignalsByCenter(loadedSignals, targetPrefix)
% FILTERSIGNALSBYCENTER Robustly removes signals not matching targetPrefix.
% Handles empty fields, non-string data, and missing fields safely.

    if nargin < 2 || isempty(targetPrefix), targetPrefix = 'Prague'; end

    sigNames = fieldnames(loadedSignals);
    if isempty(sigNames)
        stats = struct('numOriginal', 0, 'numRemoved', 0, 'numKept', 0);
        return; 
    end

    fprintf('--- Filtering by Center: "%s" ---\n', targetPrefix);

    % 1. Use structfun to create a logical mask
    % We define an anonymous function that:
    %  - Checks if 'center' exists
    %  - Checks if it's actual text
    %  - Checks if it starts with our prefix
    keepMask = structfun(@(s) isfield(s, 'center') && ...
                             (ischar(s.center) || isstring(s.center)) && ...
                             startsWith(s.center, targetPrefix, 'IgnoreCase', true), ...
                        loadedSignals);

    % 2. Identify signals to remove
    toRemove = sigNames(~keepMask);
    numRemoved = numel(toRemove);
    
    % 3. Execute removal in one batch
    if numRemoved > 0
        loadedSignals = rmfield(loadedSignals, toRemove);
    end

    % Statistics
    stats.numOriginal = numel(sigNames);
    stats.numRemoved  = numRemoved;
    stats.numKept     = sum(keepMask);
    
    fprintf('   Kept %d signals. Removed %d (missing or non-matching).\n', ...
            stats.numKept, stats.numRemoved);
end