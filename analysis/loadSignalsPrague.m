function [loadedSignals] = loadSignalsPrague(dataFolder, csvFile, signalsFolder, saveData, savePath)
    % Loads signal data and annotations from a folder into structure.
    % 
    % Inputs:
    %   - dataFolder: Path to the folder containing the CSV file and signals folder
    %   - csvFile: Name of the CSV file with metadata (e.g., '_metadataMER2020.csv')
    %   - signalsFolder: Name of the folder containing .mat files with signals
    %   - saveData (optional): Boolean indicating whether to save the data (default: false)
    %   - savePath (optional): Path to save the .mat file (default: 'loadedSignals.mat' in dataFolder)
    %
    % Outputs:
    %   - loadedSignals: Structure containing loaded signals and annotations

    if nargin < 4
        saveData = false; % Default: do not save data
    end
    if nargin < 5
        savePath = fullfile(dataFolder, 'loadedSignals.mat'); % Default save path
    end

    if saveData
        fprintf(['Loading signals with annotations has started, ' ...
            'output will be saved into: %s'],savePath)
    else 
        fprintf(['Loading signals with annotations has started, ' ...
            'output will not be saved'])
    end

    % Load the CSV file with annotations
    annotations = readtable(fullfile(dataFolder, csvFile));

    % Identify sigId values with annotations (non-NaN in any artifact column)
    artifactColumns = contains(annotations.Properties.VariableNames, 'artifacts');
    hasAnnotations = any(~isnan(table2array(annotations(:, artifactColumns))), 2);
    sigIdsWithAnnotations = annotations.sigId(hasAnnotations);

    numAnnotatedSignals = length(sigIdsWithAnnotations);
    fprintf('%d from %d signals have annotations.\n', numAnnotatedSignals, size(annotations, 1));

    % Initialize a structure to store loaded signals
    loadedSignals = struct();

    % Loop through the sigId values and load corresponding .mat files
    progressMarks = round([0.25, 0.50, 0.75, 1] * numAnnotatedSignals);
    for i = 1:numAnnotatedSignals
        sigId = sigIdsWithAnnotations{i};
        matFilePath = fullfile(dataFolder, signalsFolder, [sigId, '.mat']);

        if exist(matFilePath, 'file')
            % Load the .mat file
            signalData = load(matFilePath);

            % Generate a valid field name for the structure
            validFieldName = matlab.lang.makeValidName(['sig_', sigId]);

            % Store the signal and annotations in the structure
            loadedSignals.(validFieldName).area = signalData.area;
            loadedSignals.(validFieldName).artif = signalData.artif; % Annotation data
            loadedSignals.(validFieldName).data = signalData.data;   % Signal data
            loadedSignals.(validFieldName).sigId = signalData.sigId; % Signal ID
        else
            warning('File %s does not exist.', matFilePath);
        end

        % Report progress at 25%, 50%, and 75%
        if ismember(i, progressMarks)
            fprintf('Progress: %.0f%% signals loaded\n', (i / numAnnotatedSignals) * 100);
        end
    end

    % Optionally save the data
    if saveData
        disp('Saving the data..')
        save(savePath, 'loadedSignals', '-v7.3');
        fprintf('Data saved to %s\n', savePath);
    end

    disp('Signals with annotations successfully loaded.');
end
