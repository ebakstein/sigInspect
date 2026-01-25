function loadedSignals = convertLoadedSignalsAnnotations(loadedSignals)
% CONVERTLOADEDSIGNALSANNOTATIONS
% Converts numeric annotation fields (.artif) in loadedSignals struct
% into binary and multi-label matrices.
%
% Inputs:
%   loadedSignals - struct of signals loaded via loadSignalsFromMat()
%
% Output:
%   loadedSignals - struct with .artifBin and .artifMulti field (converted binary annotations)
%


    if ~isstruct(loadedSignals)
        error('Input must be a struct (loadedSignals).');
    end

    sigNames = fieldnames(loadedSignals);
    fprintf('Converting annotations for %d signals...\n', numel(sigNames));

    % Infer number of artifact types from all .artif vectors
    allVals = [];
    for i = 1:numel(sigNames)
        a = loadedSignals.(sigNames{i}).artif;
        if ~isempty(a)
            allVals = [allVals, a(:)']; %#ok<AGROW>
        end
    end
    if isempty(allVals)
        warning('No annotations found. Nothing to convert.');
        return;
    end
    maxN = floor(log2(max(allVals))) + 1;

    % Conversion loop
    for i = 1:numel(sigNames)
        sig = loadedSignals.(sigNames{i});
        if ~isfield(sig, 'artif') || isempty(sig.artif)
            continue;
        end

        artif = double(string(sig.artif));

        % Store new binary annotation
        loadedSignals.(sigNames{i}).artifBin = double(artif ~= 0);
        loadedSignals.(sigNames{i}).artifMulti = annotNum2Bin(double(string(artif)), maxN)'; % (Ntypes × Nwindows)
    end

    fprintf('Annotation conversion complete. Added field: .artifBin\n');
end
