function [loadedSignals, stats] = cleanSignals(loadedSignals, mode, winSec)
% CLEANSIGNALS Filters signals and synchronizes classification results.
    if nargin < 2 || isempty(mode), mode = 'signal'; end
    if nargin < 3 || isempty(winSec), winSec = 1; end 
    
    sigNames = fieldnames(loadedSignals);
    numTotal = numel(sigNames);
    stats = struct('numOriginalSignals', numTotal, 'numRemovedSignals', 0, ...
                   'samplesOriginal', 0, 'samplesRemoved', 0);
               
    fprintf('--- Cleaning Mode: %s (Window: %.2fs) ---\n', mode, winSec);
    keepSignalMask = true(numTotal, 1);

    for i = 1:numTotal
        name = sigNames{i};
        sig  = loadedSignals.(name);
        if ~isfield(sig, 'data') || isempty(sig.data), continue; end
        
        % 1. Force everything to Double for math safety
        data = double(sig.data(:)); 
        N = length(data);
        stats.samplesOriginal = stats.samplesOriginal + N;
        fs = double(sig.samplingFreq);
        hasMulti = isfield(sig, 'artifMulti') && ~isempty(sig.artifMulti);
        
        if strcmpi(mode, 'signal')
            % Logic for whole signal removal
            containsNaN = any(isnan(data));
            containsBadLabel = hasMulti && any(any(sig.artifMulti(6:end, :)));
            if containsNaN || containsBadLabel
                keepSignalMask(i) = false;
                stats.samplesRemoved = stats.samplesRemoved + N;
                continue; 
            end
            if hasMulti, loadedSignals.(name).artifMulti = sig.artifMulti(1:5, :); end

        elseif strcmpi(mode, 'window')
            sampsPerWin = round(winSec * fs);
            if sampsPerWin < 1, continue; end
            
            % 2. Reshape Logic (Fixed the Integer Mismatch)
            remainder = mod(double(N), double(sampsPerWin));
            % Use double() to allow logical multiplication
            padLen = double(remainder > 0) * (double(sampsPerWin) - double(remainder));
            
            dataPadded = [data; zeros(padLen, 1)];
            dataMat = reshape(dataPadded, sampsPerWin, []);
            nWins = size(dataMat, 2);
            
            % Identify bad windows
            badWinMask = any(isnan(dataMat), 1);
            if hasMulti
                % Handle multi-row annotation padding
                annotPadded = [double(sig.artifMulti), zeros(size(sig.artifMulti,1), padLen)];
                annotMat = reshape(annotPadded, size(sig.artifMulti,1), sampsPerWin, nWins);
                badLabelMask = squeeze(any(any(annotMat(6:end, :, :), 1), 2))';
                badWinMask = badWinMask | badLabelMask;
            end
            
            if any(badWinMask)
                % A. Clean Data
                dataCleanMat = dataMat(:, ~badWinMask);
                dataClean = dataCleanMat(:);
                wasLastKept = ~badWinMask(end);
                if wasLastKept && padLen > 0, dataClean = dataClean(1:end-padLen); end
                loadedSignals.(name).data = dataClean'; 
                
                % B. Clean Annotation Fields
                annotFields = {'artif', 'artifBin', 'artifMulti'};
                for f = 1:numel(annotFields)
                    fld = annotFields{f};
                    if isfield(sig, fld) && ~isempty(sig.(fld))
                        currAnnot = double(sig.(fld));
                        numRows = size(currAnnot, 1);
                        
                        targetTotalSamples = nWins * sampsPerWin;
                        currPadded = [currAnnot, zeros(numRows, targetTotalSamples - size(currAnnot, 2))];
                        currMat = reshape(currPadded, size(currAnnot,1), sampsPerWin, nWins);
                        cleanMat = currMat(:, :, ~badWinMask);
                        cleanVec = reshape(cleanMat, size(currAnnot,1), []);
                        if strcmpi(fld, 'artifMulti'), cleanVec = cleanVec(1:5, :); end
                        if wasLastKept && padLen > 0, cleanVec = cleanVec(:, 1:end-padLen); end
                        loadedSignals.(name).(fld) = cleanVec;
                    end
                end

                % C. Synchronize CLASSIF (The Critical Part)
                if isfield(sig, 'classif') && isstruct(sig.classif)
                    methods = fieldnames(sig.classif);
                    for m = 1:numel(methods)
                        mName = methods{m};
                        res = sig.classif.(mName);
                        
                        % Precision check: Do the classification windows match our mask?
                        if isfield(res, 'annot') && ~isempty(res.annot)
                            if size(res.annot, 2) == nWins
                                loadedSignals.(name).classif.(mName).annot = res.annot(:, ~badWinMask, :);
                            else
                                warning('Window mismatch in %s for method %s. Skipping sync.', name, mName);
                            end
                        end
                        if isfield(res, 'scores') && ~isempty(res.scores)
                            if size(res.scores, 2) == nWins
                                loadedSignals.(name).classif.(mName).scores = res.scores(:, ~badWinMask, :);
                            end
                        end
                    end
                end
                
                % D. Metadata
                loadedSignals.(name).Nsamples = length(dataClean);
                removedCount = (sum(badWinMask(1:end-1)) * sampsPerWin) + ...
                               (badWinMask(end) * (sampsPerWin - padLen));
                stats.samplesRemoved = stats.samplesRemoved + removedCount;
            end
        end
    end

    if strcmpi(mode, 'signal')
        loadedSignals = rmfield(loadedSignals, sigNames(~keepSignalMask));
        stats.numRemovedSignals = sum(~keepSignalMask);
    end
    stats = calculateFinalStats(loadedSignals, stats);
end

function stats = calculateFinalStats(loadedSignals, stats)
    stats.numRemainingSignals = numel(fieldnames(loadedSignals));
    stats.samplesRemaining    = stats.samplesOriginal - stats.samplesRemoved;
    stats.percSignalsRemoved  = (stats.numRemovedSignals / max(1,stats.numOriginalSignals)) * 100;
    stats.percSamplesRemoved  = (stats.samplesRemoved / max(1,stats.samplesOriginal)) * 100;
    fprintf('--- Scrubbing Summary: %d signals remaining from %d original count ---\n', stats.numRemainingSignals, stats.numOriginalSignals);
end