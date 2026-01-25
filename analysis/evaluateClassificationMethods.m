function [resultsTable, loadedSignals] = evaluateClassificationMethods(loadedSignals)
% EVALUATECLASSIFICATIONMETHODS Computes performance metrics for all signals.
% Robust version: Skips signals with empty classification results and warns user.

    % Configuration
    binaryMethods   = {'psd', 'cov', 'svmunion'}; 
    svmMethod       = 'svm';
    svmLabels       = {'POW','BASE','FREQ'};
    overallCenter   = 'All';
    
    signalNames = fieldnames(loadedSignals);
    allCounts = struct();
    
    % --- Helper to Accumulate Counts ---
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
    
    fprintf('Evaluating %d signals...\n', numel(signalNames));
    
    for i = 1:numel(signalNames)
        name = signalNames{i};
        sig  = loadedSignals.(name);
        
        % 1. Detect Center
        if isfield(sig, 'center') && ~isempty(sig.center)
            center = sig.center;
        elseif contains(name, 'Bra'), center = 'Bra';
        elseif contains(name, 'Brn'), center = 'Brn';
        elseif contains(name, 'Kos'), center = 'Kos';
        elseif contains(name, 'Olo'), center = 'Olo';
        else, center = 'Unknown';
        end
        
        if ~isfield(sig,'classif'), continue; end
        
        % 2. Get Ground Truth
        if isfield(sig, 'artifBin'), yTrueBin = logical(sig.artifBin); else, continue; end
        if isfield(sig, 'artifMulti'), yTrueMulti = logical(sig.artifMulti); else, yTrueMulti = []; end
        
        % 3. Evaluate Binary Methods (PSD, COV, SVMUNION)
        for m = 1:numel(binaryMethods)
            method = binaryMethods{m};
            if isfield(sig.classif, method) && isfield(sig.classif.(method), 'annot')
                pred = logical(sig.classif.(method).annot);
                
                % Skip if pred is empty
                if isempty(pred)
                    fprintf('Warning: Empty binary result for %s (Method: %s)\n', name, method);
                    continue; 
                end

                nLen = min(numel(yTrueBin), numel(pred));
                if nLen == 0, continue; end % Skip if no overlap
                
                [TP,FP,TN,FN] = getCounts(yTrueBin(1:nLen), pred(1:nLen));
                acc(method, center, 'bin', TP,FP,TN,FN);
                acc(method, overallCenter, 'bin', TP,FP,TN,FN);
            end
        end
        
        % 4. Evaluate Multilabel SVM (POW, BASE, FREQ)
        if isfield(sig.classif, svmMethod) && ~isempty(yTrueMulti)
            % Skip signals with 'Unknown' artifact type (row 5)
            if size(yTrueMulti,1) >= 5 && any(yTrueMulti(5,:) == 1), continue; end
            
            pred = sig.classif.(svmMethod).annot;
            
            % Check for empty prediction
            if isempty(pred)
                fprintf('Warning: SVM prediction is empty for signal: "%s". Skipping.\n', name);
                continue; 
            end
            
            for l = 1:numel(svmLabels)
                if size(pred,3) >= l
                    % This line caused the crash because pred was empty
                    try
                        yPred = squeeze(pred(1,:,l)); 
                        
                        if ~islogical(yPred), yPred = yPred > 0.5; end
                        
                        yTrue = yTrueMulti(l,:);
                        nLen = min(numel(yTrue), numel(yPred));
                        if nLen == 0, continue; end
                        
                        [TP,FP,TN,FN] = getCounts(yTrue(1:nLen), yPred(1:nLen));
                        acc(svmMethod, center, svmLabels{l}, TP,FP,TN,FN);
                        acc(svmMethod, overallCenter, svmLabels{l}, TP,FP,TN,FN);
                    catch ME
                        fprintf('Error evaluating %s (Label %s): %s\n', name, svmLabels{l}, ME.message);
                    end
                end
            end
        end
    end
    
    % 5. Build Final Table
    rows = [];
    keys = fieldnames(allCounts);
    for k = 1:numel(keys)
        key = keys{k};
        parts = split(string(key), "_");
        if numel(parts)<3, continue; end
        
        method = parts{1}; center = parts{2}; label = parts{3};
        counts = allCounts.(key);
        m = countsToMetrics(counts);
        
        newRow = table(string(method), string(center), string(label), ...
            counts.TP, counts.FP, counts.TN, counts.FN, ...
            m.accuracy, m.sensitivity, m.specificity, m.f1, ...
            m.precision, m.recall, m.balancedAcc, ...
            'VariableNames', {'Method','Center','Label','TP','FP','TN','FN', ...
                              'Accuracy','Sensitivity','Specificity','F1','Precision', ...
                              'Recall','BalancedAcc'});
        rows = [rows; newRow]; %#ok<AGROW>
    end
    
    % Sort for readability
    if ~isempty(rows)
        resultsTable = sortrows(rows, {'Method', 'Center', 'Label'});
    else
        resultsTable = table();
    end
end

% --- Internal Helpers ---
function [TP, FP, TN, FN] = getCounts(yTrue, yPred)
    yTrue = logical(yTrue(:)); yPred = logical(yPred(:));
    TP = sum(yTrue & yPred); TN = sum(~yTrue & ~yPred);
    FP = sum(~yTrue & yPred); FN = sum(yTrue & ~yPred);
end

function metrics = countsToMetrics(c)
    TP=c.TP; FP=c.FP; TN=c.TN; FN=c.FN;
    denom = TP+TN+FP+FN; if denom==0, denom=1; end; metrics.accuracy = (TP+TN)/denom;
    denom = TP+FN; if denom==0, denom=1; end; metrics.sensitivity = TP/denom;
    denom = TN+FP; if denom==0, denom=1; end; metrics.specificity = TN/denom;
    denom = TP+FP; if denom==0, denom=1; end; metrics.precision = TP/denom;
    denom = 2*TP+FP+FN; if denom==0, denom=1; end; metrics.f1 = 2*TP/denom;
    metrics.recall = metrics.sensitivity;
    metrics.balancedAcc = (metrics.sensitivity + metrics.specificity)/2;
end