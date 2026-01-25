function plotMetricsByCenter(loadedSignals, targetMethod)
% plotMetricsByCenter Plots ROC (top row) and PR (bottom row) curves.
% 
% Inputs:
%   loadedSignals - Struct containing classification results and ground truth.
%   targetMethod  - (Optional) String, e.g., 'svm', 'psd'. If omitted, plots all.

    sigNames = fieldnames(loadedSignals);
    if isempty(sigNames), return; end
    
    % Identify available methods
    firstSig = loadedSignals.(sigNames{1});
    if ~isfield(firstSig, 'classif')
        error('No classification results found in signals.');
    end
    allAvailableMethods = fieldnames(firstSig.classif);
    
    % Determine which methods to process
    if nargin < 2 || isempty(targetMethod)
        methodsToProcess = allAvailableMethods;
    else
        % Ensure the target method exists
        if ~ismember(targetMethod, allAvailableMethods)
            error('Method "%s" not found in loadedSignals.', targetMethod);
        end
        methodsToProcess = {targetMethod};
    end

    centers = unique(structfun(@(x) string(x.center), loadedSignals));
    colors = lines(numel(centers));

    for m = 1:numel(methodsToProcess)
        method = methodsToProcess{m};
        
        % Skip if no scores exist
        if ~isfield(firstSig.classif.(method), 'scores')
            fprintf('Skipping %s: No probability scores found.\n', method);
            continue; 
        end
        
        % Setup label logic
        currentScores = firstSig.classif.(method).scores;
        nArt = size(currentScores, 3);
        isSVM = strcmpi(method, 'svm');
        
        if isSVM
            gtField = 'artifMulti'; 
            artLabels = {'POW', 'BASE', 'FREQ'}; 
        else
            gtField = 'artifBin';
            artLabels = {'General Artifact'};
            nArt = 1; 
        end

        % Create Figure
        figName = sprintf('Performance: %s', method);
        figure('Name', figName, 'Color', 'w', 'Position', [50 100 400 * nArt, 700]);

        for a = 1:nArt
            % --- Subplot Selection ---
            % Top row: Indices 1 to nArt
            axRoc = subplot(2, nArt, a); hold on; 
            % Bottom row: Indices nArt+1 to 2*nArt
            axPr  = subplot(2, nArt, a + nArt); hold on; 
            
            title(axRoc, sprintf('%s ROC: %s', method, artLabels{a}));
            title(axPr,  sprintf('%s PR: %s',  method, artLabels{a}));

            for c = 1:numel(centers)
                centerName = centers(c);
                allScores = [];
                allTrue = [];
                
                for s = 1:numel(sigNames)
                    sig = loadedSignals.(sigNames{s});
                    if string(sig.center) ~= centerName || ~isfield(sig, gtField), continue; end
                    
                    % 1. Filter out "Unknown" category (Row 5)
                    if isSVM && size(sig.(gtField), 1) >= 5
                        if any(sig.(gtField)(5, :) == 1), continue; end
                    end
                    
                    % 2. Extract and Flatten
                    yPredProb = squeeze(sig.classif.(method).scores(:,:,a));
                    yPredProb = yPredProb(:);
                    
                    if isSVM
                        yTrue = sig.(gtField)(a, :); 
                    else
                        yTrue = sig.(gtField);
                    end
                    yTrue = yTrue(:);
                    
                    % 3. Align Dimensions
                    nLen = min(numel(yPredProb), numel(yTrue));
                    if nLen == 0, continue; end
                    s_vec = yPredProb(1:nLen);
                    t_vec = double(yTrue(1:nLen));
                    
                    % 4. Filter NaNs
                    validIdx = ~isnan(s_vec) & ~isnan(t_vec);
                    allScores = [allScores; s_vec(validIdx)];
                    allTrue   = [allTrue;   t_vec(validIdx)];
                end
                
                if isempty(allScores) || numel(unique(allTrue)) < 2, continue; end

                % --- Compute and Plot ROC ---
                [Xroc, Yroc, ~, aucRoc] = perfcurve(allTrue, allScores, 1);
                plot(axRoc, Xroc, Yroc, 'Color', colors(c,:), 'LineWidth', 2, ...
                    'DisplayName', sprintf('%s (AUC: %.2f)', centerName, aucRoc));
                
                % --- Compute and Plot PR ---
                [Xpr, Ypr, ~, aucPr] = perfcurve(allTrue, allScores, 1, 'XCrit', 'recall', 'YCrit', 'prec');
                plot(axPr, Xpr, Ypr, 'Color', colors(c,:), 'LineWidth', 2, ...
                    'DisplayName', sprintf('%s (AUC: %.2f)', centerName, aucPr));
            end
            
            % --- Polish ROC Axis ---
            plot(axRoc, [0 1], [0 1], 'k--', 'HandleVisibility', 'off');
            grid(axRoc, 'on'); box(axRoc, 'on');
            xlabel(axRoc, 'FPR'); ylabel(axRoc, 'TPR');
            legend(axRoc, 'Location', 'southeast', 'FontSize', 7);
            axis(axRoc, 'square');

            % --- Polish PR Axis ---
            grid(axPr, 'on'); box(axPr, 'on');
            xlabel(axPr, 'Recall'); ylabel(axPr, 'Precision');
            legend(axPr,  'Location', 'northeast', 'FontSize', 7);
            axis(axPr, 'square');
        end
    end
end