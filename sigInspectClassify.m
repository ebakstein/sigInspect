function annot = sigInspectClassify(signal, signalId, fs, method, varargin)
% function annot = sigInspectClassify(signal,varargin)
% classify artifacts in each second of provided micro-EEG signal
% 
% IN: 
%   signal - micro-EEG signal vector (or signal matrix with channels in rows)
%   fs     - sampling frequency in Hz
%   method - classification method:
%           'psd' - normalized PSD spectrum thresholding, based on [1]
%                  (default) 91% train /88% test set accuracy on EMBC data
%                   threshold based on [2]
%           'tree'- pre-trained decision tree, based on multiple features,
%                   trained on a multi-centric database from [2]
%           'cov' - not yet implemented
%           'svm' - 3 svm classifiers - POW, BASE, FREQ artifacts
%   params - optional parameters for some of the classification methods:
%          for 'psd': detection trheshold (default:0.01)
%          for 'cov' three parameters:
%                    - threshold (default 1.2 based on [2])
%                    - winLength: length of signal segment (default 0.25s)
%                    - aggregPerc: what proportion of one second window has
%                      to be marked as artifact to mark the whole sec. as
%                      artifact (deafault: winLength)
%           (default values based on [2])
%           for 'svm': struct containing artifact names and optionally threshold values 
% OUT:
%   annot - logical vector of annotation for each second of input signal.
%           true = artifact, false = clean signal
% 
% E. Bakstein 2015-06-29
% 
% [1] Bakstein, E. et. al.: Supervised Segmentation of Microelectrode Recording Artifacts Using Power Spectral Density, in Proceedings of IEEE EMBS, 2015
% [2] Bakstein, E.: Deep Brain Recordings in Parkinson's Disease: Processing, Analysis and Fusion with Anatomical Models, Doctoral Thesis, October 2016


% input checks
if(nargin<2)
    error('sampling frequency must be specified')
end

if(nargin<3 || isempty(method))
    method='psd';
end
    
if(isempty(signal))
    annot=[];
    return;
end

if(length(signal)<fs)
    error('signal must be at least 1s long');
end

[Nch,N] = size(signal); % number of channels + samples
Ns = ceil(N/fs);       % number of seconds


% ---- CLASSIFIER PARAMETERS ----

% features to be computed
switch(method)
    
    case {'psd','psdPrg'}
        % Default threshold values
        if strcmpi(method,'psd')
            psdThr = 0.01;   % trained on multi-center data
        else
            psdThr = 0.0085; % trained on Prague data
        end

        % If user passed a custom numeric threshold, override
        if nargin > 3 && ~isempty(varargin) && isnumeric(varargin{1}) && isscalar(varargin{1})
            psdThr = varargin{1};
        elseif nargin > 3 && ~isempty(varargin) && ~(isnumeric(varargin{1}) && isempty(varargin{1}))
            warning('Invalid PSD threshold argument. Using default value: %.3f', psdThr);
        end

        % Safety check (avoid non-scalar logic error)
        if isnumeric(psdThr) && isscalar(psdThr)
            if psdThr < 0 || psdThr > 0.03
                warning('Recommended PSD threshold range is 0.005–0.02. Value %.3f may lead to unexpected results.', psdThr);
            end
        else
            warning('PSD threshold was not a valid scalar. Using default value: %.3f', 0.01);
            psdThr = 0.01;
        end

        % features
        featNames={'maxNormPSD'}; % definition of dataset columns
        featComp = [1];           % features actually computet (for compatibility with dec. tree)        
        method = 'psd';
    
    case {'tree', 'treePrg'}
        % load classifier params
        try
            classif = load('sigInspectClassifiers.mat'); % load pre-trained classifiers
        catch err        
            error('could not load precalculated classifiers: search sigInspect root for the file sigInspectClassifers.mat (necessary for the tree-based classifiers)')
        end
        if(strcmp(method,'tree'))
            classif.tree = classif.treeAll;
        else
            classif.tree = classif.treePrg;
        end
        % features to compute
        featNames = classif.featNames; % all 19 features have to be in the set 
        featComp = setdiff(unique(classif.tree.var),0); % only some are needed (non-nan)
        method = 'tree';
    case 'cov'
        featNames = {};
        covThr = 1.2;
        winLength = .25;
        aggregPerc = winLength;

        % Unwrap any cells that contain numeric values
        if numel(varargin) == 1 && iscell(varargin{1}) && numel(varargin{1}) > 1
            args = varargin{1};
        else
            args = varargin;
        end
        
        % Unwrap each cell value safely
        args = cellfun(@(x) unwrapCell(x), args, 'UniformOutput', false);

        % Threshold (arg 1)
        if numel(args) >= 1 && isnumeric(args{1}) && isscalar(args{1}) && args{1} >= 1
            covThr = args{1};
        elseif numel(args) >= 1 && ~isempty(args{1})
            warning('Invalid COV threshold argument. Using default value: %.2f', covThr);
        end

        % Window length (arg 2)
        if numel(args) >= 2 && isnumeric(args{2}) && isscalar(args{2}) && args{2} > 0 && args{2} < 1
            winLength = args{2};
        elseif numel(args) >= 2 && ~isempty(args{2})
            warning('Invalid COV window length argument. Using default value: %.2f', winLength);
        end

        % Aggregation percentage (arg 3) 
        if numel(args) >= 3 && isnumeric(args{3}) && isscalar(args{3}) && args{3} > 0 && args{3} <= 1
            aggregPerc = args{3};
        elseif numel(args) >= 3 && ~isempty(args{3})
            warning('Invalid COV aggregation argument. Using default value: %.2f', aggregPerc);
        end

    case 'svm'
        try
            classif = load('sigInspectSVMClassifiers.mat'); % load pre-trained classifiers
            classif = classif.classifiers;
        catch err        
            error('could not load precalculated classifiers: search sigInspect root for the file sigInspectSVMClassifiers.mat (necessary for the svm classifiers)')
        end
        if nargin >= 5 && ~isempty(varargin)
            param = varargin{1};
        else
            param = struct( ...
                'POW', struct(), ...
                'BASE', struct(), ...
                'FREQ', struct()); 
        end

        % Extract artifact types from user's input struct (fields of param)
        artifactTypes = intersect(fieldnames(param), fieldnames(classif));  % make sure it's valid
    
        if isempty(artifactTypes)
            artifactTypes = fieldnames(classif);  % default to all if none provided
        end
    
        % Fill in default threshold where not provided
        thresholds = struct();
        for i = 1:numel(artifactTypes)
            art = artifactTypes{i};
            if isfield(param.(art), 'threshold')
                thresholds.(art) = param.(art).threshold;
            else
                thresholds.(art) = 0.5;
            end
        end
    
        % Collect unique relevant feature names
        featNames = {};
        for i = 1:numel(artifactTypes)
            featList = classif.(artifactTypes{i}).featNames;
            featNames = [featNames, featList];  % concatenate all
        end
        featNames = unique(featNames);  % remove duplicates
        featComp = 1:length(featNames);

        % Init annotation matrix
        Nartif = length(artifactTypes);
        method = 'svm';

    otherwise
        error('Unknown method: %s',method)
end

% ---- COMPUTE FEATURES ----
cacheFile = 'featuresCache.mat'; % or full path if needed
if ~strcmpi(method, 'cov')
    featVals = getOrComputeFeatures(signal, signalId, method, featNames, fs, cacheFile);
else
    featVals = [];
end


% ---- CLASSIFY ----
switch(method)
    case 'psd'
        % compare to preset threshold
        annot = featVals>psdThr;
        % change dims to be Nch*Ns
        annot = reshape(annot,Nch,Ns);
    case 'tree'
        % classify using decision tree
        annot = eval(classif.tree,featVals); 
        annot=strcmp(annot,'1');
        % change dims to be Nch*Ns
        annot = reshape(annot,Nch,Ns);
    case 'cov'
        annot = false(Nch,Ns);
        for chi=1:Nch
            annot(chi,:) = sigInspectClassifyCov(signal(chi,:),fs,'cov', covThr, winLength, aggregPerc,false);
        end
    case 'svm'
        artifactTypes = fieldnames(param);  % Update artifact types from param keys
        Nartif = length(artifactTypes);
        annot = false(Nch, Ns, Nartif);
        scores = zeros(Nch*Ns, Nartif);  % soft outputs
        for i = 1:Nartif
            art = artifactTypes{i};
            featureSet = classif.(art).featNames;  % features for this artifact
            [~, featureInds] = ismember(featureSet, featNames);  % indices in featVals, in correct order
            if any(featureInds == 0)
                error('Some features required by the SVM are missing in the computed features!');
            end
            X = featVals(:, featureInds);

            % Get classifier
            model = classif.(art).svmProbModel;

            % Check for NaNs and replace with zeros
            X(isnan(X)) = 0;

            % Predict using trained SVM – get scores (posterior probabilities)
            [~, score] = predict(model, X);  % score is N x 2, we take the 2nd column (class == 1)
            scores(:, i) = score(:, 2);  % probability of being artifact

            % Apply threshold to convert to binary
            annot(:,:,i) = reshape(scores(:, i) > thresholds.(art), Nch, Ns);
        end
end

function featVals = getOrComputeFeatures(signal, signalId, method, featNames, fs, cacheFile)
    % signal: matrix (channels x samples)
    % signalId: string, e.g. '#1'
    % method: string, e.g. 'svm'
    % featNames: cell array of feature names
    % fs: sampling frequency
    % cacheFile: string, path to .mat file

    % Make valid field names
    signalIdField = matlab.lang.makeValidName(signalId);
    methodField = matlab.lang.makeValidName(method);

    % Try to load cache
    featuresCache = struct();
    if exist(cacheFile, 'file')
        S = load(cacheFile);
        % If the file contains a featuresCache struct, unwrap it to top-level fields
        if isfield(S, 'featuresCache')
            % Unwrap: assign each field of featuresCache to the workspace struct
            fcFields = fieldnames(S.featuresCache);
            for k = 1:numel(fcFields)
                featuresCache.(fcFields{k}) = S.featuresCache.(fcFields{k});
            end
        else
            % Already top-level signalId fields
            featuresCache = S;
        end
    end

    % Compute expected dimensions
    Nch = size(signal, 1);
    N = size(signal, 2);
    Ns = ceil(N / fs);
    Nfeat = numel(featNames);
    expectedSize = [Nch * Ns, Nfeat];

    % Check if features already cached
    if isfield(featuresCache, signalIdField) && isfield(featuresCache.(signalIdField), methodField)
        featVals = featuresCache.(signalIdField).(methodField);

        % Verify cached size
        if isequal(size(featVals), expectedSize)
            disp(['Using cached features for ', signalIdField, ' (', method, ')']);
            return;
        else
            warning(['Cached features for ', signalIdField, ' (', method, ...
                ') have mismatched size. Expected [', num2str(expectedSize), ...
                '], got [', num2str(size(featVals)), ']. Recomputing...']);
        end
    end

    % Compute features
    featVals = nan(expectedSize);
    for si = 1:Ns
        inds = (si-1)*fs+1 : min(si*fs, N);
        fv = sigInspectComputeFeatures(signal(:,inds), featNames, fs);
        nChHere = size(fv,1);
        featVals((si-1)*Nch + (1:nChHere), :) = fv;
    end

    % Save to cache (append or create) as top-level signalId fields
    featuresCache.(signalIdField).(methodField) = featVals;
    if exist(cacheFile, 'file')
        save(cacheFile, '-struct', 'featuresCache', '-append');
    else
        save(cacheFile, '-struct', 'featuresCache');
    end
end

function val = unwrapCell(x)
    % unwrap { [value] } → value
    if iscell(x)
        val = x{1};
    else
        val = x;
    end
end

end

