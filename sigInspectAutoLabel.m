function [annotation, annotationFile] = sigInspectAutoLabel(interfSignalOrPath, pathToSave, samplingFreq, method, varargin)
% annot = sigInspectAutoLabel(iterfSignalOrPath, pathToSave, samplingFreq, method, params)  
%   label all signals from provided cell array or interface using pre-learned 
%   classifier, return/save annot
%
% IN
%   iterfSignalOrPath - signals to be classified in one of the following
%                       formats:
%           A) cell array - with (single- or multi-channel) signals in each
%           cell. Signals are row vectors, or m x N matrix, where m is number
%           of parallel channels and N number of samples
%           
%           B) interface - signal access interface of type sigInspectDataInterface, 
%           see help sigInspectDataInterface for details
%       
%           C) path to a mat-file with signals formatted according to A),
%           will use sigInspectDataBasic interface to load data - see 
%           help sigInspectDataBasic for details
% 
%           D) empty - will pop-up a gui dialog to load mat-file, formatted
%           according to C)
%   pathToSave - file path to save *.mat file with generated annotation
%                to.(only if no output parameters are set
%                default: sigInspectAutoAnnotationyyyy-mm-dd-HHMMSS.mat.
%   samplingFreq - sampling frequency in Hz
%                may be left blank if set in interface
%   method     - classifier used to label data - see sigInspectClassify for
%                details
%   params     - one or more parameters for the selected method. See
%                sigInspectClassify for details.
% 
% OUT
%   annot - annotation structure (optional)
% 
% E. Bakstein 2015-06-26
% UPDATES: 20161010 - reworked, additional classifiers (tree, cov)
% 

fprintf('----------- sigInspectAutoLabel -------------\n')
if(nargin<4)
    method=[]; % use default of sigInspectClassify
    fprintf('Classification method unspecified, using psd(default)\n')
end
fprintf('initialization ')

interface = initInterface(interfSignalOrPath);
signalIds = interface.getSignalIds;
N = length(signalIds);

if(N==0)
    error('No signals loaded')
end

% sampling freq. - from interface or so
if(isprop(interface,'settings') && isfield(interface.settings,'SAMPLING_FREQ'))
    samplingFreq = interface.settings.SAMPLING_FREQ;
    fprintf('samplingFrequency = %d Hz (from interface)\n',samplingFreq);
else
    if(nargin<3 || isempty(samplingFreq))
        samplingFreq = 24000;
        fprintf('samplingFrequency = %d Hz (DEFAULT)\n',samplingFreq);
    else
        fprintf('samplingFrequency = %d Hz (from PARAMETER)\n',samplingFreq);        
    end
end

defaultSvmArtifacts = {'POW', 'BASE', 'FREQ'};

% Get artifact types from interface or use defaults
originalArtifactTypes = {};  % Store original types for preservation
if(isprop(interface,'settings') && isfield(interface.settings,'ARTIFACT_TYPES'))
    interfaceArtifactTypes = interface.settings.ARTIFACT_TYPES;
    originalArtifactTypes = interfaceArtifactTypes;  % Keep original order
    fprintf('interface_artifact_types = %s (from interface)\n', strjoin(interfaceArtifactTypes, ', '));
    
    if strcmpi(method, 'svm')
        % For SVM, find intersection between interface types and default SVM types
        [artifactTypesToProcess, ~, ~] = intersect(interfaceArtifactTypes, defaultSvmArtifacts, 'stable');
        if isempty(artifactTypesToProcess)
            % If no intersection, use default SVM types and warn user
            artifactTypesToProcess = defaultSvmArtifacts;
            fprintf('WARNING: No common artifact types found between interface (%s) and SVM defaults (%s)\n', ...
                strjoin(interfaceArtifactTypes, ', '), strjoin(defaultSvmArtifacts, ', '));
            fprintf('Using SVM defaults: %s\n', strjoin(artifactTypesToProcess, ', '));
            % Use SVM defaults as the final artifact types (no existing annotations to preserve)
            artifactTypes = artifactTypesToProcess;
        else
            fprintf('Using intersected artifact types: %s\n', strjoin(artifactTypesToProcess, ', '));
            % Use original interface types as final types (to preserve order and existing annotations)
            artifactTypes = interfaceArtifactTypes;
        end
    else
        % For non-SVM methods, use interface types as-is
        artifactTypes = interfaceArtifactTypes;
        artifactTypesToProcess = interfaceArtifactTypes;
    end
    
    fprintf('number_of_artifact_types = %d (processed from interface)\n', length(artifactTypesToProcess));
else
    % No interface settings - use method-specific defaults
    if strcmpi(method, 'svm')
        artifactTypes = defaultSvmArtifacts;
        artifactTypesToProcess = defaultSvmArtifacts;
        fprintf('number_of_artifact_types = %d (DEFAULT SVM: %s)\n', length(artifactTypes), strjoin(artifactTypes, ', '));
    else
        artifactTypes = {'ARTIF','UNSURE'};
        artifactTypesToProcess = {'ARTIF','UNSURE'};
        fprintf('number_of_artifact_types = %d (DEFAULT NON-SVM: %s)\n', length(artifactTypes), strjoin(artifactTypes, ', '));
    end
end
Nartif = length(artifactTypes);

% default field to store automatic artifact in (for non-SVM methods)
if(isprop(interface,'settings') && isfield(interface.settings,'ARTIFACT_AUTOLABEL_WHICH'))
    artifactAutoWhich = interface.settings.ARTIFACT_AUTOLABEL_WHICH;
    if artifactAutoWhich > Nartif
        fprintf('WARNING: ARTIFACT_AUTOLABEL_WHICH (%d) exceeds number of artifact types (%d), using 1\n', artifactAutoWhich, Nartif);
        artifactAutoWhich = 1;
    end
    fprintf('auto artifact will be stored at position %d (%s)(from interface)\n',artifactAutoWhich,artifactTypes{artifactAutoWhich});
elseif(Nartif > 1 && ~strcmpi(method, 'svm'))
    artifactAutoWhich = 1;
    fprintf('auto artifact will be stored at position %d (%s)(DEFAULT)\n',artifactAutoWhich,artifactTypes{artifactAutoWhich});
else
    artifactAutoWhich = 1;
    if ~strcmpi(method, 'svm')
        fprintf('auto artifact will be stored at position %d (%s)(DEFAULT)\n',artifactAutoWhich,artifactTypes{artifactAutoWhich});
    end
end

% Check if interface has existing annotations to preserve
hasExistingAnnotations = false;
existingAnnotations = {};
if isprop(interface, 'getAnnotationsAll')
    try
        existingAnnotations = interface.getAnnotationsAll();
        if ~isempty(existingAnnotations)
            hasExistingAnnotations = true;
            fprintf('Found existing annotations in interface - will preserve non-updated types\n');
        end
    catch
        fprintf('No existing annotations found or interface does not support getAnnotationsAll()\n');
    end
end

% init empty annotation array
annotation=cell(length(signalIds),1);

% SVM threshold GUI logic
if strcmpi(method, 'svm')
    showThresholdGUI = false;
    if isempty(varargin) || isempty(varargin{1}) || ~isstruct(varargin{1})
        showThresholdGUI = true;
        initialVals = [0.5 0.5 0.5];
        param = struct();
    else
        param = varargin{1};
        % Check which of our artifact types have thresholds defined
        providedTypes = intersect(fieldnames(param), artifactTypesToProcess);
        initialVals = 0.5 * ones(1, length(artifactTypesToProcess));
        
        % Set initial values for provided types
        for k = 1:length(artifactTypesToProcess)
            artType = artifactTypesToProcess{k};
            if ismember(artType, providedTypes) && isfield(param.(artType), 'threshold')
                initialVals(k) = param.(artType).threshold;
            else
                showThresholdGUI = true;
            end
        end
        
        % If not all artifact types have thresholds, show GUI
        if length(providedTypes) < length(artifactTypesToProcess)
            showThresholdGUI = true;
        end
    end
    if showThresholdGUI
        % Pass the artifact types we're actually processing to the GUI
        result = svmThresholdGUI(initialVals, artifactTypesToProcess);
        if isempty(result)
            error('SVM threshold selection cancelled by user.');
        end
        param = struct();
        for k = 1:length(artifactTypesToProcess)
            if result.include(k)
                param.(artifactTypesToProcess{k}) = struct('threshold', result.thresholds(k));
            end
        end

        % Update artifactTypesToProcess to only include selected types
        selectedTypes = artifactTypesToProcess(result.include);
        if isempty(selectedTypes)
            error('No artifact types selected for classification.');
        end
        artifactTypesToProcess = selectedTypes;
    else
        % Update artifactTypesToProcess to only include types with parameters
        artifactTypesToProcess = fieldnames(param);
    end
    
    fprintf('Will process artifact types: %s\n', strjoin(artifactTypesToProcess, ', '));
end

fprintf(' ... initialization done\n');

% go through all signals
tic
fprintf('Auto-labelling signals (%d total) ----\n',length(signalIds))
for ii=1:N
    sigId=signalIds{ii};
    [curSignals ~]=interface.getSignalsById(sigId);
    [Nr,Nc] = size(curSignals);
    Nsec=ceil(Nc/samplingFreq);
    fprintf('   > signal %s (%d/%d), %ds',sigId, ii,N,Nsec)
        
    % Initialize annotation matrix with existing annotations if available
    allAn = false(Nr,Nsec,Nartif); 
    
    % If we have existing annotations, preserve them first
    if hasExistingAnnotations && ii <= length(existingAnnotations) && ~isempty(existingAnnotations{ii})
        existingAn = existingAnnotations{ii};
        % Make sure dimensions match, pad if necessary
        [eNr, eNsec, eNartif] = size(existingAn);
        if eNr == Nr && eNsec == Nsec && eNartif == Nartif
            allAn = existingAn;  % Use existing annotations as base
            fprintf(' (preserving existing annotations)');
        else
            fprintf(' (existing annotation dimensions mismatch - creating new)');
        end
    end
    
    if strcmpi(method, 'svm')
        % Use param and artifactTypesToProcess for classification
        curAn = sigInspectClassify(curSignals, sigId, samplingFreq, method, param);
        
        % Update only the processed artifact types in their original positions
        for a = 1:length(artifactTypesToProcess)
            artType = artifactTypesToProcess{a};
            % Find the position of this artifact type in the original list
            artIdx = find(strcmp(artifactTypes, artType));
            if ~isempty(artIdx)
                allAn(:,:,artIdx) = curAn(:,:,a);
            end
        end
        artSec = sum(any(any(curAn,3),2),1);
    else
        curAn = sigInspectClassify(curSignals, sigId, samplingFreq, method, varargin{:});
        allAn(:,:,artifactAutoWhich) = curAn;
        artSec = sum(any(curAn,3),2);
    end
    %curAn = sigInspectClassify(curSignals,samplingFreq, method,varargin{:});        
    %allAn(:,:,artifactAutoWhich)=curAn; % channel, second, artif. type
    
    % annot summary
    fprintf(' > %s artifact seconds in channels \n',sprintf('%d,',artSec))
    
    % store 
    annotation{ii} = allAn;
    % pad to artifact type count Nartif
    % annotation{ii} = padarray(allAn,[0 0 (Nartif-size(allAn,3))],0,'post');
        
end

fprintf('DONE: signals labelled in %.2f seconds----\n',toc)

% save annotation to *.mat file
annotationFile = '';
if(nargout<1)
    if(nargin<2 || isempty(pathToSave))
        pathToSave = sprintf('sigInspectAutoAnnotation%s.mat',datestr(now,'yyyy-mm-dd-HHMMSS'));
    end
    interfaceClass = class(interface);
    if strcmpi(method, 'svm')
        thresholds = struct();
        for k = 1:length(artifactTypesToProcess)
            art = artifactTypesToProcess{k};
            if isfield(param, art)
                thresholds.(art) = param.(art).threshold;
            end
        end
        save(pathToSave,'annotation','signalIds','artifactTypes','interfaceClass','thresholds','method')
    else
        save(pathToSave,'annotation','signalIds','artifactTypes','interfaceClass','method')
    end
    fprintf('ANNOTATION SAVED TO: %s\n',pathToSave)
    annotationFile = pathToSave;
else
    if(nargin<2 || isempty(pathToSave))
        pathToSave = sprintf('sigInspectAutoAnnotation%s.mat',datestr(now,'yyyy-mm-dd-HHMMSS'));
    end
    interfaceClass = class(interface);
    if strcmpi(method, 'svm')
        thresholds = struct();
        for k = 1:length(artifactTypesToProcess)
            art = artifactTypesToProcess{k};
            if isfield(param, art)
                thresholds.(art) = param.(art).threshold;
            end
        end
        save(pathToSave,'annotation','signalIds','artifactTypes','interfaceClass','thresholds','method')
    else
        save(pathToSave,'annotation','signalIds','artifactTypes','interfaceClass','method')
    end
    fprintf('ANNOTATION SAVED TO: %s\n',pathToSave)
    annotationFile = pathToSave;
end

% Return additional parameters for SVM method
if nargout > 2
    if strcmpi(method, 'svm')
        svmParams = param;
    else
        svmParams = [];
    end
end


%% ---- support functions: handling inputs & outputs

% handle input arguments - initialize signalIds
function interface = initInterface(vararg)
    % load signals using uiGetfile dialog
    interface=[];
    if(isempty(vararg))
         disp('sigInspectAutoLabel: select *.mat file with signal(s) in cell array variable signal, signals or data')         
         while(isempty(interface))
            % select signal using gui dialog
            [fileName, pathName, fileIndex] = uigetfile('*.mat','sigInspect: Select *.mat file with signal(s)');
            filePath = [pathName fileName];
            
            if(~ischar(fileName))
                h=errordlg(['sigInspectAutoLabel can not start without a signal - please select mat file with signals or provide signals or interface as an input parameter'],'can not start without signals', 'modal');
                waitfor(h)
                continue
            end
            
            try
                % load file, leave all settings at default values
                interface=sigInspectDataBasic(filePath);
            catch err
                % continue until a proper file is loaded
                bt=questdlg([err.message '> choose another file?'], err.identifier, 'Yes','Quit','Yes');
                switch(bt)
                    case 'Yes'
                        continue
                    case 'Quit'
                        return
                end
            end            
         end
    else
        % some input data provided - use interface or intialize default
        if(isa(vararg,'sigInspectDataInterface'))
            % user provided data interface
            interface=vararg;
        else
            % user hopefully provided signals as a structure
            interface=sigInspectDataBasic(vararg);
        end
    end
    