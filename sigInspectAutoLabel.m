function annotation = sigInspectAutoLabel(interfSignalOrPath, pathToSave, samplingFreq, method, varargin)
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

% default artifact types / from interface
if(isprop(interface,'settings') && isfield(interface.settings,'ARTIFACT_TYPES'))
    artifactTypes =interface.settings.ARTIFACT_TYPES;
    fprintf('number_of_artifact_types = %d (from interface)\n',length(artifactTypes));
elseif strcmpi(method, 'svm')
    % Compare with the argument which artifact Types user had selected...
    artifactTypes = {'POW', 'BASE', 'FREQ'};
    fprintf('number_of_artifact_types = %d (DEFAULT: POW, BASE, FREQ)\n', length(artifactTypes));
else
    artifactTypes= {'ARTIF','UNSURE'};
    fprintf('number_of_artifact_types = %d (DEFAULT)\n',length(artifactTypes));
end
Nartif = length(artifactTypes);

% default field to store automatic artifact in
if(isprop(interface,'settings') && isfield(interface.settings,'ARTIFACT_AUTOLABEL_WHICH'))
    artifactAutoWhich =interface.settings.ARTIFACT_AUTOLABEL_WHICH;
    fprintf('auto artifact will be stored at position %d (%s)(from interface)\n',artifactAutoWhich,artifactTypes{artifactAutoWhich});
elseif(Nartif > 1)
    fprintf('saving to multiple artifact types')
else
    artifactAutoWhich = 1;
    fprintf('auto artifact will be stored at position %d (%s)(DEFAULT)\n',artifactAutoWhich,artifactTypes{artifactAutoWhich});
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
        % Only check the artifact types the user wants
        artifactTypes = intersect(fieldnames(param), {'POW','BASE','FREQ'});
        initialVals = [0.5 0.5 0.5];
        for k = 1:numel(artifactTypes)
            art = artifactTypes{k};
            idx = find(strcmp({'POW','BASE','FREQ'}, art));
            if isfield(param.(art), 'threshold')
                initialVals(idx) = param.(art).threshold;
            else
                showThresholdGUI = true;
            end
        end
        % If user didn't provide all three, you can decide if you want to show GUI or not
        if numel(artifactTypes) < 3
            showThresholdGUI = true;
        end
    end
    if showThresholdGUI
        result = svmThresholdGUI(initialVals);
        if isempty(result)
            error('SVM threshold selection cancelled by user.');
        end
        artifactNames = {'POW','BASE','FREQ'};
        param = struct();
        for k = 1:3
            if result.include(k)
                param.(artifactNames{k}) = struct('threshold', result.thresholds(k));
            end
        end
    end
    artifactTypes = fieldnames(param);  % Update artifact types from param keys
end

Nartif = length(artifactTypes); % <-- Set Nartif here, after artifactTypes is finalized

% After initializing interface and signalIds, compute loaderHash
loaderHash = getDataLoaderHash(interface);
disp(loaderHash);

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
        
    % init empty annot
    % annotation matrix - rows=channels, columns=seconds, slices=artifact types
    allAn = false(Nr,Nsec,Nartif); 
    if strcmpi(method, 'svm')
        % Use param and artifactTypes as set above
        curAn = sigInspectClassify(curSignals, sigId, samplingFreq, method, param, loaderHash);
        for a = 1:length(artifactTypes)
            allAn(:,:,a) = curAn(:,:,a);
        end
        artSec = sum(any(any(curAn,3),2),1);
    else
        curAn = sigInspectClassify(curSignals, sigId, samplingFreq, method, varargin{:}, loaderHash);
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
if(nargout<1)
    if(nargin<2 || isempty(pathToSave))
        pathToSave = sprintf('sigInspectAutoAnnotation%s.mat',datestr(now,'yyyy-mm-dd-HHMMSS'));
    end
    interfaceClass = class(interface);
    if strcmpi(method, 'svm')
        thresholds = struct();
        for k = 1:length(artifactTypes)
            art = artifactTypes{k};
            thresholds.(art) = param.(art).threshold;
        end
        save(pathToSave,'annotation','signalIds','artifactTypes','interfaceClass','thresholds')
    else
        save(pathToSave,'annotation','signalIds','artifactTypes','interfaceClass')
    end
    fprintf('ANNOTATION SAVED TO: %s\n',pathToSave)
end
    
% % bw-compatibility with old matlab
% if(nargout<1)
%     clear annot
% end


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

function loaderHash = getDataLoaderHash(interface)
    % Create a unique string from loader settings, file list, and file modification times
    % Works with different data loader types: sigInspectDataCsv, sigInspectDataBasic, sigInspectDataMatDir
    
    % Get settings if available
    if isprop(interface, 'settings')
        s = jsonencode(interface.settings);
        fprintf('\nSettings: %s\n', s);
    else
        s = '';
    end
    
    % Get file list and their modification times
    fileInfo = '';
    f = '';
    
    % Check for different possible file list properties
    if isprop(interface, 'fileList')
        fileList = interface.fileList;
        f = jsonencode(interface.fileList);
        fprintf('File list: %s\n', f);
    elseif isprop(interface, 'dataPath')
        % For sigInspectDataCsv, dataPath might contain the directory
        fileList = {interface.dataPath};
        f = jsonencode(interface.dataPath);
        fprintf('Data path: %s\n', f);
    elseif isprop(interface, 'pathOrMatrix')
        % For sigInspectDataBasic, might have pathOrMatrix
        if ischar(interface.pathOrMatrix)
            fileList = {interface.pathOrMatrix};
            f = jsonencode(interface.pathOrMatrix);
            fprintf('Path or matrix: %s\n', f);
        else
            fileList = {};
            f = '';
        end
    else
        fileList = {};
        f = '';
    end
    
    % Get file modification times if we have file paths
    if ~isempty(fileList)
        if iscell(fileList) || ischar(fileList)
            if ischar(fileList)
                fileList = {fileList};
            end
            for i = 1:length(fileList)
                if exist(fileList{i}, 'file')
                    fileDetails = dir(fileList{i});
                    if ~isempty(fileDetails)
                        % Include filename and modification time
                        fileInfo = [fileInfo, fileList{i}, ':', num2str(fileDetails(1).datenum), ';'];
                        fprintf('File %s modified at %s\n', fileList{i}, datestr(fileDetails(1).datenum));
                    end
                elseif exist(fileList{i}, 'dir')
                    % For directories, get modification time of the directory
                    fileDetails = dir(fileList{i});
                    if ~isempty(fileDetails)
                        fileInfo = [fileInfo, fileList{i}, ':', num2str(fileDetails(1).datenum), ';'];
                        fprintf('Directory %s modified at %s\n', fileList{i}, datestr(fileDetails(1).datenum));
                    end
                end
            end
        end
    end
    
    % Add data path explicitly if available
    if isprop(interface, 'dataPath')
        dataPath = interface.dataPath;
    else
        % Try to extract path from first file if available
        if ~isempty(fileList)
            if iscell(fileList)
                dataPath = fileparts(fileList{1});
            else
                dataPath = fileparts(fileList);
            end
        else
            dataPath = '';
        end
    end
    
    % Include class type in the hash to distinguish between different loader types
    classType = class(interface);
    
    % Combine all components in the hash
    loaderString = [s f fileInfo dataPath classType];
    fprintf('Loader string: %s\n', loaderString);
    loaderHash = string2hash(loaderString);

function hash = string2hash(str)
    % Simple hash function for a string (returns hex string)
    str = double(str);
    hash = dec2hex(mod(sum(str .* (1:numel(str))), 2^32), 8);
