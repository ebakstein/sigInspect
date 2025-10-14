function [annotation, annotationFile, params] = sigInspectAutoLabel(interfSignalOrPath, pathToSave, samplingFreq, method, varargin)
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

% Artifact types
if strcmpi(method,'svm')
    artifactTypes = {'POW','BASE','FREQ'};
elseif isprop(interface,'settings') && isfield(interface.settings,'ARTIFACT_TYPES')
    artifactTypes = interface.settings.ARTIFACT_TYPES;
else
    artifactTypes = {'ARTIF','UNSURE'};
end
Nartif = length(artifactTypes);
fprintf('Artifact types: %s\n', strjoin(artifactTypes, ', '));


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


% Placeholder for potential annotation copying logic - if the logic of
% matching artifact type names is implemented, TODO this part
% (Currently disabled — left here intentionally)

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

% init empty annotation array
annotation=cell(length(signalIds),1);

% params definition
if isempty(varargin)
    params = [];
else
    % if user passed a single struct as last arg, take it
    if numel(varargin) == 1
        params = varargin{1};
    else
        % otherwise follow original behaviour (varargin{:}) but safer
        params = varargin;
    end
end

% SVM threshold GUI logic (minimal: show GUI if no params passed)
if strcmpi(method, 'svm')
    if isempty(params) || ~isstruct(params)
        params = svmThresholdGUI();
        if isempty(params)
            return
        end
    end
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
    
    curAn = sigInspectClassify(curSignals, sigId, samplingFreq, method, params);
    % ensure curAn has 3 dimensions (Nch × Nsec × Nart)
    if size(curAn,3) > 1
        allAn(:,:,:) = curAn;
    else
        allAn(:,:,artifactAutoWhich) = curAn;
    end
    
    artSec = sum(any(curAn,3),2);
    % annot summary
    fprintf(' > %s artifact seconds in channels \n',sprintf('%d,',artSec))
    
    % store 
    annotation{ii} = allAn;
end

fprintf('DONE: signals labelled in %.2f seconds----\n',toc)

% save annotation 
% if(nargout<1)
if(nargin<2 || isempty(pathToSave))
    pathToSave = sprintf('sigInspectAutoAnnotation_%s_%s.mat',...
        lower(method), datestr(now,'yyyy-mm-dd-HHMMSS'));
end
interfaceClass = class(interface);
save(pathToSave,'annotation','signalIds','artifactTypes','interfaceClass','method', 'params')
fprintf('Saved annotation + params to: %s\n', pathToSave);
annotationFile = pathToSave;



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
    