function [param, ok] = svmThresholdGUI(initialVals, artifactNames, saveToFile, okCallback)
% svmThresholdGUI - GUI for setting SVM artifact thresholds (POW, BASE, FREQ)
%
% OUTPUT:
%   param.(artifactType) = struct('threshold', selected_value);
%
% INPUTS:
%   initialVals  - Initial threshold values (scalar or 1×N vector)
%   artifactNames - Cell array of artifact type names
%   saveToFile   - Logical, whether to prompt saving thresholds
%   okCallback   - Optional function handle for user-defined actions on OK

if nargin < 1 || isempty(initialVals)
    initialVals = [0.5 0.5 0.5];
end
if nargin < 2 || isempty(artifactNames)
    artifactNames = {'POW','BASE','FREQ'};
end

if nargin < 3
    saveToFile = false;
end

if nargin < 4
    okCallback = []; 
end

% Ensure vectors are the same length
numArtifacts = numel(artifactNames);
if isscalar(initialVals)
    initialVals = repmat(initialVals, 1, numArtifacts);
elseif numel(initialVals) ~= numArtifacts
    error('Length of initialVals must match artifactNames.');
end

% Create GUI
windowHeight = 120 + 60 * numArtifacts;

f = figure('Name','SVM Artifact Thresholds','NumberTitle','off',...
    'MenuBar','none','ToolBar','none','Position',[500 500 350 windowHeight],...
    'Resize','off','WindowStyle','modal');

sliderHandles = gobjects(1,numArtifacts);
editHandles = gobjects(1,numArtifacts);
% checkHandles = gobjects(1,numArtifacts);

% Create controls for each artifact type
for i = 1:numArtifacts
    yPos = windowHeight - 60 - 60*(i-1);
    
    % Artifact name label
    uicontrol('Style','text','String',artifactNames{i},...
        'Position',[30 yPos-10 60 30],'FontSize',12,'HorizontalAlignment','left');
    
    % Threshold slider
    sliderHandles(i) = uicontrol('Style','slider','Min',0,'Max',1,'Value',initialVals(i),...
        'Position',[100 yPos 150 20],'Tag',artifactNames{i});
    
    % Threshold value edit box
    editHandles(i) = uicontrol('Style','edit','String',sprintf('%.2f',initialVals(i)),...
        'Position',[270 yPos-10 50 30],'FontSize',12);
    
    % % Include checkbox
    % checkHandles(i) = uicontrol('Style','checkbox','String','Include','Value',1,...
    %     'Position',[330 yPos-10 70 30],'FontSize',12);
end

% Set up callbacks after all handles are created
for i = 1:numArtifacts
    % Slider updates edit box
    sliderHandles(i).Callback = @(src,~) set(editHandles(i),'String',sprintf('%.2f',src.Value));
    
    % Edit box updates slider (with validation)
    editHandles(i).Callback = @(src,~) validateAndUpdateSlider(src, sliderHandles(i));
end

uicontrol('Style','pushbutton','String','OK','FontSize',12,...
    'Position',[80 20 80 40],'Callback',@okCallbackFcn);
uicontrol('Style','pushbutton','String','Cancel','FontSize',12,...
    'Position',[220 20 80 40],'Callback',@cancelCallback);

% Store result variable
param = struct();

% Wait for user interaction
uiwait(f);

if nargout > 1
    ok = ~isempty(fieldnames(param));  % true if user clicked OK, false if canceled
end


% ---------------------- Nested Functions ---------------------- %
    function validateAndUpdateSlider(editHandle, sliderHandle)
        % Validate edit box input and update corresponding slider
        value = str2double(editHandle.String);
        if isnan(value) || value < 0 || value > 1
            % Invalid input - revert to slider value
            editHandle.String = sprintf('%.2f', sliderHandle.Value);
            warndlg('Threshold must be between 0 and 1', 'Invalid Input', 'modal');
        else
            sliderHandle.Value = value;
        end
    end
    
    function okCallbackFcn(~,~)
        thresholds = zeros(1, numArtifacts);
        for k = 1:numArtifacts
            thresholds(k) = sliderHandles(k).Value;
        end
        
        % Convert to structured param format
        param = array2paramStruct(artifactNames, thresholds);

        if saveToFile
            [file, path] = uiputfile('svm_thresholds.mat','Save thresholds as');
            if ischar(file)
                save(fullfile(path,file),'param');
                fprintf('Thresholds saved to: %s\n', fullfile(path,file));
            end
        end
        
        result = struct('thresholds', thresholds, 'artifactNames', {artifactNames});
        
        % Call the user-supplied callback if provided
        if exist('okCallback','var') && ~isempty(okCallback) && isa(okCallback, 'function_handle')
            okCallback(result);
        end
        
        uiresume(f);
        if ishandle(f)
            close(f);
        end
    end
    
    function cancelCallback(~,~)
        param = struct();
        uiresume(f);
        if ishandle(f)
            close(f);
        end
        fprintf('SVM threshold selection canceled by user.')
    end

    function s = array2paramStruct(names, values)
        s = struct();
        for k = 1:numel(names)
            s.(names{k}) = struct('threshold', values(k));
        end
    end
end