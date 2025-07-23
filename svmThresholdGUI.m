function result = svmThresholdGUI(initialVals, artifactNames, saveToFile, okCallback)
% svmThresholdGUI - GUI for setting SVM artifact thresholds (POW, BASE, FREQ)
% Returns a struct: result.thresholds (1x3), result.include (1x3 logical)

if nargin < 1 || isempty(initialVals)
    initialVals = [0.5 0.5 0.5];
end
if nargin < 2 || isempty(artifactNames)
    artifactNames = {'POW','BASE','FREQ'};
end

if nargin < 3
    saveToFile = false;
end

% Ensure vectors are the same length
numArtifacts = length(artifactNames);
if length(initialVals) ~= numArtifacts
    if isscalar(initialVals)
        initialVals = repmat(initialVals, 1, numArtifacts);
    else
        error('Length of initialVals (%d) must match length of artifactNames (%d)', ...
            length(initialVals), numArtifacts);
    end
end

% Calculate window height based on number of artifact types
windowHeight = 120 + 60 * numArtifacts;

f = figure('Name','SVM Artifact Thresholds','NumberTitle','off',...
    'MenuBar','none','ToolBar','none','Position',[500 500 450 windowHeight],...
    'Resize','off','WindowStyle','modal');

sliderHandles = gobjects(1,numArtifacts);
editHandles = gobjects(1,numArtifacts);
checkHandles = gobjects(1,numArtifacts);

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
    
    % Include checkbox
    checkHandles(i) = uicontrol('Style','checkbox','String','Include','Value',1,...
        'Position',[330 yPos-10 70 30],'FontSize',12);
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

% Initialize result variable
result = [];

% Wait for user interaction
uiwait(f);

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
        include = false(1, numArtifacts);
        
        for j = 1:numArtifacts
            thresholds(j) = sliderHandles(j).Value;
            include(j) = logical(checkHandles(j).Value);
        end
        
        % Check if at least one artifact type is selected
        if ~any(include)
            warndlg('Please select at least one artifact type to include.', 'No Types Selected', 'modal');
            return;
        end
        
        if saveToFile
            [file, path] = uiputfile('svm_thresholds.mat','Save thresholds as');
            if ischar(file)
                save(fullfile(path,file),'thresholds','include','artifactNames');
                fprintf('Thresholds saved to: %s\n', fullfile(path,file));
            end
        end
        
        result = struct('thresholds', thresholds, 'include', include, 'artifactNames', {artifactNames});
        
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
        result = [];
        uiresume(f);
        if ishandle(f)
            close(f);
        end
    end
end