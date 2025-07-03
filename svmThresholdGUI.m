function result = svmThresholdGUI(initialVals, saveToFile, okCallback)
% svmThresholdGUI - GUI for setting SVM artifact thresholds (POW, BASE, FREQ)
% Returns a struct: result.thresholds (1x3), result.include (1x3 logical)

if nargin < 1 || isempty(initialVals)
    initialVals = [0.5 0.5 0.5];
end
if nargin < 2
    saveToFile = false;
end

f = figure('Name','SVM Artifact Thresholds','NumberTitle','off',...
    'MenuBar','none','ToolBar','none','Position',[500 500 450 250],...
    'Resize','off','WindowStyle','modal');

artifactNames = {'POW','BASE','FREQ'};
sliderHandles = gobjects(1,3);
editHandles = gobjects(1,3);
checkHandles = gobjects(1,3);

% First, create all controls
for i = 1:3
    uicontrol('Style','text','String',artifactNames{i},...
        'Position',[30 200-60*(i-1) 60 30],'FontSize',12,'HorizontalAlignment','left');
    sliderHandles(i) = uicontrol('Style','slider','Min',0,'Max',1,'Value',initialVals(i),...
        'Position',[100 210-60*(i-1) 150 20],'Tag',artifactNames{i});
    editHandles(i) = uicontrol('Style','edit','String',sprintf('%.2f',initialVals(i)),...
        'Position',[270 200-60*(i-1) 50 30],'FontSize',12);
    checkHandles(i) = uicontrol('Style','checkbox','String','Include','Value',1,...
        'Position',[330 200-60*(i-1) 70 30],'FontSize',12);
end

% Now set up the callbacks (after all handles are valid)
for i = 1:3
    sliderHandles(i).Callback = @(src,~) set(editHandles(i),'String',sprintf('%.2f',src.Value));
    editHandles(i).Callback = @(src,~) set(sliderHandles(i),'Value',str2double(src.String));
end

uicontrol('Style','pushbutton','String','OK','FontSize',12,...
    'Position',[80 20 80 40],'Callback',@okCallbackFcn);
uicontrol('Style','pushbutton','String','Cancel','FontSize',12,...
    'Position',[220 20 80 40],'Callback',@cancelCallback);

uiwait(f);

function okCallbackFcn(~,~)
    thresholds = zeros(1,3);
    include = false(1,3);
    for j = 1:3
        thresholds(j) = sliderHandles(j).Value;
        include(j) = logical(checkHandles(j).Value);
    end
    if saveToFile
        [file, path] = uiputfile('svm_thresholds.mat','Save thresholds as');
        if ischar(file)
            save(fullfile(path,file),'thresholds','include');
        end
    end
    result = struct('thresholds', thresholds, 'include', include);
    % Call the user-supplied callback
    if exist('okCallback','var') && ~isempty(okCallback) && isa(okCallback, 'function_handle')
        okCallback(result);
    end
    uiresume(f);
    close(f);
end

function cancelCallback(~,~)
    result = [];
    uiresume(f);
    close(f);
end
end