function thresholds = svmThresholdGUI(initialVals, saveToFile, okCallback)
% svmThresholdGUI - GUI for setting SVM artifact thresholds (POW, BASE, FREQ)
% Usage:
%   thresholds = svmThresholdGUI([0.5 0.5 0.5]);
%   thresholds = svmThresholdGUI([0.5 0.5 0.5], true); % also saves to .mat
% Returns a 1x3 vector: [POW BASE FREQ]

if nargin < 1 || isempty(initialVals)
    initialVals = [0.5 0.5 0.5];
end
if nargin < 2
    saveToFile = false;
end

f = figure('Name','SVM Artifact Thresholds','NumberTitle','off',...
    'MenuBar','none','ToolBar','none','Position',[500 500 350 250],...
    'Resize','off','WindowStyle','modal');

artifactNames = {'POW','BASE','FREQ'};
sliderHandles = gobjects(1,3);
editHandles = gobjects(1,3);

% First, create all controls
for i = 1:3
    uicontrol('Style','text','String',artifactNames{i},...
        'Position',[30 200-60*(i-1) 60 30],'FontSize',12,'HorizontalAlignment','left');
    sliderHandles(i) = uicontrol('Style','slider','Min',0,'Max',1,'Value',initialVals(i),...
        'Position',[100 210-60*(i-1) 150 20],'Tag',artifactNames{i});
    editHandles(i) = uicontrol('Style','edit','String',sprintf('%.2f',initialVals(i)),...
        'Position',[270 200-60*(i-1) 50 30],'FontSize',12);
end

% Now set up the callbacks (after all handles are valid)
for i = 1:3
    sliderHandles(i).Callback = @(src,~) set(editHandles(i),'String',sprintf('%.2f',src.Value));
    editHandles(i).Callback = @(src,~) set(sliderHandles(i),'Value',str2double(src.String));
end

uicontrol('Style','pushbutton','String','OK','FontSize',12,...
    'Position',[60 20 80 40],'Callback',@okCallbackFcn);
uicontrol('Style','pushbutton','String','Cancel','FontSize',12,...
    'Position',[200 20 80 40],'Callback',@cancelCallback);

uiwait(f);

function okCallbackFcn(~,~)
        thresholds = zeros(1,3);
        for j = 1:3
            thresholds(j) = sliderHandles(j).Value;
        end
        if saveToFile
            [file, path] = uiputfile('svm_thresholds.mat','Save thresholds as');
            if ischar(file)
                save(fullfile(path,file),'thresholds');
            end
        end
        % Call the user-supplied callback
        if exist('okCallback','var') && ~isempty(okCallback) && isa(okCallback, 'function_handle')
            okCallback(thresholds);
        end
        uiresume(f);
        close(f);
    end

    function cancelCallback(~,~)
        thresholds = [];
        uiresume(f);
        close(f);
    end
end