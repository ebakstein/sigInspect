function varargout = sigInspectOverview(varargin)
% SIGINSPECTOVERVIEW MATLAB code for sigInspectOverview.fig
%      SIGINSPECTOVERVIEW, by itself, creates a new SIGINSPECTOVERVIEW or raises the existing
%      singleton*.
%
%      H = SIGINSPECTOVERVIEW returns the handle to a new SIGINSPECTOVERVIEW or the handle to
%      the existing singleton*.
%
%      SIGINSPECTOVERVIEW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SIGINSPECTOVERVIEW.M with the given input arguments.
%
%      SIGINSPECTOVERVIEW('Property','Value',...) creates a new SIGINSPECTOVERVIEW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before sigInspectOverview_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to sigInspectOverview_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help sigInspectOverview

% Last Modified by GUIDE v2.5 01-Apr-2019 17:07:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @sigInspectOverview_OpeningFcn, ...
                   'gui_OutputFcn',  @sigInspectOverview_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before sigInspectOverview is made visible.
function sigInspectOverview_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to sigInspectOverview (see VARARGIN)


% Choose default command line output for sigInspectOverview
handles.output = hObject;

if(isempty(varargin) || ~isstruct(varargin{1}) || ~isfield(varargin{1},'sigInspectMainWindow') || ~ishandle(varargin{1}.sigInspectMainWindow))
    % improper call
    erh=errordlg('sigInspectOverview GUI called with improper parameter - should be sigInspectOverview(parentHandles)','sigInspectOverview ERROR','modal');
    waitfor(erh);
    close(hObject);
%     warning('sigInspectOverview GUI called with improper parameter - should be sigInspectOverview(parentHandles)')
    return;
end
parentHandles = varargin{1};

% Set the click detection on the axis - does not work on axes child objects - lines, patches etc
%  set(handles.sigAxes,'buttondownfcn',{@axesClicked,handles}) 

 % set click callback for whole figure
 set(gcf,'WindowButtonDownFcn',{@figureClicked,handles});


% set annot mode selection change callback
set(handles.annotPanel, 'SelectionChangeFcn',@annotPlotSelection);

 
% x=(0:100)/100*4*pi;
% y=sin(x);
% plot(handles.sigAxes,x,y);

% update parent handles - add links to this figure
parentHandles.overviewFig = handles.overviewFigure;
parentHandles.overviewSigAxes = handles.sigAxes;

% restore window parameters from last session
if(isfield(parentHandles,'overviewPlotAnnotMode') && ~isempty(parentHandles.overviewPlotAnnotMode))
    set(handles.annotPanel, 'SelectedObject',handles.(['annotRadio' parentHandles.overviewPlotAnnotMode]));
end
if(isfield(parentHandles,'overviewPosition') && ~isempty(parentHandles.overviewPosition))
    set(handles.overviewFigure, 'Position',parentHandles.overviewPosition);
end
    

guidata(parentHandles.sigInspectMainWindow, parentHandles);
% uiresume(parentHandles.figure);

handles.parentFig = parentHandles.sigInspectMainWindow;
% Update handles structure
guidata(hObject, handles);
% uiwait();

% UIWAIT makes sigInspectOverview wait for user response (see UIRESUME)
% uiwait(handles.overviewFigure);
% 
% % handle clicking on the plot - i.e. second selection
% function axesClicked(hObject, eventdata,handles)
% 
% % get click coordinates
% crd = get(hObject,'currentpoint');
% sec=ceil(crd(1));
% % disp(crd(1,1:2));
% 
% % get handle to parent setSecond function, call it
% handles=guidata(handles.overviewFigure);
% parentHandles=guidata(handles.parentFig);
% cb=get(parentHandles.secondSelect,'Callback');
% 
% % set second
% set(parentHandles.secondSelect,'Value',sec)
% % call callback of parent figure
% cb(parentHandles.secondSelect,[]);

function figureClicked(hObject, eventdata,handles)

% get click coordinates
crd = get(hObject,'currentpoint');

% axis position
axPos = get(handles.sigAxes,'Position');

relPosX = (crd(1)-axPos(1))/axPos(3);
relPosY = (crd(2)-axPos(2))/axPos(4);

% clicked out of axes? > leave
inAx = relPosX>0 && relPosX<1 && relPosY>0 && relPosY<1;
if(~inAx)
    return;
end

% convert to axis coordinates
axLim = axis(handles.sigAxes);
pointX = relPosX*(axLim(2)-axLim(1))+axLim(1);
% pointY = relPosY*(axLim(4)-axLim(3))+axLim(3); % unused

sec=ceil(pointX);
% disp(crd(1,1:2));

% get handle to parent setSecond function, call it
handles=guidata(handles.overviewFigure);
parentHandles=guidata(handles.parentFig);
cb=get(parentHandles.secondSelect,'Callback');

% set second
set(parentHandles.secondSelect,'Value',sec)
% call callback of parent figure
cb(parentHandles.secondSelect,[]);



% --- Outputs from this function are returned to the command line.
function varargout = sigInspectOverview_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = []; %handles.output;


% --- Executes when user attempts to close overviewFigure.
function overviewFigure_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to overviewFigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% set properties of parent figure
if(~isempty(handles) && isfield(handles,'parentFig') && ishandle(handles.parentFig))
    parentHandles=guidata(handles.parentFig);
    parentHandles.overviewFig = []; % overview window figure handle
    parentHandles.overviewSigAxes = [];  % overview window signal axis handles
    set(parentHandles.overviewChck,'Value',0);
    parentHandles.overviewPosition = get(handles.overviewFigure,'Position');
    guidata(parentHandles.sigInspectMainWindow,parentHandles);
end

% Hint: delete(hObject) closes the figure
delete(hObject);


function annotPlotSelection(hObject, eventdata)

newLabel = get(eventdata.NewValue,'String');

% set new label to main window's handles
handles=guidata(hObject);
parentHandles = guidata(handles.parentFig);
parentHandles.overviewPlotAnnotMode = newLabel;
guidata(handles.parentFig, parentHandles);

% redrawOverview - from main window
parentHandles.redrawOverviewFun(parentHandles);


% --- Executes on key press with focus on overviewFigure and none of its controls.
function overviewFigure_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to overviewFigure (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
parentHandles = guidata(handles.parentFig);
parentHandles.keyPressFun(eventdata,parentHandles);


% --- Executes on key release with focus on overviewFigure and none of its controls.
function overviewFigure_KeyReleaseFcn(hObject, eventdata, handles)
% hObject    handle to overviewFigure (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was released, in lower case
%	Character: character interpretation of the key(s) that was released
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) released
% handles    structure with handles and user data (see GUIDATA)
parentHandles = guidata(handles.parentFig);
parentHandles.keyReleaseFun(eventdata,parentHandles);
    
