function varargout = sigInspect(varargin)
% SIGINSPECT - a GUI for signal inspection and annotation
%
% USE:
% 1 - single multi-channel signal (matrix as input)
% C x N matrix (C = channels, N = samples )
%   sigInspect(signal, samplingFreq); % signal: chan. in rows, samples in columns
%       or
%   sigInspect(signal, samplingFreq, sAnnotTimes); % signal annotation times (a numeric vector)
%       or
%   sigInspect(signal, samplingFreq, sAnnotTimes, sAnnotTexts); % signal annotation times and texts
%       % (times as a numeric vector, texts as a cell array)
%
% 2 - multiple signals (cell array as input)
%   s={signal1,signal2,signal3};
%   sigInspect(s, samplingFreq);
%       or
%   sigInspect(s, samplingFreq, sAnnotTimes); % signal annotation times (a numeric vector)
%       or
%   sigInspect(s, samplingFreq, sAnnotTimes, sAnnotTexts); % signal annotation times and texts
%       % (times as a numeric vector, texts as a cell array)
%
% 3 - with pre-initialized data interface
%   intf = sigInspectDataCsv('/home/data-path/')
%   intf.settings.SAMPLING_FREQ = 8000;
%   intf.settings.ARTIFACT_TYPES = {'type 1','type 2','other'};
%   sigInspect(intf)
% 
%
% Find current version and PDF manual at https://github.com/ebakstein/sigInspect  
% 
% See also: SIGINSPECTAUTOLABEL
%
% Author: 
%   Eduard Bakstein, eduard.bakstein@felk.cvut.cz
%   and the neuro.felk.cvut.cz research group
%   

% Edit the above text to modify the response to help sigInspect

% Last Modified by GUIDE v2.5 20-Jan-2017 13:21:44


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @sigInspect_OpeningFcn, ...
                   'gui_OutputFcn',  @sigInspect_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);

 % initialization: user may provide 1-2 input params (filename, fs), during
 % initialization called with 4 input params
 if nargin > 2 && ischar(varargin{1}) % initialization (internal calls)
     gui_State.gui_Callback = str2func(varargin{1});
 else % user-prompted
    % check singleton
    h = findall(0,'tag','sigInspectMainWindow');   
    if(~isempty(h))
        eh = errordlg('Another running instance of sigInspect detected. Only a single copy of sigInspect is allowed, please, close the other instance first.','Another sigInspect running','modal');
        waitfor(eh);
        return  
    end
     
 end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before sigInspect is made visible.
function sigInspect_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to sigInspect (see VARARGIN)
global soundIsPlaying;
soundIsPlaying=false;
% Choose default command line output for sigInspect
handles.output = hObject;

handles.settings=struct();

% basic defaults
handles.settings.SAMPLING_FREQ = 24000;                                     % sampling frequency in Hz - default value 24kHz

% signal handling
handles.settings.DECIMATE_FACTOR = 1;                                       % decimate signals by this factor upon load

handles.settings.ADAPT_GAIN_TO_SIGNAL = 1;                                  % automatically adapt gain slider to signal
handles.settings.ADAPT_GAIN_QUANTILE = .002;                                % use this quantile for gain adaptation
handles.settings.ADAPT_GAIN_REFERENCE_AMPLITUDE = 20;                       % reference signal amplification at gain=1

handles.settings.NORMALIZE_SIGNAL_PER_CHANNEL=1;                            % normalize each channel (parallel signal) separately
handles.settings.NORMALIZE_SIGNAL_PER_CHANNEL_QUANTILE=0.1;                 % use this quantile for channel normalization

% overview window
handles.settings.OVERVIEW_DECIMATE_FACTOR = 20;                             % decimate signals by this factor before view in the overview window - this is multipliet with the basic DECIMATE_FACTOR
handles.settings.OVERVIEW_GAIN = 5;                                         % increase gain in overview window (to adapt for decimation)
handles.settings.OVERVIEW_ALWAYS_ON_TOP = 0;                                % keeps overview window always on top
handles.settings.OVERVIEW_CHANNEL_COLOR = [0,0.447,0.741];                  % signal color in overview window
handles.settings.OVERVIEW_ARTIFACT_COLOR='r';                               % color to plot artifacts (shown in overview only)

% automatic artifact type buttons - TODO:
handles.settings.ARTIFACT_TYPES={'ARTIF','UNSURE'};                         % cell array with artifact type abbreviations (max TODO types)
handles.settings.ARTIFACT_AUTOLABEL_WHICH=1;                                % which artifact type shall be use in sigInspectAutoLabel (just for compatibility)

% spectrogram
handles.settings.SPECTROGRAM_NFFT = 1024;                                   % points to be used in spectrogram FFT
handles.settings.SPECTROGRAM_FREQ_LIMS=[0 3000];                            % limit spectrogram to this frequency range
handles.settings.LINK_X_AXIS = 0;                                           % link x axis in spectrogram to signal view (useful for zoom)
handles.settings.DISABLE_SPECTROGRAM = 0;                                   % disable spectrogram (useful for long signals)
handles.settings.ENABLE_WHOLE_SPECTROGRAM = 0;                              % adds additional checkbox to toggle between current second and whole signal in spectrogram
handles.settings.ENABLE_HIDE_UNCHECKED = 1;                                 % enable checkbox to hide unchecked channels
handles.settings.WHOLE_SPECTROGRAM_SHOW_RECT = 1;                           % shows current second in whole spectrogram (works only if ENABLE_WHOLE_SPECTROGRAM=1)

% time series plot
handles.settings.SHOW_SIG_INFO = 1;                                         % show signal info from interface
handles.settings.PLOT_STEP = 150;                                           % distance on the y axis between channels
handles.settings.PLOT_CHANNELS = 5;                                         % max. number of prallel channels to plot (max. 9)
handles.settings.REVERSE_CHANNEL_ORDER = 1;                                 % 0: first channel at the top, last at the bottom 1: first channel at the bottom
handles.settings.CHANNELS_SAME_COLOR = 0;                                   % 0: use different color for each channel (see CHANNEL_COLORS) 1: use same color for all channels (first color from CHANNEL_COLORS)
handles.settings.CHANNEL_COLORS = {[0,0.447,0.741],...                      % cell array of channel colors in time series plot, if there are more channels than colors, colors are repeated. First color used if CHANNELS_SAME_COLOR=1
    [0.85,0.325,0.098],[0.929,0.694,0.125],[0.494,0.184,0.556],[0.466,0.674,0.188],[0.301,0.745,0.933],[0.635,0.078,0.184]};

% audio
handles.settings.START_WITH_SOUND_ON = 0;                                   % start with sound turned on
handles.settings.USE_BEEP = 0;                                              % beep when going to another signal
handles.settings.USE_AUDIOPLAYER = 1;                                       % use audio playback by audioplayer (matlab audioplayer instead of play())
handles.settings.PLAY_VIA_INTERNAL_PLAYER = 1;                                     % use more sophisticated audio
handles.settings.INTERNAL_PLAYER_DBG=0;                                            % show debug messages from INTERNAL_PLAYER (works only if USE_AUDIOPLAYER=1)
handles.settings.PLAY_SOUND_IN_SYNC = 0;
handles.settings.CHECK_CONST_SAMPLES=1;                                     % check signal for constant samples - visualize them in color


% annotation - loading checks
handles.settings.ANNOT_DEFAULT_FILENAME = 'sigInspectAnnot##.mat';          % default name for sigInspect annotation when save button is hit. ## is replaced by current date and time

% signal annotation - marking specific time positions with vertical lines and text
handles.settings.ANNOT_SIGNAL_TIMES = []; % a vector of times to annotate
handles.settings.ANNOT_SIGNAL_TEXTS = {}; % a cell array of texts to display at individual ANNOT_SIGNAL_TIMES (can be empty)

% debug settings
handles.settings.ANNOT_FILE_CHECK_INTERFACE_CLASS = 1;                      % should loaded annotation be checked for the same type of data interface?  DEBUG ONLY
handles.settings.ANNOT_FILE_CHECK_ARTIFACT_TYPES = 1;                       % should loaded annotation be checked for the same type of artifacts        DEBUG ONLY

% deprecated - to be removed
% handles.settings.SHOW_LOADING_LABEL=1;                                    % show "Loading" over the plot when loading signals

%-----------  TODO: ----------------------


%-------------------------  END TODO ---------------------------------



% ----------------- INITIALIZATION --------------------

handles = sigInspectInit(handles,varargin);

% Update handles structure
guidata(hObject, handles);

if(handles.quitNow)
    close(handles.sigInspectMainWindow);
end





% global figure initialization
function handles=sigInspectInit(handles, vararg)
    % ------ load interface from input arguments or data from file

    % addPath
    if(~exist('sigInspectDataBasic.m','file'))
        sigInspectAddpath;
    end
    
    % data interface
    handles.interface=[];
    handles.quitNow = 0;
    % -- loading data from artifactCatalog - transferred from catalogAnnotator
    if(isempty(vararg))
         disp('sigInspect: *.mat file with signal(s) in cell array variable signal, signals or data')         
         while(isempty(handles.interface))
            % select signal using gui dialog
            [fileName, pathName, fileIndex] = uigetfile('*.mat','sigInspect: Select *.mat file with signal(s)');
            filePath = [pathName fileName];
            
            if(~ischar(fileName))            
                h=errordlg(['sigInspect can not start without a signal - please select mat file with signals or provide signals or interface as an input parameter'],'can not start without signals', 'modal');
                waitfor(h)
                handles.quitNow=1;
                return
            end
            
            try
                % load file, leave all settings at default values
                handles.interface=sigInspectDataBasic(filePath);
            catch err
%                 % continue until a proper file is loaded
%                 h=errordlg([err.message '\n choose another file'], err.identifier, 'modal');
%                 waitfor(h)
%                 continue
                  % continue until a proper file is loaded
                bt=questdlg([err.message '> choose another file?'], err.identifier, 'Yes','Quit','Yes');
                switch(bt)
                    case 'Yes'
                        continue
                    case 'Quit'
                        handles.quitNow=1;
                        return
                end
 
            end            
         end

    else
        tmp = vararg{1};
        if(isa(tmp,'sigInspectDataInterface'))
            % user provided data interface
            handles.interface=tmp;
        else
            handles.interface=sigInspectDataBasic(tmp);
            if(length(vararg)>1 && ~isempty(vararg{2}) && isnumeric(vararg{2}))
                handles.interface.settings.SAMPLING_FREQ = vararg{2};
            end
            if(length(vararg)>2 && ~isempty(vararg{3}) && isnumeric(vararg{3}))
                handles.settings.ANNOT_SIGNAL_TIMES = vararg{3};
                if(length(vararg)>3 && ~isempty(vararg{4}) && iscell(vararg{4}))
                    handles.settings.ANNOT_SIGNAL_TEXTS = vararg{4};
                else
                    warning('unrecognized argument #4: a cell array expected for annotation texts');
                end
            else
                warning('unrecognized argument #3: a numeric vector expected for annotation times');
            end
        end
    end
    
    % TODO ---- -------------------------------------------------------------
    % automatic artifact button count
    % change initialization manual redraw to setSignal(handles,1,1)
    % END TODO ---------------------------------------------------------------
    
    if(isempty(handles.interface))
        h=errordlg('No signals nor interface loaded - quitting','sigInspect initialization error');
        waitfor(h);
        handles.quitNow = 1;
        return
    end

     % ------ load signal info from interface
    handles.signalIds=handles.interface.getSignalIds();
    handles.signalIds=handles.signalIds(:);
    handles.N = length(handles.signalIds);

    % check minimum signal count
    if(handles.N<1)
        dh=errordlg('Interface provided 0 signals - closing','sigInspect load error','modal');
        waitfor(dh);
        handles.quitNow=1;
        return
    end

    
    
    % internal state properties - connected with gui - 
    % DO NOT EDIT -----------
    handles.internal.INTERNAL_PLAYER_INITIALIZE=0;
    handles.internal.INTERNAL_PLAYER_REFRESH=1;
    handles.internal.INTERNAL_PLAYER_MARK_PLAYER_DESTROYABLE=2;
    handles.internal.INTERNAL_PLAYER_ADD_SIGNAL=3;
    handles.internal.INTERNAL_PLAYER_PLAY_SIGNAL=4;

    handles.internal.MaxChannels=10;
    handles.internal.MaxArtifactTypes=6;
    handles.internal.CheckboxLimitYPos=[.265 .947];
    handles.internal.CheckboxPos=[.96 NaN .034 .042]; % checkbox x-pos and size - for automatic generation of checkboxes
    % -----------------------

    
      
    % load settings from interface
    handles=copySettings(handles, handles.interface); 

    handles.artifactTypeN = length(handles.settings.ARTIFACT_TYPES);

    handles.annotation = cell(size(handles.signalIds));
    handles.seenSignals = false(size(handles.signalIds));
    
    % initial gain and other settings
    set(handles.gainSlider,'Value',1);
    set(handles.gainEdit,'Value', get(handles.gainSlider,'Value'));
    set(handles.gainEdit,'String', get(handles.gainSlider,'Value'));

    set(handles.thrSlider,'Value',30);
    set(handles.thrEdit,'Value', get(handles.thrSlider,'Value'));
    set(handles.thrEdit,'String', get(handles.thrSlider,'Value'));

    handles.anyChanges = false;
    handles.audioplayer = [];
    handles.curSignals = []; % array for storing current raw signals
    handles.curSigLenSec = 0;
    handles.curSigInfo=''; % current signal info, provided by interface
    handles.spectroWholePatch=[]; % handle to second-showing patch in whole spectrogram
    
    % color settings
    if(isempty(handles.settings.CHANNEL_COLORS) || ~iscell(handles.settings.CHANNEL_COLORS))
        warning('Invalid color format in handles.settings.CHANNEL_COLORS. Should be a nonempty cell array of color specifiers compatible with plot function (RGB triples or character shortcuts). Using default (blue)')
        handles.settings.CHANNEL_COLORS = {[0,0.447,0.741]};
    end
    if(handles.settings.CHANNELS_SAME_COLOR)
        handles.internal.channelColors = handles.settings.CHANNEL_COLORS(1);
    else
        handles.internal.channelColors = handles.settings.CHANNEL_COLORS(:);
    end
    %           repeat for all possible channels
    handles.internal.channelColors = repmat(handles.internal.channelColors,ceil(handles.internal.MaxChannels/length(handles.internal.channelColors)),1);
    handles.internal.channelColors = handles.internal.channelColors(1:handles.internal.MaxChannels);

    % for displaying time positions in the signals axes:
    % left mouse click followed by right mouse click results in displaying
    % the corresponding time interval in console
    handles.internal.lastButtonPressedInSignalAxes=0; % which mouse button was pressed when mouse in the signals axes?
    handles.internal.lastPositionInSignalAxes=0; % last mouse position in the signals axes
    
    % for overview window
    handles.overviewFig = []; % overview window figure handle
    handles.overviewSigAxes = [];  % overview window signal axis handles
    handles.curSignalsOverview = [];
    handles.overviewPlotAnnotMode = 'all';
    
    % function handles - for overview
    handles.redrawFun = @redraw;
    handles.redrawOverviewFun = @redrawOverview;
    handles.keyPressFun = @keyPressHandler;
    
    
    handles.samplingFreq = handles.settings.SAMPLING_FREQ/handles.settings.DECIMATE_FACTOR;
    handles.samplingFreqOverview = handles.samplingFreq/handles.settings.OVERVIEW_DECIMATE_FACTOR;

    set(handles.soundOnChck,'Value',handles.settings.START_WITH_SOUND_ON);
    
    % default threshold = 1/5 of y-step
    thr=handles.settings.PLOT_STEP/5;
    set(handles.thrSlider,'Value',thr);
    set(handles.thrEdit,'Value',thr);
    
    % threshold maximum - .5 * y-step
    thrMax=handles.settings.PLOT_STEP/3;
    set(handles.thrSlider,'Max',thrMax);
        
    % initialize spectrogram    
    redrawSpectrogram(handles,[],true);
    
    % initialize INTERNAL_PLAYER
    internalPlayer(handles,handles.internal.INTERNAL_PLAYER_INITIALIZE);
    
    handles=initArtifButtons(handles);
    handles=initChannelCheckboxes(handles);
    handles=initOverview(handles);

    % update id selector content
    redrawSignalSelect(handles);
    handles=setSignal(handles,1,1);
    handles = playSound(handles);
    
    % zoom callbacks
%     set(zoom(handles.signalAxes),'ActionPostCallback',@(obj,eventdata)showSecond(guidata(obj)));
%     set(zoom(handles.signalAxes),'ActionPostCallback',@(obj,eventdata)showSecond(guidata(obj)));
     set(zoom(handles.signalAxes),'ActionPostCallback',@postZoom);
    

    % check call without parameters
    nrgChck(nargout,'init')

% updates annotation buttons to match current artifact types
function handles = initArtifButtons(handles)
    for ai=1:handles.internal.MaxArtifactTypes
        btnName=sprintf('annot%dBtn',ai);
        if(ai>handles.artifactTypeN)
            delete(handles.(btnName));
%             set(handles.(btnName),'Visible','off')
        else
            set(handles.(btnName),'String',sprintf('F%d - %s',ai,handles.settings.ARTIFACT_TYPES{ai}))
            set(handles.(btnName),'Visible','on')
        end
    end

    nrgChck(nargout,'initArtifButtons')

function handles = initChannelCheckboxes(handles)
    sigAxPos = get(handles.signalAxes,'Position') ;%handles.internal.CheckboxLimitYPos;
    defaultPos = handles.internal.CheckboxPos;
    minMax=[sigAxPos(2) sigAxPos(2)+sigAxPos(4)]-defaultPos(4)/2;
    
    step = (minMax(2)-minMax(1))/(handles.settings.PLOT_CHANNELS+1);
    
    % add checkbox for each channel
    for ci=1:handles.settings.PLOT_CHANNELS
        chi=ci;
        if(~handles.settings.REVERSE_CHANNEL_ORDER)
            chi=handles.settings.PLOT_CHANNELS-ci+1;
        end
        
        % visible, used checkboxes
        chckName=sprintf('ch%dChck',ci);        
        % 
        ypos=minMax(1)+chi*step;
        pos=defaultPos;
        pos(2)=ypos;
        handles.(chckName)= uicontrol(handles.sigInspectMainWindow,...
                                  'Style','checkbox',...
                                  'Tag',chckName,...
                                  'Value',0,...
                                  'Units','normalized',...
                                  'Position',pos,...
                                  'Visible','on',...
                                  'String',num2str(ci),...
                                  'Callback',@chanChckCallback);
    end        

nrgChck(nargout,'initChannelCheckboxes')


% function voidCallBack

% checks settings field of interface, overwrites existing settings with
% those provided
function handles = copySettings(handles, interface)
    
    % check wheteher interface is of appropriate type
    if(~isa(interface,'sigInspectDataInterface'))
        error('sigInspect/addSettings: Not an object of class sigInspectDataInterface')
    end
    
    % copy all relevant settings fields (if any)
    fld=fieldnames(interface.settings);
    if(~isempty(fld))
        for fi=1:length(fld)            
            fName=fld{fi};
            
            % check if settings field exists
            if(~isfield(handles.settings, fName))
                warning('nonexistent settings field in interface: %s - skipping',fName)
                continue;
            end
            
            intVal=interface.settings.(fName);
            setVal=handles.settings.(fName);
  
            % check if new field is not empty
            if(isempty(intVal))
                warning('empty interface settings field: %s - skipping',fName)
                continue;
            end

            % check type of both fields
            if(~strcmp(class(setVal),class(intVal)))
                warning('existing settings field %s and new value in interface are of different types (%s vs %s) - skipping ',fName,class(setVal),class(intVal));
                continue;
            end

            % specific checks for individual settings
            switch(fName)
                case 'ARTIFACT_TYPES'
                    if(length(intVal) > handles.internal.MaxArtifactTypes)
                        warning('only at most %d artifact types are supported - trimming ARTIFACT_TYPES',handles.internal.MaxArtifactTypes)
                        intVal = intVal(1:handles.internal.MaxArtifactTypes);
                    end
                case 'PLOT_CHANNELS'
                    if(intVal > handles.internal.MaxChannels)
                        warning('only %d channels are supported - setting PLOT CHANNELS to %d',handles.internal.MaxChannels,handles.internal.MaxChannels)
                        intVal = handles.internal.MaxArtifactTypes;
                    end
            end

            % everything ok > copy
            handles.settings.(fName)=intVal;
        end
    end
    
% --- generates empty annotation for current signal
function annot = emptyCurAnnot(handles)
annot = false(getCurChan(handles),handles.curSigLenSec,handles.artifactTypeN); % annotation matrix - rows=channels, columns=seconds, slices=artifact types
        
% --- Outputs from this function are returned to the command line.
function varargout = sigInspect_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if(isfield(handles,'output'))
    varargout{1} = handles.output;
end

function doAdaptGainToSignal(handles)
    
    elCount = getCurChan(handles);

    qntl=[];
    for ci=1:elCount
        sig = handles.curSignals(ci,:);
        sig = sig-mean(sig);
        q=max(abs(quantile(sig,[handles.settings.ADAPT_GAIN_QUANTILE,1-handles.settings.ADAPT_GAIN_QUANTILE])));
        qntl=[qntl q];
    end
    qntl=median(qntl);
       
    %thr=get(handles.thrSlider,'Value');
    %gain=min(thr/qntl, 10);
    gain=min(handles.settings.ADAPT_GAIN_REFERENCE_AMPLITUDE/qntl, 10);
    
    set(handles.gainEdit,'String', gain);
    set(handles.gainEdit,'Value', gain);
    set(handles.gainSlider,'String', gain);
    set(handles.gainSlider,'Value', gain);

% redraws plot and figure
% params:
%  adaptGainToSignal:
%     0 - do not adapt gain according to signal
%     1 (the default) - do adapt gain according to signal, if enabled
%           by handles.settings.ADAPT_GAIN_TO_SIGNAL
%     2 - do adapt gain according to signal (even if not enabled 
%           by handles.settings.ADAPT_GAIN_TO_SIGNAL)
function redraw(handles,adaptGainToSignal)
    if nargin<2
        adaptGainToSignal=1;
    end
    if adaptGainToSignal<2
        adaptGainToSignal=handles.settings.ADAPT_GAIN_TO_SIGNAL*adaptGainToSignal;
    end

    axes(handles.signalAxes);
    
    % hide during redraw - useless: works only for background
%     set(handles.signalAxes,'Visible','off')
    
    cla() % clear signal axes
    sec=getCurSec(handles);
    xlim(sec+[-1 0])
    ylim([0,(1 + handles.settings.PLOT_CHANNELS)*handles.settings.PLOT_STEP])
    
    
    hold on    

    % adapt gain automatically?
    if adaptGainToSignal
        doAdaptGainToSignal(handles);
    end
            
    elCount=getCurChan(handles);
    selChans = getSelectedChannels(handles);    
    hideUnchecked = get(handles.hideChck, 'Value');
    
    for ci=1:elCount
                
        
        % plot all available channels        
        sig = handles.curSignals(ci,:);                
        t = linspace(0,(length(sig)-1)/handles.samplingFreq,length(sig));
        
        if(handles.settings.REVERSE_CHANNEL_ORDER)                        
            yShift = handles.settings.PLOT_STEP*ci;
        else
            yShift = handles.settings.PLOT_STEP*(handles.settings.PLOT_CHANNELS-ci+1);
        end
                    
        sig = sig-mean(sig);
        
        
        if(hideUnchecked && ~ismember(ci,selChans))
            % skip unchecked signals if selected
            line(t([1 end]),yShift*[1 1],'Color',.8*[1 1 1]);
            continue
        end        
        
        % amplify
        sig = sig*get(handles.gainEdit,'Value');        
        
        % plot        
        plot(t,sig+yShift,'color',handles.internal.channelColors{ci});
        
        % constant samples
        if(handles.settings.CHECK_CONST_SAMPLES)
            emptyInds= constSamples(sig,handles.curSigLenSec,0.01);
            if(any(emptyInds))
                plot(t(emptyInds),sig(emptyInds)+yShift,'r.','MarkerSize',1)
            end
        end
        
        % lines - threshold
        thr = get(handles.thrSlider,'Value');

        xlm = [0 handles.curSigLenSec];
        line(xlm,yShift+thr*[1 1],'Color','black')
        line(xlm,yShift-thr*[1 1],'Color','black')                
    end
    
%     xlim(sec+[-1 0])
%     ylim([0,(1 + handles.settings.PLOT_CHANNELS)*handles.settings.PLOT_STEP])


    if ~isempty(handles.settings.ANNOT_SIGNAL_TIMES)
        % display signal annotation
        tmpn=length(handles.settings.ANNOT_SIGNAL_TIMES); % # of times
        % a trick to make plotting fast: we plot multiple vertical lines as a single line polygon
        tmpx=reshape(handles.settings.ANNOT_SIGNAL_TIMES,1,tmpn); % reshape to form a row vector
        tmpx=[tmpx; tmpx; repmat(NaN,1,tmpn)]; % 3 rows: 1st: x position, 2nd: x position, 3rd: NaN
        tmpx=reshape(tmpx,prod(size(tmpx)),1); % reshape to column vector
        yl=ylim();
        yl(2)=yl(2)-diff(yl)/10; % annotation lines do not spand thw whole ylim()
        tmpy=repmat([yl NaN]',tmpn,1); % form a column vector
        plot(tmpx,tmpy,'color','k','LineStyle',':','LineWidth',1);

        if ~isempty(handles.settings.ANNOT_SIGNAL_TEXTS)
            tmpn=length(handles.settings.ANNOT_SIGNAL_TEXTS);
            yl=ylim();
            tmpy=yl(2)-diff(yl)/20; % the text appear above the lines
            tmpy=repmat(tmpy,1,tmpn);
            text(handles.settings.ANNOT_SIGNAL_TIMES(1:tmpn),tmpy,handles.settings.ANNOT_SIGNAL_TEXTS,'HorizontalAlignment','center');
        end
    end

    
    hold off
%     set(handles.signalAxes,'Visible','on')
    
%     tic
    redrawSpectrogram(handles);
%    toc    
    
    if(handles.settings.LINK_X_AXIS)
        linkaxes([handles.signalAxes,handles.spectroAxes],'x')
    end
        
    redrawSignalSelect(handles)
        
    % show true sig info from DAO
    showSigInfo(handles)
    
    % redraw overview window
    redrawOverview(handles)

% redraw plot in overview figure    
function redrawOverview(handles, justPatch) 
    overviewOn=get(handles.overviewChck,'Value')==1;    
    persistent patchHandle;
    if(overviewOn && ~isempty(handles.overviewFig)) 
        if(nargin<2)
            justPatch=0;
        end
        
%         disp('DBG: redrawOverview')
        % switched on and overviewFig initialized
        if(~isempty(handles.curSignalsOverview) && ~isempty(handles.overviewSigAxes))       
            if(~justPatch)
                patchHandle=[];
                cla(handles.overviewSigAxes);
                hold(handles.overviewSigAxes,'on');
                [Nch N]= size(handles.curSignalsOverview);
                t=(1:N)/handles.samplingFreqOverview;
                for ii=1:Nch
                   sig=handles.curSignalsOverview(ii,:);
                       
                   if(handles.settings.REVERSE_CHANNEL_ORDER)                        
                     yShift = handles.settings.PLOT_STEP*ii;
                   else
                     yShift = handles.settings.PLOT_STEP*(Nch-ii+1);
                   end
%                    yShift = handles.settings.PLOT_STEP*ii;
                   sig = sig-mean(sig);
                   % amplify
                   sig = sig*get(handles.gainEdit,'Value')*handles.settings.OVERVIEW_GAIN;        
                   % plot
                   plot(handles.overviewSigAxes,t,sig+yShift,'color',handles.settings.OVERVIEW_CHANNEL_COLOR);

                   % overplot annotation
                   plotInds = getOverviewAnnotInds(handles, ii);
                   if(any(plotInds))
                       sigTmp = sig;
                       sigTmp(~plotInds)= nan;
                       plot(handles.overviewSigAxes,t,sigTmp+yShift,'Color',handles.settings.OVERVIEW_ARTIFACT_COLOR)
                   end
                end                  
                ylim(handles.overviewSigAxes,[0 (Nch+1)*handles.settings.PLOT_STEP]);
                % plot signal edge-to-edge
                xlim(handles.overviewSigAxes,[t(1) t(end)])
            end 
            % patch - current second
            sec=getCurSec(handles);
            ylm=ylim();
            
            if(~isempty(patchHandle))
                delete(patchHandle);
            end
            patchHandle=patch([sec-1 sec sec sec-1],[ylm(1) ylm(1) ylm(2) ylm(2)],'r','FaceColor','r','FaceAlpha',.15,'EdgeColor','r','Parent',handles.overviewSigAxes);
            
            
            
            % return control to main figure
%             figure(handles.sigInspectMainWindow)
            
        else
            warning('curSignalsOverview not initialized')
        end        
    end

% indices of artifact signal samples to plot - based on overview annot
% checkbox
function plotInds = getOverviewAnnotInds(handles,chan)

    sigN = size(handles.curSignalsOverview,2);
    plotInds = false(1,sigN);    
    % no overview figure present, return zeros - no artif to plot
    if(isempty(handles.overviewFig) || isempty(handles.annotation{getCurSig(handles)})) 
        return;
    end
    
    % check type correctness
    if(~any(strcmp(handles.overviewPlotAnnotMode,{'none','all','maxCleanSegment'})))
        error('unknown annotation display method in sigOverview: %s',handles.overviewPlotAnnotMode);
    end

    
    annot = handles.annotation{getCurSig(handles)};
    annot = sum(annot(chan,:,:),3);
    
    if(strcmp(handles.overviewPlotAnnotMode,'none') || all(annot==0))
        return;
    end
    
    if(strcmp(handles.overviewPlotAnnotMode,'maxCleanSegment'))
        annot = ~maxNonZeroSegment(~annot);
    end
    
    anT = (1:size(annot,2))-.5;
    sigT=(0:sigN-1)/handles.samplingFreqOverview;
    plotInds = logical(interp1(anT,double(annot),sigT,'nearest','extrap'));
    

    
    % interpolate to signal dimensions
    
function ids=maxNonZeroSegment(v)
    % ids=maxNonZeroSegment(v)
    %   returns logical indices of the longest nonzero (and non-nan) 
    %   segment in vector v
    nz=~isnan(v) & v~=0;    % nonzero inds
    tr=diff([0 nz 0]);      % transitions

    trSt=find(tr==1); % add start and end
    trEn=find(tr==-1)-1; % add start and end

    trLn = trEn-trSt+1; % length

    [~, id] = max(trLn);

    startseg = trSt(id);
    endseg   = trEn(id); 

    ids=false(size(v));
    ids(startseg:endseg)=1;

% initialize overview figure + handles
function handles = initOverview(handles)
    % is the checkbox on?
    isOn=get(handles.overviewChck,'Value')==1;
    if(~isOn)
        % overview switched off - close if open
        if(~isempty(handles.overviewFig))
            
            close(handles.overviewFig)
            handles.overviewFig=[];
            
%             if(handles.settings.OVERVIEW_ALWAYS_ON_TOP)
%                 set(handles.sigInspectMainWindow, 'WindowButtonMotionFcn', '')
%             end

        end
    else
        % initialize overview window - will write its handles to 
%         uiwait()
        % uiresume called by sigInspectOverview initialization
        sigInspectOverview(handles);
        
        % changes sigInspect handles! - sets handles.overviewFig, handles.overviewSigAxes        
        handles=guidata(handles.sigInspectMainWindow);        
        if(handles.settings.OVERVIEW_ALWAYS_ON_TOP)
            %set(handles.sigInspectMainWindow, 'WindowButtonMotionFcn', @overviewOnTop)
            WinOnTop(handles.overviewFig);
        end    
    end
    
    nrgChck(nargout,'initOverview')

% function overviewOnTop(hObject,varargin)
% handles=guidata(hObject);
% if(~isempty(handles.overviewFig))
%     set(handles.overviewFig,'Visible','on')
% end

    
    
function initSpectrogram(handles)
    % sets proper visibility of spectrogram controls
    if(handles.settings.DISABLE_SPECTROGRAM)
        set(handles.spectroAxes,'Visible','off');
        set(handles.spectroChck,'Visible','off');        
        
        % resize signal window to use the free space
        spPos=get(handles.spectroAxes,'Position');
        sigPos=get(handles.signalAxes,'Position');
        sigPos(2) = spPos(2)+.02; sigPos(4) = sigPos(4)+spPos(4)-.02;
        set(handles.signalAxes,'Position',sigPos);
    end
    
    if(~handles.settings.ENABLE_WHOLE_SPECTROGRAM || handles.settings.DISABLE_SPECTROGRAM)
        % hide whole spectogram toggle
        set(handles.spectroWholeChck,'Visible','off')        
    end        
    if(~handles.settings.ENABLE_HIDE_UNCHECKED)
        set(handles.hideChck,'Visible','off')
    end


function redrawSpectrogram(handles,forceRedraw,initialize)
    %redrawSpectrogram(handles) - redraws spectrogram plot (if only 1
    %channel is selected), or clears the plot(otherwise)    

    persistent chanN;
    persistent sigN;
    
    % initialization only
    if(nargin>2&&initialize)
        % initialize persistent vars. + hide spectrogram
        chanN=nan;
        sigN=nan;
        initSpectrogram(handles);     
        return
    end
    
    % spectrogram disabled        
    if handles.settings.DISABLE_SPECTROGRAM        
        return                       
    end

    
    if(nargin<2 || isempty(forceRedraw))
        forceRedraw=false;
    end
 
    
    ah=handles.spectroAxes;
    th=handles.spectroTxt;
    sh=handles.signalAxes;
    
    ch=getSelectedChannels(handles);
    sn=getCurSig(handles);
    sc=getCurSec(handles);
    
    % do redraw
    if(~isempty(ch) && get(handles.spectroChck,'Value')) % at least one sig selected + spectro on

        % get signal
        sig=handles.curSignals(ch(1),:);
        for ii=ch(2:end)
            sig=sig+handles.curSignals(ii,:);
        end
        
%         % crop only the respective part
%         sec=getCurSec(handles);
%         sig=sig((sec-1)*handles.settings.SAMPLING_FREQ+1 : sec*handles.settings.SAMPLING_FREQ);
         %axes(ah);
%          set(ah,'Visible','off');
         if(length(chanN)~=length(ch) || any(chanN~=ch) || sigN~=sn || forceRedraw)
            % plotting itself (if source has changed - not just threshold etc)
            NFFT = handles.settings.SPECTROGRAM_NFFT;
            nOverlap = round(3*NFFT/4);

            % extend the signal in order to make the spectrogram to
            % span the whole time range
            nOverlapHalf = ceil(nOverlap/2);
            nFilling = length(sig)-NFFT*floor(length(sig)/NFFT);
            sig = [fliplr(sig(1:nOverlapHalf)) sig fliplr(sig(end-nOverlap-nFilling:end))];
           
            hold off;
            [~,F,T,P] = spectrogram(sig,NFFT,nOverlap,NFFT,handles.samplingFreq);
            T = T - T(1);

    %         T = T*1000; % time to ms
            imagesc(T,F,log(abs(P)),'Parent',ah);
            set(ah,'ydir','normal')
            set(th,'String',num2str(ch));
            %set(ah,'Visible','off');
            ylim(ah,handles.settings.SPECTROGRAM_FREQ_LIMS)

            
            chanN=ch;
            sigN=sn;
        end
        % show x axis on the top
        %set(ah,'Box','off','XAxisLocation','top');
        set(ah,'Box','off','XAxisLocation','top','XTick',[]);
        
        % hide signal xtic and x axis
        %set(sh,'XTick',[]);%,'XTickLabel',[]);
            
%         xlabel('time [ms]')
%         ylabel('frequency [Hz]')
%         title('spectrogram')
        % x axis - depends on whether wholeSpectrogram is on
        if(handles.settings.ENABLE_WHOLE_SPECTROGRAM && get(handles.spectroWholeChck,'Value')==1)
            xlim(ah,[0 handles.curSigLenSec]);
        else
            xlim(ah,sc+[-1 0]);
        end
%         set(ah,'Visible','on');

    else
        % clear spectro axes
        cla(ah);
        set(ah,'Visible','off');
        set(th,'String','');
        
        % show x axis in signal plot
        set(sh,'XTickMode','auto','XTickLabelMode','auto');
    end
    
function showSigInfo(handles)
        
        hndl = handles.sigInfoTxt;
                
        % toggle visibility
        if(handles.settings.SHOW_SIG_INFO)
            set(hndl,'Visible','on')
        else
            set(hndl,'Visible','off');
            return;
        end
        
%         row = getCurSig(handles);
%         an = handles.annotation{row};
%         
%         els = '';
%         for ii=1:getCurChan()
%             els = [els num2str(ii) ':' an.electrodes{ii}(1:3) ':' an.areas{ii} ' - '];
%         end
%         els = els(1:end-3);
%          str = sprintf('#%d-%s (%s)',an.patientId,an.position,els);
         sigId=handles.signalIds{getCurSig(handles)};
         str = sprintf('%s: %s',sigId,handles.curSigInfo);              
         
         set(hndl,'String',str)

% --- Sets second selector values to 1: length of current signal    
function handles = redrawSecondSelect(handles)
set(handles.secondSelect,'string',arrayfun(@num2str,1:handles.curSigLenSec,'uniform',0))
nrgChck(nargout,'redrawSecondSelect')


function redrawSignalSelect(handles)
    % redraws signal select according to annotation
    str = {};
    for ii = 1:handles.N        
        an = handles.annotation{ii};
        if(~isempty(an) && any(an(:)))
            ast = '*';        
        elseif(handles.seenSignals(ii))
            % seen signals annotated by a +
            % '+':
            ast = '.';
        else
            ast='';
        end
        
        str{ii} = sprintf('%s%s',ast,handles.signalIds{ii});
    end    
    set(handles.signalSelect,'String',char(str));

function showSecond(handles, doRedrawOverview)        
    persistent spectroWholePatch;
    
    sec = getCurSec(handles);
    xlim(handles.signalAxes,sec+[-1 0]);
    
    if(~isempty(spectroWholePatch))
        delete(spectroWholePatch);
        spectroWholePatch=[];
    end
    
    if(~(handles.settings.ENABLE_WHOLE_SPECTROGRAM && get(handles.spectroWholeChck,'Value')))
        xlim(handles.spectroAxes,sec+[-1 0]);
    elseif(get(handles.spectroChck,'Value'))
        xlim(handles.spectroAxes,[0 handles.curSigLenSec]);
        if(handles.settings.WHOLE_SPECTROGRAM_SHOW_RECT)
            ylm=ylim(handles.spectroAxes);
            spectroWholePatch=patch([sec-1 sec sec sec-1],[ylm(1) ylm(1) ylm(2) ylm(2)],'r','FaceColor','b','FaceAlpha',.15,'EdgeColor','b','Parent',handles.spectroAxes);
        end
    end    
    
    if(nargin<2 || doRedrawOverview)
        redrawOverview(handles,true);
    end
    dispAnnotation(handles);
%     playSound(handles);

function internalPlayer_stop_fun(source,event)
    userdata=get(source,'UserData');
    handles=userdata{2};
    if handles.settings.INTERNAL_PLAYER_DBG, fprintf(1,'>>> player %d is finished\n',userdata{1});end
    internalPlayer(handles,handles.internal.INTERNAL_PLAYER_MARK_PLAYER_DESTROYABLE,userdata{1});
    % unfortunatelly, can't start any new players from within here

function rv = internalPlayer(handles,cmd,arg,arg2,arg3)
    persistent playerId;
    persistent playerQueue;
    persistent signalQueue;
    persistent destroyablePlayerId;

    if nargin<5
        arg3=[];
        if nargin<4
            arg2=[];
            if nargin<3
                arg=[];
            end
        end
    end
    rv=[];
    switch cmd
        case handles.internal.INTERNAL_PLAYER_INITIALIZE
            playerId=0;
            playerQueue=[];
            signalQueue=[];
            destroyablePlayerId=[];
        case handles.internal.INTERNAL_PLAYER_REFRESH
            if handles.settings.INTERNAL_PLAYER_DBG, fprintf(1,'internalPlayer: refreshing\n');end
            if handles.settings.INTERNAL_PLAYER_DBG, fprintf(1,'  playerQueue size: %d\n',length(playerQueue));end
            if handles.settings.INTERNAL_PLAYER_DBG, fprintf(1,'  signalQueue size: %d\n',length(signalQueue));end
            % remove any finished players from the queue
            if ~isempty(destroyablePlayerId)
                p=playerQueue(1);
                tmp=get(p,'UserData');
                currentPlayerId=tmp{1};
                if destroyablePlayerId~=currentPlayerId
                    error('the current player is NOT detroyable');
                end
                destroyablePlayerId=[];
                if ~isplaying(p)
                    if handles.settings.INTERNAL_PLAYER_DBG, fprintf(1,'  removing unused player %d\n',currentPlayerId);end
                    p=[];
                    playerQueue=playerQueue(2:end);
                else
                    error('SHOULD NOT HAPPEN!!!');
                    %break
                end
            end
            % start a new player for next signal, if any
            if isempty(playerQueue) && ~isempty(signalQueue)
                s=signalQueue{1}{1};
                fs=signalQueue{1}{2};
                info=signalQueue{1}{3};
                if handles.settings.INTERNAL_PLAYER_DBG, fprintf(1,'  playing next signal (%s) via player %d\n',info,playerId);end
                p=audioplayer(s,fs);
                set(p,'StopFcn',@internalPlayer_stop_fun);
                set(p,'UserData',{playerId,handles});
                rv=playerId;
                playerId=playerId+1;
                signalQueue=signalQueue(2:end);
                play(p);
                %tic
                playerQueue=p;
            end
        case handles.internal.INTERNAL_PLAYER_MARK_PLAYER_DESTROYABLE
            if handles.settings.INTERNAL_PLAYER_DBG, fprintf(1,'internalPlayer: player %d finished\n',arg);end
            destroyablePlayerId=arg;
            %toc
        case handles.internal.INTERNAL_PLAYER_ADD_SIGNAL
            if handles.settings.INTERNAL_PLAYER_DBG, fprintf(1,'internalPlayer: attempting to queue a new signal\n');end
            if handles.settings.INTERNAL_PLAYER_DBG, fprintf(1,'  playerQueue size: %d\n',length(playerQueue));end
            if handles.settings.INTERNAL_PLAYER_DBG, fprintf(1,'  signalQueue size: %d\n',length(signalQueue));end
            if length(signalQueue)>0
                if handles.settings.INTERNAL_PLAYER_DBG, fprintf(1,'!!! queue not empty, removing older signals\n');end
                signalQueue=[];
            end
            signalQueue=[signalQueue {{arg,arg2,arg3}}];
        case handles.internal.INTERNAL_PLAYER_PLAY_SIGNAL
            if handles.settings.INTERNAL_PLAYER_DBG, fprintf(1,'internalPlayer: play signal\n');end
            internalPlayer(handles,handles.internal.INTERNAL_PLAYER_ADD_SIGNAL,arg,arg2,arg3);
            internalPlayer(handles,handles.internal.INTERNAL_PLAYER_REFRESH);
        otherwise
            error(['internalPlayer: unknown command ' num2str(cmd)]);
    end

function handles = playSound(handles)
    global soundIsPlaying;
    %TODO play sound of filtered selected signal here
%     playCh = str2num(get(get(handles.chSelect,'SelectedObject'),'String'));
    

    if(~getSoundState(handles))
        % sound off
        return
    end
    if(~soundIsPlaying) % is used only when PLAY_SOUND_IN_SYNC is 1
    sig = getMixedFiltSignal(handles);

    if(any(sig))
        % play if filtered signal not only 0
%         
        if(handles.settings.USE_AUDIOPLAYER)
            if handles.settings.PLAY_VIA_INTERNAL_PLAYER
                internalPlayer(handles,handles.internal.INTERNAL_PLAYER_PLAY_SIGNAL,sig,handles.samplingFreq,[num2str(getCurSig(handles)) ' ' num2str(getCurSec(handles))]);
            else
                %stopSound(handles);
                handles.audioplayer = audioplayer(sig, handles.samplingFreq);
                if(handles.settings.PLAY_SOUND_IN_SYNC)
                   % FUJ FUJ HACK, but it works
                   disp('Playing sound');
                   soundIsPlaying=true;
                   playblocking(handles.audioplayer);
                   soundIsPlaying=false;
                    disp('Playing ended');
                else
                    play(handles.audioplayer);
                end
                % guidata(handles.sigInspectMainWindow, handles);
            end
        else
            soundsc(sig, handles.samplingFreq);  
        end
    end
     
    nrgChck(nargout,'playSound');
    
    end
    

function stopSound(handles)
if(~isempty(handles.audioplayer) && isplaying(handles.audioplayer))
    if(handles.settings.PLAY_SOUND_IN_SYNC==0)
        stop(handles.audioplayer);
    end
end

% change state of overview window (e.g. by keypress)
function handles=toggleOverview(handles,state)
if(nargin>1)
    isOn=state;
else
    isOn=~get(handles.overviewChck,'Value');
end 
set(handles.overviewChck,'Value',isOn);

handles=initOverview(handles);
if(isOn)
    handles=computeOverviewSignals(handles);
    guidata(handles.sigInspectMainWindow,handles);
    redrawOverview(handles);
end
       
function handles=toggleSound(handles)
    % stop sound (if playing)
    st = getSoundState(handles);

    st = ~st;

    set(handles.soundOnChck,'Value',st);

    if(~st)
        stopSound(handles)
    else
        handles=playSound(handles);
    end
    
   nrgChck(nargout,'toggleSound');

function toggleSpectrogram(handles)
    % switch spectrogram on/off (for keyboard)
    
    hndl=handles.spectroChck;
    st=~get(hndl,'Value');            
    
    if(strcmp(get(hndl,'Enable'),'on')) %only if enabled

        set(hndl,'Value',st);

        % redraw (plot/clear if necessary)
        ch=getSelectedChannels(handles);
        if(length(ch)>0) % at least 1 ch selected (avoid unnecessary redraws)
            redrawSpectrogram(handles,true);
%             showSecond(handles);
        end
    end
    
function toggleHideUnselected(handles, hide)
    % toggle between "show all channels" and "hide unchecked channels"
    if(~handles.settings.ENABLE_HIDE_UNCHECKED)
        return;
    end
    
    if(nargin<2)
        % invert state
        hide = ~ get(handles.hideChck,'Value')==1;
    end
    set(handles.hideChck,'Value',hide)    
    redraw(handles,0);

    
function toggleSpectrogramWhole(handles,whole)
    % Toggle between spectrogram of whole signal and current second
    % works only if ENABLE_WHOLE_SPECTROGRAM is set to 1
    if(~handles.settings.ENABLE_WHOLE_SPECTROGRAM)
        return;
    end    
        
    if(nargin<2)
        % invert state
        whole = ~ get(handles.spectroWholeChck,'Value')==1;
        set(handles.spectroWholeChck,'Value',whole)
    end
    if(whole)
%         xlm=[0 handles.curSigLenSec];
        set(handles.spectroWholeChck,'BackgroundColor',[1 0 0]);
    else
%         xlm=getCurSec(handles)+[-1 0];
        set(handles.spectroWholeChck,'BackgroundColor',get(handles.spectroChck,'BackgroundColor'));
    end     
%     xlim(handles.spectroAxes,xlm);
    showSecond(handles);
%     figure(handles.sigInspectMainWindow)



function enableDisableSpectroChck(handles,noSpectRedraw)
    %         enable/disable spectrogram checkbox based on number of selected channels
    if(nargin<2)
        noSpectRedraw=false;
    end
    
    forceRedraw=false;
    
    hndl=handles.spectroChck;
    on=strcmp(get(hndl,'Enable'),'on');            

    nch=length(getSelectedChannels(handles));
    
    if(nch>0 && ~on)
        set(hndl,'Enable','on');
        forceRedraw=true;
    elseif(nch==0 && on)
        set(hndl,'Enable','off');        
    end
    
    if(~noSpectRedraw)
        redrawSpectrogram(handles,forceRedraw);
    end


function handles=playBeep(handles)

    if(handles.settings.USE_BEEP && handles.settings.PLAY_SOUND_IN_SYNC==0)  
        if(handles.settings.USE_AUDIOPLAYER)
            % synthesize beep
            Fs=1000;toneFreq=350;nSeconds=.05;
            sig = 0.5*sin(linspace(0, nSeconds*toneFreq*2*pi, round(nSeconds*Fs)));
%             sig = sig.*hamming(length(sig))'; % make it smooth
            
            if handles.settings.PLAY_VIA_INTERNAL_PLAYER
                internalPlayer(handles,handles.internal.INTERNAL_PLAYER_PLAY_SIGNAL,sig,Fs);
                pause(.2);
            else
            
                % make sure no sound is playing
                stopSound(handles);
            
                %play it
                handles.audioplayer = audioplayer(sig, Fs);
                play(handles.audioplayer);
            end
        else
            beep();
        end
    end
  nrgChck(nargout,'playBeep');
  
function nrgChck(nrg,funName)
    % if nargout==0, returns warning
    if(nrg==0)
        warning('%s should be called with output parameter handles!',funName)
    end
  
function sig = getMixedFiltSignal(handles)    
    % returns mixed channels filtered for sound playing - channels
    % identified from channel selectors

    % signals to play
    channels = getSelectedChannels(handles);

    sigLen = handles.samplingFreq*1;
    sig = zeros(1,sigLen);
    for ii = 1:length(channels)
        curSig = getFiltSignal(handles,channels(ii));
        curSig = curSig(1:min(length(curSig),sigLen));
        % add to signal        
        sig(1:length(curSig)) = sig(1:length(curSig)) + curSig;
    end
    
        % normalize to -1:1
%     sig = sig+min(sig);
%     sig = 2*sig/max(sig)-1;
    
    sig = sig/length(channels);

function sig = getFiltSignal(handles,chan)
    % returns filtered signal for playing for given channel       

    % get data
    sig = handles.curSignals(chan,:);
    
    % select second
    xlm = getCurSec(handles)+[-1 0];
    t = linspace(0,(length(sig)-1)/handles.samplingFreq,length(sig));
    sig = sig(t>=xlm(1) & t<xlm(2));

    % amplify and filter (thresholding)
    sig = sig * get(handles.gainEdit,'Value');
    thr =  get(handles.thrSlider,'Value');
    sig(sig<thr & sig>-thr) = 0;    

    
%  loadSignals: loads signals for current position to
%  handles.curSignals using current interface
function handles = loadSignals(handles)
    row = getCurSig(handles);
    sigId=handles.signalIds{row};
    
    % loading sign - not debugged yet
%     if(handles.settings.SHOW_LOADING_LABEL)
%        lm = axis(handles.signalAxes);
% %        ptchH = patch([lm(1) lm(2) lm(2) lm(1)],[lm(3) lm(3) lm(4) lm(4)],[1 1 1],'FaceColor',[1 1 1],'FaceAlpha',.5,'EdgeColor','w','Parent',handles.signalAxes);
% %        uistack(ptchH,'top')
%        txtH = text((lm(1)+lm(2))/2, (lm(3)+lm(4))/2, 'Loading...','BackgroundColor','w','HorizontalAlignment','center','VerticalAlignment','middle','Margin',40,'EdgeColor','k','Parent',handles.signalAxes);
%     end
    
    [curSignals handles.curSigInfo]=handles.interface.getSignalsById(sigId);      
        
    [Nr,Nc] = size(curSignals);
        
    if(Nr>handles.settings.PLOT_CHANNELS)
        errMsg=sprintf('Signal has more channels (%d) than maximum channels set in the settings: handles.settings.PLOT_CHANNELS=%d',Nr,handles.settings.PLOT_CHANNELS);        
        h=errordlg(errMsg,'ERROR: too many channels');
        waitfor(h);
    end
    
    curSignals = double(curSignals);    
    
    % replace NaNs with zeros
    if(any(isnan(curSignals(:))))
        curSignals(isnan(curSignals))=0;
    end
    
    % not an integer number of seconds - pad with zeros
    lenSecTmp = Nc/handles.samplingFreq;
    if(lenSecTmp<ceil(lenSecTmp))
        IntNc = ceil(lenSecTmp)*handles.samplingFreq;
        curSignals = padarray(curSignals,[0 IntNc-Nc],0,'post');
        Nc=IntNc;
    end
        
    
    % decimate signals, store to handles
    % downsample (if required)
    handles.curSignals=zeros(Nr,Nc/handles.settings.DECIMATE_FACTOR);
    if(handles.settings.DECIMATE_FACTOR>1)
        for ci=1:Nr
 			handles.curSignals(ci,:)=decimate(curSignals(ci,:), handles.settings.DECIMATE_FACTOR);
        end        
    else
        handles.curSignals = curSignals;
    end
    % compute seconds from decimated signal
    handles.curSigLenSec = ceil(size(handles.curSignals,2)/handles.samplingFreq);
    clear curSignals    
    
    if(handles.settings.NORMALIZE_SIGNAL_PER_CHANNEL)
        for ci=1:getCurChan(handles)
 			handles.curSignals(ci,:)=normalizeSignalByQuantile(handles.curSignals(ci,:),handles.settings.NORMALIZE_SIGNAL_PER_CHANNEL_QUANTILE);   
        end
        % correct normalization to account for PLOT_STEP
        handles.curSignals = handles.curSignals*handles.settings.PLOT_STEP/150;
    end 
    
    % compute decimated signal for overview
    if(get(handles.overviewChck,'Value')==1)
        handles=computeOverviewSignals(handles);
    end
    
%     % clear loading patch
%      if(handles.settings.SHOW_LOADING_LABEL)
% %        if(~isempty(ptchH)&&ishandle(ptchH))
% %            delete(ptchH);
% %        end
%        if(~isempty(txtH)&&ishandle(txtH))
%            delete(txtH);
%        end       
%     end
    
    
    % store data to gui
    guidata(handles.sigInspectMainWindow, handles)
    nrgChck(nargout,'loadSignals');

% decimate signals for overview window    
function handles = computeOverviewSignals(handles)    

    if(isempty(handles.curSignals))
        return
    end
%     if(~isempty(handles.curSignalsOverview))
%         return
%     end
    
    handles.curSignalsOverview=[];
    
    if(handles.settings.OVERVIEW_DECIMATE_FACTOR>1)
        handles.curSignalsOverview = zeros(size(handles.curSignals,1),ceil(size(handles.curSignals,2)/handles.settings.OVERVIEW_DECIMATE_FACTOR));
        for ii=1:size(handles.curSignals,1)
            handles.curSignalsOverview(ii,:) = decimate(handles.curSignals(ii,:), handles.settings.OVERVIEW_DECIMATE_FACTOR);
        end
    else
        handles.curSignalsOverview = handles.curSignals;
    end
    nrgChck(nargout,'computeOverviewSignals')

    
function curSec = getCurSec(handles)
    curSec = get(handles.secondSelect,'Value');
    
function curSig = getCurSig(handles)
    % gets current signal number from the selectbox
    curSig = get(handles.signalSelect,'Value');
%     [st ed] = regexp(str,'\d*'); 
%     curSig = num2str(str(st:ed));
 
% returns number of channels in currently loaded signal
function curChan = getCurChan(handles)
if(isempty(handles.curSignals))
    warning('empty curSignals. getCurChan probably loaded before setSignal')
    curChan=[];
else
    curChan=size(handles.curSignals,1);
end


function handles = saveAll(handles)
    annotation = handles.annotation;
    filedatetime = datestr(now, 'yyyy-mm-dd-HHMMss');
    
    % save also metadata identification
    interfaceClass = class(handles.interface);
    signalIds= handles.signalIds;
    artifactTypes= handles.settings.ARTIFACT_TYPES;
    
    
    filepath = strrep(handles.settings.ANNOT_DEFAULT_FILENAME,'##', filedatetime);
    save(filepath,'annotation','signalIds','interfaceClass','artifactTypes') ;
    msg = ['annotation saved to ' filepath];
    h = msgbox(msg,'modal');
    disp(msg)
    
    % no changes compared to saved version
    handles.anyChanges = false;
%     guidata(handles.sigInspectMainWindow,handles);
    nrgChck(nargout,'saveAll');


function handles = annotateCurrent(handles, annotNum)
%     ch = str2num(get(get(handles.chSelect,'SelectedObject'),'String'));
    % check if no out-of-bounds annotation is mistakenly requested
    if(annotNum>handles.artifactTypeN)
        %warning(sprintf('requested annotation by type: %d, which is out of bounds (see settings.ARTIFACT_TYPES)',annotNum))
        return
    end

    chs = getSelectedChannels(handles);
    row = getCurSig(handles); %get(handles.signalSelect,'Value');
    sec = getCurSec(handles);
    
    if(isempty(chs))
        handles=playBeep(handles);
        disp('No channels selected => nothing to annotate')
        return
    end
        
    annot = handles.annotation{row};
    if(isempty(annot))
        error('uninitialized annotation. annotateCurrent probably called before setSignal')
    end
    
    
    cur = mean(annot(chs,sec,annotNum));
    if(cur>0.5)
        an = 0;
    else
        an = 1;
    end
    
    annot(chs,sec,annotNum) = an;
    handles.annotation{row} = annot;
    handles.anyChanges = true;
    
    dispAnnotation(handles)
    
    if(~strcmp(handles.overviewPlotAnnotMode,'none'))
        redrawOverview(handles);
    end

    nrgChck(nargout,'annotateCurrent');

 
function handles = markCurrentUnsure(handles, annot)
    warning('unsure deprecated, will be removed soon')
% %     ch = str2num(get(get(handles.chSelect,'SelectedObject'),'String'));
%     chs = getSelectedChannels(handles);
%     row = getCurSig(handles); %get(handles.signalSelect,'Value');
%     sec = getCurSec(handles);
%     
%     if(nargin<2)
%         % negate current state
%         curSel = mean(handles.annotation(row).unsure(chs,sec));
%         if(curSel <0.5)
%             handles.annotation(row).unsure(chs,sec) = 1;
%         else
%             handles.annotation(row).unsure(chs,sec) = 0;
%         end
%     else
%         handles.annotation(row).unsure(chs,sec) = annot;
%     end
%     
%     handles.anyChanges = true; % some changes done
% %     guidata(handles.sigInspectMainWindow, handles)
%     dispAnnotation(handles);
%     nrgChck(nargout,'markCurrentUnsure');


function dispAnnotation(handles)
    % displays current annotation from annot    
    chs = getSelectedChannels(handles);
    
    row = getCurSig(handles);% get(handles.signalSelect,'Value');
    sec = getCurSec(handles);

    annot = handles.annotation{row};
    if(isempty(annot))
        error('uninitalized annotation. dispAnnot probably called before setSignal')
    end
    
    try
        annot = permute(annot(chs,sec,:),[1 3 2]);
    catch e
        error('Annotation dimension mismatch - annotation could not be loaded')
    end
        
%     if(length(chs)==1)
%         % force artifacts in columns, channels in 
%         annot = annot(:)';
%     end
    fntSizes = [8 10 12]; % none some all
    fntWeights = {'normal','bold','bold'};


    for bi = 1:handles.artifactTypeN;
        fntweight = 'normal';
        fntsize = 8;
        
        prop = mean(annot(:,bi));        
        if(prop==1)
            % all selected channels share the settings
            tp = 3;
        elseif(prop>0)
            % some share
            tp = 2;
        else
            % none
            tp = 1;
        end

        % update info on which channels contain this artifact type:
        % get the current label of the artifact button, ...
        currStr = get(handles.(sprintf('annot%dBtn',bi)),'String');
        % ... remove the " (x,y,z)" postfix from it, if any
        idx = find('('==currStr);
        if ~isempty(idx)
            currStr = currStr(1:(idx(1)-2));
        end
        % ... and construct and append a new channel annotation
        if sum(annot(:,bi))>0
            tmp = annot(:,bi)';
            tmp = tmp.*chs;
            tmp = tmp(tmp>0);
            currStr = [currStr ' (' num2str(tmp,'%d') ')'];
        end
        set(handles.(sprintf('annot%dBtn',bi)),'String',currStr);
        set(handles.(sprintf('annot%dBtn',bi)),'FontWeight',fntWeights{tp})
        set(handles.(sprintf('annot%dBtn',bi)),'FontSize',fntSizes(tp))
    end
    
    

function toggleChannel(handles,channel)
    % change state of given channel selector
        
    if(channel>handles.settings.PLOT_CHANNELS)
        % ignore num keys above currently set number of channels
        return
    end
   
    hndl = handles.(sprintf('ch%dChck',channel));
    val =  get(hndl,'Value');
    
    if(strcmp(get(hndl,'Enable'),'on'))
        set(hndl,'Value',~val);
    else
        % not possible (disabled)
%         playBeep(handles);
    end
    
    if(get(handles.hideChck, 'Value'))
        redraw(handles,0)
    end
    dispAnnotation(handles);
    enableDisableSpectroChck(handles);

    

function selectAllChannels(handles,nval,noSpectRedraw)

    nChans = getCurChan(handles);
  
    if(nargin<3)
        noSpectRedraw=false;
    end
    
    if(nargin<2||isempty(nval))
        nval = 1;
        if(length(getSelectedChannels(handles))==nChans)
          % all selected
          nval = 0;
        end
    end
      
    
  for ii=1:nChans
    hndl = handles.(sprintf('ch%dChck',ii));
%     if(strcmp(get(hndl,'Enable'),'on'))
        % channel enabled
        set(hndl,'Value',nval);
        if(nval)
            set(hndl,'Enable','on');
        end
%     end
  end    
   if(get(handles.hideChck, 'Value'))
        redraw(handles,0)
   end
   dispAnnotation(handles);
   enableDisableSpectroChck(handles,noSpectRedraw);
 
function invertChannelSelection(handles)
% invert channel selection
%   row = getCurSig(handles); %get(handles.signalSelect,'Value');   
  nChans = getCurChan(handles);      
    
  for ii=1:nChans
    hndl = handles.(sprintf('ch%dChck',ii));
%     if(strcmp(get(hndl,'Enable'),'on'))
        % channel enabled
        set(hndl,'Value',~get(hndl,'Value'));
%     end
  end    
   if(get(handles.hideChck, 'Value'))
        redraw(handles,0)
   end
   dispAnnotation(handles);
   enableDisableSpectroChck(handles);
  
function sndOn = getSoundState(handles)
    sndOn = get(handles.soundOnChck,'Value');
      
function channels = getSelectedChannels(handles)
% channels = getSelectedChannels(handles) - returns ids of selected channels
fixChannelSelection(handles)
channels = [];
for ii=1:handles.settings.PLOT_CHANNELS
    if(get(handles.(sprintf('ch%dChck',ii)),'Value'))
        % channel ii on
        channels = [channels ii];
    end
end

function fixChannelSelection(handles)

nChans = getCurChan(handles);

for ii=nChans+1:handles.settings.PLOT_CHANNELS
    set(handles.(sprintf('ch%dChck',ii)),'Value',0)
    set(handles.(sprintf('ch%dChck',ii)),'Enable','off')
end

function enableAllChannels(handles)
        % enable all available channel selectors
    for ci = 1:handles.settings.PLOT_CHANNELS
        onoff = 'off';
        val = 0;
        if(ci<=elCount)
            onoff = 'on';
            val=1;
        end
        set(handles.(sprintf('ch%dChck',ci)),'Value',val)
        set(handles.(sprintf('ch%dChck',ci)),'Enable',onoff)
    end                        

% --- Executes on selection change in signalSelect.
function signalSelect_Callback(hObject, eventdata, handles)
% hObject    handle to signalSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% loadDbsSignals(handles);
% selectAllChannels(handles,1);
% uicontrol(handles.sigInspectMainWindow);
% redraw(handles)
handles=setSignal(handles,get(handles.signalSelect,'Value'),1);
guidata(handles.sigInspectMainWindow,handles);


% --- Executes during object creation, after setting all properties.
function signalSelect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to signalSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in nextBtn.
function nextBtn_Callback(hObject, eventdata, handles)
% hObject    handle to nextBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = nextSignal(handles);
guidata(handles.sigInspectMainWindow,handles);

function handles = nextSignal(handles)
    current = getCurSig(handles); %get(handles.signalSelect,'Value');
    if(current<handles.N)
%         str = get(handles.signalSelect,'String');
%         next = str{current+1};
        next = current+1;
        handles = setSignal(handles,next,1);        
    else
        handles=playBeep(handles);
    end
    
    nrgChck(nargout,'nextSignal');


% setSignal - main signal setting function
% sets all necessary internal variables after current
% signal has changed, loads dbsData and redraws
function handles = setSignal(handles,sigId,sec)
    % handles = setSignal(handles,sigId,sec)
    
%     prevHndl = handles;    
%     prevSigId = getCurSig(handles);
%     
    set(handles.signalSelect,'Value',sigId);
%     try
        handles = loadSignals(handles);    
%     catch
%         errordlg(['Error loading signal: ' handles.signalIds{sigId}],'Data load error')
%         set(handles.signalSelect,'Value',prevSigId);
%         return;
%     end
    
    % check whether annot is initialized, if not, do it
    if(isempty(handles.annotation{sigId}))
        handles.annotation{sigId} = emptyCurAnnot(handles);
    end
    
    % last second (-1) - recompute according to current signal length
    if(sec==-1)
        sec = handles.curSigLenSec;
    end
    handles = redrawSecondSelect(handles);
    set(handles.secondSelect,'Value',sec);
    
    % set current signal as seen
    handles.seenSignals(sigId)=1;
    
    selectAllChannels(handles,1,1);
    cla(handles.spectroAxes);
    redraw(handles);
%     playSound(handles)

    nrgChck(nargout,'setSignal');    
    

% --- Executes on button press in prevBtn.
function prevBtn_Callback(hObject, eventdata, handles)
% hObject    handle to prevBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prevSignal(handles)

function handles=prevSignal(handles)
current = get(handles.signalSelect,'Value');
if(current>1)
    prev = current-1;
    handles=setSignal(handles,prev,-1);
else
    % beep
    handles=playBeep(handles);
end
nrgChck(nargout,'prevSignal');

% --- Executes on button press in replayBtn.
function replayBtn_Callback(hObject, eventdata, handles)
% hObject    handle to replayBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = playSound(handles);
guidata(handles.sigInspectMainWindow,handles);

% --- Executes on button press in stopBtn.
function stopBtn_Callback(hObject, eventdata, handles)
% hObject    handle to stopBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
stopSound(handles);

% --- Executes on button press in annot5Btn.
function annot1Btn_Callback(hObject, eventdata, handles)
% hObject    handle to annot5Btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = annotateCurrent(handles,1);
guidata(handles.sigInspectMainWindow,handles);


function annot2Btn_Callback(hObject, eventdata, handles)
% hObject    handle to annot5Btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = annotateCurrent(handles,2);
guidata(handles.sigInspectMainWindow,handles);


% --- Executes on button press in annot5Btn.
function annot3Btn_Callback(hObject, eventdata, handles)
% hObject    handle to annot5Btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = annotateCurrent(handles,3);
guidata(handles.sigInspectMainWindow,handles);



% --- Executes on button press in annot4Btn.
function annot4Btn_Callback(hObject, eventdata, handles)
% hObject    handle to annot4Btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = annotateCurrent(handles,4);
guidata(handles.sigInspectMainWindow,handles);


% --- Executes on button press in annot5Btn.
function annot5Btn_Callback(hObject, eventdata, handles)
% hObject    handle to annot5Btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = annotateCurrent(handles,5);
guidata(handles.sigInspectMainWindow,handles);


% --- Executes on button press in annot6Btn.
function annot6Btn_Callback(hObject, eventdata, handles)
% hObject    handle to annot6Btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = annotateCurrent(handles,6);
guidata(handles.sigInspectMainWindow,handles);


% --- Executes on slider movement.
function thrSlider_Callback(hObject, eventdata, handles)
% hObject    handle to thrSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

set(handles.thrEdit,'String', get(handles.thrSlider,'Value'));
set(handles.thrEdit,'Value', get(handles.thrSlider,'Value'));
redraw(handles,0); % 0=do not adapt gain automatically


function changeThreshold(handles,multip)
% set threshold to a multiple of current value
thr = multip*get(handles.thrSlider,'Value');

% check bounds
thr = min(thr,get(handles.thrSlider,'Max'));
thr = max(thr,get(handles.thrSlider,'Min'));

set(handles.thrSlider,'Value',thr)

set(handles.thrEdit,'String', thr);
set(handles.thrEdit,'Value', thr);
redraw(handles,0); % 0=do not adapt gain automatically

function changeGain(handles,multip)
% set threshold to a multiple of current value
gain = multip*get(handles.gainSlider,'Value');

% check bounds
gain = min(gain,get(handles.gainSlider,'Max'));
gain = max(gain,get(handles.gainSlider,'Min'));

set(handles.gainSlider,'Value',gain)

set(handles.gainEdit,'String', gain);
set(handles.gainEdit,'Value', gain);
redraw(handles,0); % 0=do not adapt gain automatically


% --- Executes during object creation, after setting all properties.
function thrSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thrSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function gainSlider_Callback(hObject, eventdata, handles)
% hObject    handle to gainSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.gainEdit,'String', get(handles.gainSlider,'Value'));
set(handles.gainEdit,'Value', get(handles.gainSlider,'Value'));
redraw(handles,0); % 0=do not adapt gain automatically
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function gainSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gainSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function gainEdit_Callback(hObject, eventdata, handles)
% hObject    handle to gainEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gainEdit as text
%        str2double(get(hObject,'String')) returns contents of gainEdit as a double


% --- Executes during object creation, after setting all properties.
function gainEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gainEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function thrEdit_Callback(hObject, eventdata, handles)
% hObject    handle to thrEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thrEdit as text
%        str2double(get(hObject,'String')) returns contents of thrEdit as a double


% --- Executes during object creation, after setting all properties.
function thrEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thrEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in secondSelect.
function secondSelect_Callback(hObject, eventdata, handles)
% hObject    handle to secondSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns secondSelect contents as cell array
%        contents{get(hObject,'Value')} returns selected item from secondSelect
showSecond(handles);
% uicontrol(handles.sigInspectMainWindow);



% --- Executes during object creation, after setting all properties.
function secondSelect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to secondSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in nextSecBtn.
function nextSecBtn_Callback(hObject, eventdata, handles)
% hObject    handle to nextSecBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = nextSecond(handles);
guidata(handles.sigInspectMainWindow,handles);

function handles = nextSecond(handles)
curSec = getCurSec(handles);
if(curSec<handles.curSigLenSec)
    set(handles.secondSelect,'Value',curSec+1);
    showSecond(handles)
else
    % beep and go to previous signal
    handles=playBeep(handles);
    if(getCurSig(handles)<handles.N)
        % if next sig. available - move to it
        handles = nextSignal(handles);
        set(handles.secondSelect,'Value',1);
    end    
end
handles = playSound(handles);

nrgChck(nargout,'nextSecond');


% --- Executes on button press in prevSecBtn.
function prevSecBtn_Callback(hObject, eventdata, handles)
% hObject    handle to prevSecBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = prevSecond(handles);
guidata(handles.sigInspectMainWindow,handles)


function handles = prevSecond(handles)
    curSec = getCurSec(handles);
    if(curSec>1)
        set(handles.secondSelect,'Value',curSec-1);
    else
        % beep and go to previous signal
        handles=playBeep(handles);
        if(getCurSig(handles)>1)
            handles = prevSignal(handles);
            set(handles.secondSelect,'Value',handles.curSigLenSec);
        end
    end
    showSecond(handles) 
    handles = playSound(handles);
    
    nrgChck(nargout,'prevSecond');
    

        


% --- Executes on button press in SaveBtn.
function SaveBtn_Callback(hObject, eventdata, handles)
% hObject    handle to SaveBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = saveAll(handles);
guidata(handles.sigInspectMainWindow,handles);


function keyPressHandler(evt,handles)        
% handles=guidata(hObject);
switch(evt.Key)
    case 'rightarrow'
        handles=nextSecond(handles);
    case 'leftarrow'
        handles=prevSecond(handles);
   
    case 'pagedown'
        handles=nextSecond(handles);
    case 'pageup'
        handles=prevSecond(handles);
     
    case 'numpad1'
       toggleChannel(handles,1)
    case 'numpad2'
        toggleChannel(handles,2)
    case 'numpad3'
        toggleChannel(handles,3)
    case 'numpad4'
        toggleChannel(handles,4)
    case 'numpad5'
        toggleChannel(handles,5)
    case 'numpad6'
        toggleChannel(handles,6)
    case 'numpad7'
        toggleChannel(handles,7)
    case 'numpad8'
        toggleChannel(handles,8)
    case 'numpad9'
        toggleChannel(handles,9)
    case 'numpad0'
        toggleChannel(handles,10)
    case '1'
       toggleChannel(handles,1)
    case '2'
        toggleChannel(handles,2)
    case '3'
        toggleChannel(handles,3)
    case '4'
        toggleChannel(handles,4)
    case '5'
        toggleChannel(handles,5)
    case '6'
        toggleChannel(handles,6)
    case '7'
        toggleChannel(handles,7)
    case '8'
        toggleChannel(handles,8)
    case '9'
        toggleChannel(handles,9)
    case '0'
        toggleChannel(handles,10)
    case 'a'
        selectAllChannels(handles)
    case 'i'
        invertChannelSelection(handles)   
    
    case 'u'
        handles=markCurrentUnsure(handles);
    
    case 'f1'
        handles=annotateCurrent(handles,1);
    case 'f2'
        handles=annotateCurrent(handles,2);
    case 'f3'
        handles=annotateCurrent(handles,3);
    case 'f4'
        handles=annotateCurrent(handles,4);
    case 'f5'
        handles=annotateCurrent(handles,5);
    case 'f6'
        handles=annotateCurrent(handles,6);
%     case 'f7'
%         handles=annotateCurrent(handles,7);
%     case 'f8'
%         handles=annotateCurrent(handles,7);
%     case 'f9'
%         handles=annotateCurrent(handles,7);
    
    case 'r'
        set(handles.soundOnChck,'Value',true);
        handles=playSound(handles);
    case 's'
        handles=saveAll(handles);
    case 'o'
        handles=toggleSound(handles);
    case 'v'
        handles=toggleOverview(handles);
    case 'g'
        redraw(handles,2)
    case 'f'
        toggleSpectrogram(handles);
    case 'w'
        toggleSpectrogramWhole(handles)
    case 'h'
        toggleHideUnselected(handles)

        
    case 'add'
        changeThreshold(handles,1.05)
    case 'subtract'
        changeThreshold(handles,0.95)
    case 'uparrow'
        changeGain(handles,1.05)
    case 'downarrow'
        changeGain(handles,0.95)
    
        
        
    otherwise
%         disp(['key pressed: ' evt.Key]) % DEBUG:display code of key
%         pressed

end

guidata(handles.sigInspectMainWindow,handles);



    
% --- Executes when user attempts to close sigInspectMainWindow.
function sigInspectMainWindow_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to sigInspectMainWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
% any changes
if(exist('handles','var') && isstruct(handles) && isfield(handles,'anyChanges')&& handles.anyChanges)

    
    % remove any finished players from the internalPlayer queue such that matlab
    % frees system resources associated with the players and can exit cleanly
    internalPlayer(handles,handles.internal.INTERNAL_PLAYER_REFRESH);
    
    r = questdlg('There are unsaved changes. Save before closing?','Save annotation?','Save','Discard','Cancel','Save');

    switch(r)
        case 'Save'
            handles=saveAll(handles);
        case 'Cancel'
            return;
    end
end

% close overview + main fig.
if(isstruct(handles)&&isfield(handles,'overviewFig')&&~isempty(handles.overviewFig) && ishandle(handles.overviewFig))
    delete(handles.overviewFig) % close overview window
end
delete(hObject)


% --- universal channel checkbox callback
function chanChckCallback(hObject, callbackdata)
handles=guidata(hObject);
dispAnnotation(handles)
enableDisableSpectroChck(handles);


% --- Executes on key press with focus on sigInspectMainWindow or any of its controls.
function sigInspectMainWindow_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to sigInspectMainWindow (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
keyPressHandler(eventdata, handles)


% --- Executes on button press in soundOnChck.
function soundOnChck_Callback(hObject, eventdata, handles)
% hObject    handle to soundOnChck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of soundOnChck
snd = getSoundState(handles);
if(snd)
    handles = playSound(handles);
    guidata(handles.sigInspectMainWindow,handles);
else
    stopSound(handles);
end
    


% --------------------------------------------------------------------
% function SaveMenuItem_Callback(hObject, eventdata, handles)
% % hObject    handle to SaveMenuItem (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% handles = saveAll(handles);
% guidata(handles.sigInspectMainWindow,handles);


% --- Executes on button press in spectroChck.
function spectroChck_Callback(hObject, eventdata, handles)
% hObject    handle to spectroChck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of spectroChck
redrawSpectrogram(handles,true); % force redraw
showSecond(handles);


% --------------------------------------------------------------------
function saveTool_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to saveTool (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = saveAll(handles);
guidata(handles.sigInspectMainWindow,handles);


% --------------------------------------------------------------------
function openTool_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to openTool (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file path] = uigetfile('*.mat','open dbsArtifactAnnotation*.mat file (saved annotation)');
% if ~isequal(file, 0)
%     open(file);
% end
filePath = [path file];

if(~isempty(filePath) && any(filePath))
    data = load(filePath);
    
    % annotation variable must be present
    if(~isfield(data,'annotation'))
        errordlg('No annotation found in file. Annotation not loaded','Annotation load error');
        return;
    end
    
    % signalIds variable must be present & match loaded signalIds
    if(~isfield(data,'signalIds'))
        errordlg('No signalIds variable found in the file. Annotation could not be loaded','Annotation load error');
        return
    else
        data.signalIds=data.signalIds(:); % convert to col vector
        idsOk = length(data.signalIds)==length(handles.signalIds);
        idsOk = idsOk && iscell(data.signalIds) && all(~cellfun(@isempty,data.signalIds)) && all(strcmp(data.signalIds,handles.signalIds));
        if(~idsOk)
            errordlg('SignalIds in annotation file do not match loaded ones. Annotation not loaded','Annotation load error');
            return
        end                               
    end
    
    % signalIds and annotation must be of the same length
    if(length(data.signalIds)~=length(data.annotation))
        errordlg('Annotation and signalIds in the file are of different length. Annotation not loaded','Annotation load error');
        return;
    end 

    % metadata check: interfaceClass, artifactTypes - may be joined with current on request if these are not present or are empty
    if( (~isfield(data,'interfaceClass') || isempty(data.interfaceClass))||...
        (~isfield(data,'artifactTypes')  || isempty(data.artifactTypes)))
        
    % no metadata present in file
        res = questdlg('Annotation does not contain metadata (artifactTypes or interfaceClass). Do you want to convert the loaded annotation to conform to current data structure?','Loading annotation');
        switch(res)
            case 'Yes'
                % join with current, extend annotation
                convAnnot = convertAnnot(data.annotation,handles.settings.ARTIFACT_TYPES);
                
                % check if conversion successful
                if(isempty(convAnnot))
                    errordlg('Annotation could not be converted','Annotation load error');
                    return;
                end
                setAnnot(handles,convAnnot)
                msgbox('Conversion ok. Annotation loaded.','Annotation loaded','modal')
                return
            otherwise
                % cancel loading
                msgbox('Annotation not loaded','Annotation load error','warn','modal')
                return
        end
    end 
    
    % check metadata if present
    % check interface type
    if((handles.settings.ANNOT_FILE_CHECK_INTERFACE_CLASS && ~strcmp(data.interfaceClass,class(handles.interface))))
        res = questdlg(sprintf('Interface type in annotation file (%s) does not match current data (%s). Do you want to ignore the mismatch?',data.interfaceClass,class(handles.interface)),'Loading annotation');
        switch(res)
            case 'Yes'
            otherwise
                msgbox('User abort. Annotation NOT loaded','User abort','warn','modal')
                return
        end
    end
    % check artifact types
    if(handles.settings.ANNOT_FILE_CHECK_ARTIFACT_TYPES  && ~isequal(data.artifactTypes,handles.settings.ARTIFACT_TYPES))        
        res = questdlg(sprintf('Annotation file artifact types differ from currently loaded ones ({%s} vs {%s}). Do you want to attempt conversion to current data? This will join all artifact types into one',sprintf('%s,',data.artifactTypes{:}),sprintf('%s,',handles.settings.ARTIFACT_TYPES{:})),'Loading annotation');
%         res = questdlg('Annotation file artifact types differ from currently loaded ones. Do you want to attempt conversion to current data? This will join all artifact types into one','Loading annotation');
        switch(res)
            case 'Yes'
                % join with current, extend annotation
                convAnnot = convertAnnot(data.annotation,handles.settings.ARTIFACT_TYPES);
                
                % check if conversion successful
                if(isempty(convAnnot))
                    errordlg('Annotation could not be converted','Annotation load error');
                    return;
                end
                setAnnot(handles,convAnnot)
                msgbox('Conversion ok. Annotation loaded.','Annotation loaded','modal')
                return
            otherwise
                % cancel loading
                msgbox('User abort. Annotation not loaded','User abort','warn','modal')
                return
        end
    end  
    setAnnot(handles, data.annotation);
    msgbox('Annotation loaded','Annotation loaded','modal')
    
end

function setAnnot(handles,newAnnot)
    % sets new annotation to handles, does no error checking
    handles.annotation = newAnnot;
    guidata(handles.sigInspectMainWindow, handles);
    dispAnnotation(handles);    
    redrawSignalSelect(handles);
    redrawOverview(handles,0); % 0=not just patch


% convert annotation structure to current number of artifacts etc
function annot2 = convertAnnot(annot,artifactTypes)
    annot2=cell(size(annot));
    if(isempty(annot))
        return
    end

    fprintf('-- annotation conversion + check started\n')
    % check annot one by one
    for ii=1:length(annot)
        an = annot{ii};

        if(isempty(an))
            warning('   empty annotation at row #%d - SKIPPING',ii);
            annot2{ii}=[];
            continue;
        end

        if(size(an,3)>length(artifactTypes))
            an = any(an,3);
            fprintf('  row #%d - more artifact types than current - joining artifact types into one! \n',ii);
        end

        % zero-pad annotation to current artifactTypes count
        dimsToPad= length(artifactTypes)-size(an,3);
        an = padarray(an,[0 0 dimsToPad],0,'post');
        annot2{ii}=an;
    end
    fprintf(' > conversion ok')
    

function normSignal=normalizeSignalByQuantile(signal,qntl)
signal = signal-mean(signal);
q=max(abs(quantile(signal,[qntl,1-qntl])));
if iqr(signal)>0
    normSignal=signal/q*10;
else
    % do not normalize all-zeros signal
    normSignal=signal;
end


function [ind] = constSamples(d, minLen,eps)
% ind = constSamples(d, minLen,eps) - identify constant regions of length
% at least minLen samples in the sequence d
% OUT: ind - logical indices of constant regions in sequence d
% E. Bakstein 10.1.2013

ind = false(size(d));

if(nargin<3 || isempty(eps))
    eps = 0;
end


constLen = 1;
for i = 2:length(d)  
%     if(d(i-1)*(1-eps) <= d(i) && d(i) <= d(i-1)*(1+eps))
    if(abs(d(i-1)-d(i))<=eps)
        constLen = constLen + 1;
    else
        if(constLen>=minLen)
            ind(i-constLen:i-1) = true;
        end
        constLen = 1;
    end
end

% handle constant regions at the end:
if(constLen>=minLen)
    ind(i-constLen+1:end) = true;
end


% --- Executes on button press in overviewChck.
function overviewChck_Callback(hObject, eventdata, handles)
% hObject    handle to overviewChck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
isOn=get(hObject,'Value')==1;
toggleOverview(handles,isOn);


% --- Executes on button press in spectroWholeChck.
function spectroWholeChck_Callback(hObject, eventdata, handles)
% hObject    handle to spectroWholeChck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
toggleSpectrogramWhole(handles,get(hObject,'Value'));


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over sigInfoTxt.
function sigInfoTxt_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to sigInfoTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clipboard('copy', get(handles.sigInfoTxt,'str'));
 msgbox('SignalId copied to clipboard','SignalId copied to clipboard','modal')
               


% --- Executes on button press in hideChck.
function hideChck_Callback(hObject, eventdata, handles)
% hObject    handle to hideChck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of hideChck
toggleHideUnselected(handles, get(hObject, 'Value'))


% --- Executes on mouse press over signals axes background.
function signalAxes_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to signalAxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.internal.lastButtonPressedInSignalAxes==1 && eventdata.Button==3
    % display the time interval of the last and the current mouse position in the signals axes
    % (in seconds, aligned to 10ms such that the interval reported contains the real interval)
    fprintf(1,' %6.2f, %6.2f, \n',...
        floor(handles.internal.lastPositionInSignalAxes*100)/100,...
        ceil(eventdata.IntersectionPoint(1)*100)/100);
end
handles.internal.lastPositionInSignalAxes=eventdata.IntersectionPoint(1);
handles.internal.lastButtonPressedInSignalAxes=eventdata.Button;
guidata(handles.sigInspectMainWindow, handles);

function postZoom(obj,e)
xl = get(e.Axes,'xlim');
if(diff(xl)==1) % restored to 1s x-axis
    showSecond(guidata(obj))
end

