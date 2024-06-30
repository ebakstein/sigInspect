% Example sigInspectDataInterface for loading data from CSV files in a
% directory
% E. Bakstein 2015-07-07
% 
classdef sigInspectDataCsv < sigInspectDataInterface 
    % define class, inherit sigInspectDataInterface abstract class
    properties
        dirPath='';
    end
    methods
        % constructor - just store the path, set settings
        function obj=sigInspectDataCsv(dirPath, varargin)
            obj.dirPath = dirPath;
            %obj.settings.PLOT_STEP=1.5;      % distance between channels on the y-axis
            obj.settings.ARTIFACT_TYPES={'Type A','Type B','Unsure'};

            samplingFreq = [];
            % Check if the second argument is provided and it is numeric
            if ~isempty(varargin) && isnumeric(varargin{1})
                samplingFreq = varargin{1};
            end

            % If sampling frequency is still empty, prompt the user
            if isempty(samplingFreq)
                % If the second argument is not provided or not numeric, prompt the user
                prompt = {'Enter sampling frequency:'};
                dlgtitle = 'Input';
                dims = [1 35];
                definput = {'12000'}; % Default value
                samplingFreqInput = inputdlg(prompt, dlgtitle, dims, definput);
                if isempty(samplingFreqInput)
                    disp('No sampling frequency provided. Using default 12000 Hz.');
                    samplingFreq = 12000; % Default sampling frequency
                else
                    samplingFreq = str2double(samplingFreqInput{1}); % Convert from cell to double
                end
            end
            % Set the sampling frequency
            obj.settings.SAMPLING_FREQ = samplingFreq;
        end        
        % return list of signal ids - load all csv files from a directory,
        % use filenames as signalIds
        function signalIds = getSignalIds(obj)
            if exist(obj.dirPath, 'file') == 2 % If it's a file
                [~, name, ext] = fileparts(obj.dirPath);
                signalIds = {[name ext]}; % Return the file name as the only signalId
            elseif exist(obj.dirPath, 'dir') == 7 % If it's a directory
                lst = dir([obj.dirPath '/*.csv']);
                signalIds = {lst(:).name}';
            else
                error('Path %s does not exist', obj.dirPath);
            end
        end
        
        % read signals based on signalId (=filename)
        function [signals, chInfo] = getSignalsById(obj, signalId)
            if exist(obj.dirPath, 'file') == 2 % If dirPath is actually a file path
                filePath = obj.dirPath; % Use it directly
            else
                filePath = fullfile(obj.dirPath, signalId); % Construct the file path
            end
            if ~exist(filePath, 'file')
                error('File %s does not exist', filePath);
            end
            chInfo = '';
            signals = readmatrix(filePath)'; % Load and transpose to have signals in rows
        end
    end
end