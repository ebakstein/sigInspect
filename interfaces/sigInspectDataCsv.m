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
        function obj=sigInspectDataCsv(dirPath)
            obj.dirPath = dirPath;
            obj.settings.SAMPLING_FREQ=6000; % 6kHz sampling rate
            obj.settings.PLOT_STEP=1.5;      % distance between channels on the y-axis
            obj.settings.ARTIFACT_TYPES={'Type A','Type B','Unsure'};            
        end        
        % return list of signal ids - load all csv files from a directory,
        % use filenames as signalIds
        function signalIds = getSignalIds(obj)
            if(isempty(obj.dirPath) || ~exist(obj.dirPath,'dir'))
                error('Directory %s empty or does not exist',obj.dirPath)
            end
            lst=dir([obj.dirPath '/*.csv']);
            signalIds = {lst(:).name}';
        end
        
        % read signals based on signalId (=filename)
        function [signals chInfo]= getSignalsById(obj,signalId)
            filePath = [obj.dirPath '/' signalId];
            if(~exist(filePath,'file'))
                error('File %s does not exist',filePath)
            end
            chInfo='';
            signals=csvread(filePath)'; % load + transpose to have signals in rows
        end
    end
end