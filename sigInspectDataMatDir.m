% sigInspectDataInterface for loading data from MAT files in a
% directory
% E. Bakstein 2016-01-25
% 
classdef sigInspectDataMatDir < sigInspectDataInterface 
    % define class, inherit sigInspectDataInterface abstract class
    properties
        dirPath='';
        matTmpl='';
        sigVar='';
    end
    methods
        % constructor - just store the path/mat template, set settings
        function obj=sigInspectDataMatDir(dirPath,sigVar)
            % mat file template or directory path
            matTmpl = dirPath;
            if(isempty(strfind(matTmpl,'\*')))
                matTmpl = [matTmpl '/*.mat'];
            else
                [dirPath,~,~] = fileparts(dirPath);
            end
            obj.matTmpl = matTmpl;
            obj.dirPath = dirPath;
            
            % variable with signal
            if(nargin<2)
                sigVar = 'data';
            end            
            obj.sigVar = sigVar;
        end        
        % return list of signal ids - load all csv files from a directory,
        % use filenames as signalIds
        function signalIds = getSignalIds(obj)
            if(isempty(obj.dirPath) || ~exist(obj.dirPath,'dir'))
                error('Directory %s empty or does not exist',obj.dirPath)
            end            
            lst=dir(obj.matTmpl);
            signalIds = {lst(:).name}';
        end
        
        % read signals based on signalId (=filename)
        function [signals chInfo]= getSignalsById(obj,signalId)
            filePath = [obj.dirPath '/' signalId];
            if(~exist(filePath,'file'))
                error('File %s does not exist',filePath)
            end
            
            signals=load(filePath)'; % load + transpose to have signals in rows
            
            if(~isfield(signals,obj.sigVar))
                error('File %d does not contain signal variable %s',filePath,obj.sigVar);                
            end
            
            signals=signals.(obj.sigVar);
            chInfo='';            
        end
    end
end