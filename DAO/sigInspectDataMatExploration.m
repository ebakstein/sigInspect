% implements data loading for DAO exploration patients stored in Mat file
% by explDataAndAnnotToMat
% 
%  sigInspectDataDaoExploration(filename.mat) - 
%   IN: filename.mat > matfile with signals and signalIds
% 
% E. Bakstein 2015-11-26
% 
classdef sigInspectDataMatExploration < sigInspectDataInterface
    properties
        signalIds={};  % exploration ids such as 192t
%         settings={}; % default from abstract class
%         dbs;
    end
        methods
        function obj=sigInspectDataMatExploration(filepath)
            if(nargin<1)
                error('provide filename of *.mat file with signals and signalIds');
            end
            
            if(~exist(filepath,'file'))
                error('file >%s< not found',filepath)
            end
            L = load(filepath);
            if(~isfield(L,'data'))
                error('no "data" variable in matfile')
            end
            if(~isfield(L,'signalIds'))
                error('no "signalIds" variable in matfile')
            end
            
            % constructor with list of signal ids
            obj.signalIds=signalIds; 
            
            global dbs;
%             obj.dbs=dbs;
            obj.settings.SAMPLING_FREQ=dbs.samplingFreq;
%            obj.settings.ARTIFACT_TYPES=dbs.artifactTypes(2:end,3)';
            obj.settings.ARTIFACT_TYPES={'POW'  'BASE'  'FREQ'  'IRIT'  'OTHR'  'ARTIF'};
            obj.settings.REVERSE_CHANNEL_ORDER = 1;
            
        end
        
        % return list of signal ids
        function signalIds = getSignalIds(obj)
            if(isempty(obj.signalIds))
                error('signal Ids not initialized yet')
            end
            signalIds = obj.signalIds;            
        end
        
        function [signals chInfo]= getSignalsById(obj,signalId)
            % loads raw signal data by exploration Id
            data=daoGet(dbsExplSignalIDParse(signalId),rawSignalDataType);
            
            % add all channels
            signals=[]; % data matrix: signals in rows
            chInfo='';  % channel labels
            for si=1:length(data)
                signals(si,:)=data(si).data;
                chInfo = sprintf('%s %d:%s:%s',chInfo,si,data(si).electrode(1:3),data(si).area);
            end
            % DEBUG
            signals=signals/10;
            if(nargout<2)
                clear chInfo;
            end
        end
    end
end