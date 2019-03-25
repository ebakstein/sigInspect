% implements data loading for DAO exploration patients
% 
% E. Bakstein 2015-05-29
% 
classdef sigInspectDataDaoExploration < sigInspectDataInterface
    properties
        signalIds={};  % exploration ids such as 192t
%         settings={}; % default from abstract class
%         dbs;
    end
        methods
        function obj=sigInspectDataDaoExploration(signalIds)
            if(nargin<1)
                error('provide list of signalIds');
            end
            
            % check input format
            if(~iscell(signalIds))
                if(ischar(signalIds))
                    signalIds = {signalIds};
                else
                    error('signalIds must be a cell array of string signalIds')
                end
            end
            
            % checkIds format
            isOK=cellfun(@(x) ~isempty(x) && ischar(x) && ~isempty(regexp(x,'^[0-9]+t[0-9]+p[0-9]+[a-z]$')),signalIds);
            if(any(~isOK))
                probIdsStr=sprintf('%s, ',signalIds{~isOK});
                error('some (%d) signal ids are empty or in invalid format: %s',sum(~isOK),probIdsStr);
            end                       
            
            % constructor with list of signal ids
            obj.signalIds=signalIds; 
            
            global dbs;
%             obj.dbs=dbs;
            obj.settings.SAMPLING_FREQ=dbs.samplingFreq;
            obj.settings.ARTIFACT_TYPES=dbs.artifactTypes(2:end,3)';
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