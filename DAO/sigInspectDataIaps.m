% Implementation of the 'sigInspectDataInterface' operating with
% IAPS data via DAO.
% 
% Author: T.Sieger 2015-07-09
% 
classdef sigInspectDataIaps < sigInspectDataInterface
    properties
        signalIds={};  % resting state signal IDs such as '30rp1'
        % settings={}; % defined in 'sigInspectDataInterface' already
    end

    methods
        function this=sigInspectDataIaps()
            global dbs;

            this.signalIds={};

            for iPat=1:length(dbs.iapsPatients)
                patient=dbs.iapsPatients{iPat};
                for iPos=1:length(patient.positions)
                    id=sprintf('%d%s',patient.id,patient.positions{iPos});
                    this.signalIds=[this.signalIds {id}];
                end
            end
            %this.signalIds=this.signalIds(1:2);

            %this.settings.DECIMATE_FACTOR=1;
            %this.settings.SPECTROGRAM_FREQ_LIMS=[0 min(3000,dbs.samplingFreq/2/this.settings.DECIMATE_FACTOR)];
            %this.settings.SPECTROGRAM_NFFT=128;
            %this.settings.OVERVIEW_DECIMATE_FACTOR=16;
            %this.settings.ARTIFACT_TYPES={'AUTO'};
            %this.settings.ANNOT_FILE_CHECK_ARTIFACT_TYPES=0;
            this.settings.ENABLE_WHOLE_SPECTROGRAM = 1;
            %this.settings.ARTIFACT_TYPES=dbs.artifactTypes(2:end,3);
            this.settings.REVERSE_CHANNEL_ORDER=0;

        end
        
        % return list of signal ids
        function signalIds = getSignalIds(this)
            signalIds = this.signalIds;
        end
        
        function [signals chInfo]= getSignalsById(this,signalId)
            % loads raw signal data by exploration Id
            data=daoGet(dbsPatientSignalIdParse(signalId),rawSignalDataType);
            
            % add all channels
            signals=[]; % data matrix: signals in rows
            chInfo='';  % channel labels
            for si=1:length(data)
                tmp=reshape(data(si).data',1,prod(size(data(si).data)));
                signals(si,:)=tmp;
                chInfo = sprintf('%s %d:%s:%s',chInfo,si,data(si).electrode(1:3),data(si).area);
            end
            size(signals)
            any(any(isnan(signals)))
            signals(1,1:10)
            if(nargout<2)
                clear chInfo;
            end
        end
    end
end
