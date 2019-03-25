% Implementation of the 'sigInspectDataInterface' operating with
% resting state data via DAO.
% 
% Author: T.Sieger 2015-07-05
% 
classdef sigInspectDataResting < sigInspectDataInterface
    properties
        signalIds={};  % resting state signal IDs such as '30rp1'
        % settings={}; % defined in 'sigInspectDataInterface' already
    end

    methods
        function this=sigInspectDataResting()
            global dbs;

            this.signalIds={};

            for iPat=1:length(dbs.restingPatients)
                patient=dbs.restingPatients{iPat};
                for iPos=1:length(patient.positions)
                    id=sprintf('%drp%d',patient.id,iPos);
                    this.signalIds=[this.signalIds {id}];
                end
            end
            this.signalIds=this.signalIds;

            this.settings.DECIMATE_FACTOR=8;
            this.settings.SPECTROGRAM_FREQ_LIMS=[0 dbs.samplingFreq/2/this.settings.DECIMATE_FACTOR];
            this.settings.SPECTROGRAM_NFFT=128;
            this.settings.OVERVIEW_DECIMATE_FACTOR=16;
            this.settings.OVERVIEW_GAIN=20;
            %this.settings.ARTIFACT_TYPES={'AUTO'};
            %this.settings.ANNOT_FILE_CHECK_ARTIFACT_TYPES=0;
            this.settings.ENABLE_WHOLE_SPECTROGRAM=1;
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
                signals(si,:)=data(si).data;%(1:240000);%
                chInfo = sprintf('%s %d:%s:%s',chInfo,si,data(si).electrode(1:3),data(si).area);
            end
            if(nargout<2)
                clear chInfo;
            end
        end
    end
end
