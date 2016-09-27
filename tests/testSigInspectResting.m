% testSigInspectResting
% Testing sigInspect with long resting-state data.
%
% Author: T. Sieger 2015-06-09

    global dbs;

    dataSel=daoDataSelector();
    dataSel.missingDataPolicy=dataSel.missingDataPolicy_ALLOW_INCOMPLETE_EPISODES;
    dataSel.experimentPart=1;
    patient=dbs.restingPatients{1};
    dataSel.patientId=patient.id;
    positionIdx=1;
    position=patient.positions{positionIdx};
    dataSel.position=position;
    dataSel=dbsFillDataSelSide(dataSel);
    electrodeIdx=1;
    %area=patient.areas{positionIdx}{electrodeIdx};
    %electrode=patient.electrodes{positionIdx}{electrodeIdx};
    %dataSel.electrode=electrode;
    x=daoGet(dataSel,rawSignalDataType());
    x
    s=zeros(length(x),length(x(1).data));
    for i=1:length(x);s(i,:)=x(i).data;end
    clear x
    sigInspect(s)
