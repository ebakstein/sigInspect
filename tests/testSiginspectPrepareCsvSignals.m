% testSiginspectPrepareCsvSignals - prepare CSV signals for sigInspect demo
% E. Bakstein 2015-07-07
global dbs;
daoIds={'18t245p3x','20t220p21x'};
newnames={'signal1','signal2'};
pth='artifacts/sigInspect/csvDemo/%s.csv';
for ii=1:length(daoIds)
    sigId = daoIds{ii};
    d = daoGet(dbsExplSignalIDParse(sigId),rawSignalDataType);
    signals=[];
    for si=1:3 %length(d)-1
        sig = decimate(d(si).data(1:dbs.samplingFreq*4),4)/100;
        signals(:,si)=sig';
    end        
    csvwrite(sprintf(pth,newnames{ii}),signals)
end