% testSigInspect
% testing sigInspect initialization with different properties
% E. Bakstein 2015-06-04

d=daoGet(dbsExplSignalIDParse('15t218p21l'),rawSignalDataType);
sig1=d.data;

d2=daoGet(dbsExplSignalIDParse('15t218p22x'),rawSignalDataType);
sig2=[];for ii=1:length(d2);sig2=[sig2;d2(ii).data];end
sig2 = [sig2 sig2]; % 20s - twice the same signal

d3=daoGet(dbsExplSignalIDParse('15t218p23x'),rawSignalDataType);
sig3=[];for ii=1:3;sig3=[sig3;d3(ii).data];end

sig4=[sig1;sig3;sig3;sig3];


% loading single signal (implicit call to sigInspectDataBasic)
sigInspect({sig1,sig2,sig3,sig4})

%% different sampling rate - basic interface
d=daoGet(dbsExplSignalIDParse('15t218p21l'),rawSignalDataType);
sig1=d.data;

d2=daoGet(dbsExplSignalIDParse('15t218p22x'),rawSignalDataType);
sig2=d2.data;

sigInspect({sig1,sig2},12000)

%% shorter signals - automatically pad NaNs
d=daoGet(dbsExplSignalIDParse('15t218p21l'),rawSignalDataType);
sig1=d.data;
sig1=sig1(1:end-3571);

sigInspect({sig1})


%% loading sith dbs DAO
intf=sigInspectDataDaoExploration({  '15t218p21x','20t220p21x','18t245p3x','61t233p22m'});
intf.settings.ARTIFACT_TYPES={'TYPE1','TYPE2','UNKNOWN'};
intf.settings.REVERSE_CHANNEL_ORDER=0;
sigInspect(intf);