% artifAnnot2016
% E. Bakstein 2016-02-02

% addpath(genpath('sigInspect'))

%% LOAD SIGNALS FOR ANNOTATION

dataH = sigInspectDataMatDir('c:/Users/Eda/Projects/DBS/Data/artifacts/artifAnnot2016/EB/Data/','d');
dataH.settings.ARTIFACT_TYPES = {'POW', 'BASE','FREQ','IRIT','OTHR'}; % DBS DAO artif. types
dataH.settings.SAMPLING_FREQ = 24000; 
sigInspect(dataH);

% + SAVE/LOAD ANNOTATION MANUALLY IN SIGINSPECT (OPEN or open icon)


% %% KOSICE
% 
% dataH = sigInspectDataMatDir('Data/Kos_*.mat','d');
% dataH.settings.ARTIFACT_TYPES = {'POW', 'BASE','FREQ','IRIT','OTHR'}; % DBS DAO artif. types
% dataH.settings.SAMPLING_FREQ = 24341; 
% sigInspect(dataH);
