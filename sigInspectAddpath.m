% sigInspectAddpath - add sigInspect directories to path
% loaded automatically during sigInspect initialization
% E. Bakstein 20161209

% subdirectories to be added
dirs = {'interfaces'};

% reference - path to sigInspectInit.m file
[pth,~,~] = fileparts(which(mfilename));
for ii=1:length(dirs)
    addpath(genpath([pth '/' dirs{ii}]));
end
