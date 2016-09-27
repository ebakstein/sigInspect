function annot = sigInspectClassify(signal,fs,method)
% function annot = sigInspectClassify(signal,varargin)
% classify artifacts in each second of provided micro-EEG signal
% 
% IN: 
%   signal - micro-EEG signal vector
%   fs     - sampling frequency in Hz
%   method - classification method:
%           'psd' - normalized PSD spectrum thresholding, based on [1]
%                  (default) 91% train /88% test set accuracy on EMBC data
%           UNFINISHED: 'tree'     - pre-trained decision tree, based on multiple features - 
% OUT:
%   annot - logical vector of annotation for each second of input signal.
%           true = artifact, false = clean signal
% 
% E. Bakstein 2015-06-29
% 
% [1] Bakstein, E. et. al.: Supervised Segmentation of Microelectrode Recording Artifacts Using Power Spectral Density, in Proceedings of IEEE EMBS, 2015


% input checks
if(nargin<2)
    error('sampling frequency must be specified')
end

if(nargin<3 || isempty(method))
    method='psd';
end
    
if(isempty(signal))
    annot=[];
    return;
end

if(length(signal)<fs)
    error('signal must be at least 1s long');
end

N = length(signal); % number of samples
Ns = ceil(N/fs);    % number of seconds

% features to be computed
switch(method)
    case 'psd'
        featNames={'maxNormPSD'};
    case 'tree'
        error('tree method not yet implemented')
        featNames={'maxNormPSD','stdNormPSD','power'};
    otherwise
        error('Unknown method: %s',method)
end

% ---- COMPUTE FEATURES ----
Nfeat = length(featNames);
featVals = nan(Ns,Nfeat);
for si=1:Ns
    inds = (si-1)*fs+1 : min(si*fs, N);    
    featVals(si,:) = sigInspectComputeFeatures(signal(inds),featNames,fs);
end


% ---- CLASSIFY ----
switch(method)
    case 'psd'
        % compare to preset threshold
        annot = featVals(:,1)'>0.0123;
    case 'tree'
        % classify using decision tree
        error('Tree not yet implemented')
        featNames={'maxNormPSD','stdNormPSD','power'};
end

