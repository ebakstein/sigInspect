% implements data loading from *.mat file or from a cell array or matrix
%
% USE
%   MATRIX/CELL ARRAY: sigInspectDataBasic(signals, samplingFreq)
%       signals - signals as a matrix with parallel recordings in rows.
%                 multiple signals can be provided as a cell array
%       samplingFreq - sampling frequency in Hz (default:24000Hz)
% 
%   MAT-FILES: sigInspectDataBasic(pathToMat, samplingFreq)
%       Searches for a variable named 'signal','data' or 'signals', 
%       which is a matrix (signals in rows) or a cell array of such matrices
%       if there is a variable samplingFreq or fs in the file, it is used as
%       sampling frequency. Otherwise, default sampling freq. 24000 is used
%
% E. Bakstein 2015-06-01
% 
classdef sigInspectDataBasic < sigInspectDataInterface
  
    
    % -------------- properties --------------
  
    properties
        signalIds={};
        signals={};
    end
                      

    methods
        % constructor 
        
        function obj=sigInspectDataBasic(pathOrMatrix, samplingFreq, signalIds)
            if(nargin<1)
                error('provide path to a mat file to load or cell array of signals as arguments');
            end
            
            if(ischar(pathOrMatrix))
                 % path to a mat file
                 obj.loadMat(pathOrMatrix);                 
            elseif(iscell(pathOrMatrix))                
                if(~obj.checkDataCell(pathOrMatrix))
                    error('first input must be a cell array of non-empty numeric vectors/matrices')                        
                end
                obj.signals=pathOrMatrix;                
            
            elseif(isnumeric(pathOrMatrix))
                % encapsulate the input into 
                obj.signals={pathOrMatrix};
            end
            
            % 2nd argument - samplingFreq
            if(nargin>1 && ~isempty(samplingFreq))
                if(~isnumeric(samplingFreq) && samplingFreq>0)
                    warning('Sampling frequency (second argument) must be a (positive) number - ignoring!')
                else
                    obj.settings.SAMPLING_FREQ=samplingFreq;
                end
            end
            
            % 3rd argument - signalIds
            if(nargin>2)
                if(iscell(signalIds))
                    obj.signalIds=signalIds; 
                else
                    warning('signal ids must be a cell array of strings - ingnoring, using default numbering')
                    obj.generateIds(); % use default ids
                end
            else
                obj.generateIds();
            end 
            
            % determine max. channel count from data
            obj.settings.PLOT_CHANNELS = max(cellfun(@(x)size(x,1),obj.signals));
            
        end
        
        
        % return a signal stored in signals variable
        function [signals chInfo] = getSignalsById(obj,signalId)            
            ind=find(strcmp(obj.signalIds,signalId));
            if(isempty(ind))
                error('invalid/unknown signal id: %s',signalId)
            end
            signals = obj.signals{ind};            
            
            % channel labels: no info, returning empty chars
            if(nargout>1)
                chInfo='';
            end
        end
        
        % return list of signal ids
        function signalIds = getSignalIds(obj)
            if(isempty(obj.signalIds))
                error('signal Ids not initialized yet')
            end
            signalIds = obj.signalIds;            
        end
        
        
        % load signals from a mat file
        function loadMat(obj,filePath)
            if(~exist(filePath,'file'))
                error('file does not exist: %s',filePath)
            end
            
            if(isempty(strfind(filePath,'.mat')))
                error('only files *.mat supported by sigInspectDataBasic, for *.csv or *.txt run sigInspectDataCsv')
            end

            tmp=load(filePath);
            chckNames={'signal','data','signals'}; % variable names to be checked
            fnd=0;
            % check possible names of data/signals variable
            for ci=1:length(chckNames)
                if(isfield(tmp,chckNames{ci})) % field exists
                    var=tmp.(chckNames{ci});
                    if(iscell(var)) % it is a cell(array?)
                        ok = obj.checkDataCell(var);
                        if(~ok)
                            warning('file %s contains cell array %s, which is not compliant > must be a cell array of nonempty matrices/vectors - skipping',filePath,chckNames{ci})
                            continue
                        end
                        obj.signals=var;   
                        fnd=1;
                        break
                    elseif(ismatrix(var) && isnumeric(var))
                        obj.signals={var};
                        fnd=1;
                        break;
                    else                       
                        warning('file %s contains cell array %s, is not a matrix/not numeric - skipping',filePath,chckNames{ci})
                        continue
                    end
                end
            end
            if(~fnd)
                error('No signals variable (data, signal or signals) could be found in the file %s',filePath)
            end
        end
                                                                        
        % generate signalIds (cell array of strings '1' to length(signals)),
        % put them to signalIds property
        function generateIds(obj)
            if(isempty(obj.signals))
                warning('No signals loaded > no ids')
                return
            end
            obj.signalIds = arrayfun(@(x)['#' num2str(x)], 1:length(obj.signals), 'unif', 0);
        end
    end

    
    
    % -------------- static methods --------------
    
    
    
    methods(Static)
       % check data cell array - should contain only non-empty numerical matrices
       function ok = checkDataCell(dataCell)
           if(~iscell(dataCell))
               ok=false;
               return
           end
           allOk = cellfun(@(x) ismatrix(x) && isnumeric(x) && ~isempty(x),dataCell);                
           ok=all(allOk);
      end
      function ok = checkIdsCell(idsCell)
           if(~iscell(idsCell))
               ok=false;
               return
           end      
          allOk = cellfun(@(x) ischar(x) && ~isempty(x),idsCell);
          ok=all(allOk) && length(unique(idsCell))==length(idsCell);
      end   
    end
end