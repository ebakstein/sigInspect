% sigInspectDataInterface
% base interface for loading data incl. metadata to mer signal inspection
% tool
%
% abstract handle-type class
% E. Bakstein 2015-05-29

classdef sigInspectDataInterface < handle 
   properties(Abstract)
%       signalIds  % replaced by getSignalIds - can be private
       
%        settings;
       % settings for sigInspect GUI - any set field will overwrite
       % defaults
   end   
   methods(Abstract)      
%       [signals, labels] = getSignalsById(signalId)
      [signals chInfo] = getSignalsById(obj,signalId)
      % returns a matrix of corresponding signals for given signalId (e.g. 
      % parallel signals from a multielectrode etc.)
      % signals are in matrix in rows, time samples are columns
      % chInfo (string) may provide info about channels

     [signalIds] = getSignalIds(obj)
      % returns list of signal identifiers - a cell array of strings, each signalId identifying
      % all signals/channels to be viewed at once (e.g. parallel signals from a multielectrode)   


   end
   properties
       settings=struct();
       % settings for sigInspect GUI - any set field will overwrite
       % defaults
   end      
end