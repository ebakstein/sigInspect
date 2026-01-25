function pId = getPatientId(signalId)
    % Extracts Patient ID from signal name using specific patterns.
    % Updated to match the logic from "getPatientIds"
    
    % Check for 'ch' followed by digits at the end
    chxPattern = regexp(signalId, 'ch\d+$', 'once');
    
    if ~isempty(chxPattern)
        % Case 1: 'chx' pattern found (e.g., 'sig_Kos_JurM_Dex_1_1_ch1')
        % Extract characters from 5th to 12th (inclusive)
        if length(signalId) >= 12
            pId = signalId(5:12);
        else
            pId = ''; 
            warning('Signal ID "%s" is too short to extract patient ID.', signalId);
        end
    else
        % Case 2: 'chx' pattern not found (e.g., 'sig_59t146p4x')
        % Extract substring before the first 't'
        if contains(signalId, 't')
            pId = extractBefore(signalId, 't');
        else
            % Fallback if 't' is missing to prevent crash
            pId = signalId;
        end
    end
end

