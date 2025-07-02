function f = compPeakToRMS(segment)
    % Peak to peak vs RMS 
    % high for peaks
    peakToPeak = max(segment,[],2) - min(segment,[],2);
    rmsValue = sqrt(mean(segment.^2,2));
    f = peakToPeak ./ rmsValue;
end