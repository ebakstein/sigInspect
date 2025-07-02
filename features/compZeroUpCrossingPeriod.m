function zup = compZeroUpCrossingPeriod(segment, fs, smoothWindowDuration)
    % Computes the Zero Up-Crossing Period (ZUCP)
    % signal - input signal
    % fs - sampling frequency
    % smoothWindow - number of samples for smoothing (set to 1 for no smoothing)
    
    if smoothWindowDuration > 1
        smoothWindow = round(smoothWindowDuration * fs);
        b = ones(1, smoothWindow) / smoothWindow; % MA
        segment = filtfilt(b, 1, segment);
    end

    % Zero crossings - from negative to positive
    zeroCrossings = find(segment(1:end-1) < 0 & segment(2:end) >= 0);

    % Compute periods between consecutive zero-crossings
    if length(zeroCrossings) < 2
        zup = NaN;
        return;
    end
    periods = diff(zeroCrossings) / fs; % seconds
    zup = mean(periods);
end