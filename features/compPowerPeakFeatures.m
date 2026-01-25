function [numPeaks, meanPeakHeight, peakFreq, peakRMSRatio, avgPeakWidth] = compPowerPeakFeatures(signal, fs, smoothWindowDuration, plotOption)
    % signal: channels x time

    [numChannels, numSamples] = size(signal);
    powerSignal = signal.^2;

    % Gaussian smoothing
    windowSize = round(smoothWindowDuration * fs);
    gaussKernel = gausswin(windowSize);
    gaussKernel = gaussKernel(:) / sum(gaussKernel);  % column vector

    if numSamples < 3 * windowSize
        % Signal too short for filtfilt, fall back to conv2 (same as V1 fallback)
        smoothSignal = conv2(powerSignal, gaussKernel', 'same');
    else
        % Use filtfilt on transposed signal (Time x Channels) to match V1
        % We transpose powerSignal so filtfilt works along the time dimension (columns)
        smoothSignal = filtfilt(gaussKernel, 1, powerSignal')'; 
    end
    
    % Initialize outputs
    numPeaks        = zeros(numChannels, 1);
    meanPeakHeight  = zeros(numChannels, 1);
    peakFreq        = zeros(numChannels, 1);
    peakRMSRatio    = zeros(numChannels, 1);
    avgPeakWidth    = zeros(numChannels, 1);

    minDistance = round(fs * 0.004); % in samples

    for ch = 1:numChannels
        smSig = smoothSignal(ch, :);

        if numSamples <= minDistance
            warning('Signal too short for peak detection in channel %d.', ch);
            continue;
        end

        % Adaptive thresholding
        baselineNoise = median(abs(smSig - median(smSig)));
        minProm = 1.5 * baselineNoise;
        minHeight = median(smSig) + 8 * baselineNoise;

        % Find peaks
        [pks, locs, width, ~] = findpeaks(smSig, ...
            'MinPeakProminence', minProm, ...
            'MinPeakHeight', minHeight, ...
            'MinPeakDistance', minDistance);

        if isempty(pks)
            continue;
        end

        % Peak-based features
        numPeaks(ch)       = numel(pks);
        meanPeakHeight(ch) = mean(pks, 'omitnan');
        avgPeakWidth(ch)   = mean(width, 'omitnan');

        if numel(locs) > 1
            intervals = diff(locs) / fs;
            peakFreq(ch) = 1 / mean(intervals);
        end

        % RMS ratio
        peakRMSRatio(ch) = meanPeakHeight(ch) / rms(smSig);
    end

    % Plot only the first channel (if requested)
    if plotOption
        timeVector = (0:numSamples-1) / fs;
        smSig = smoothSignal(1, :);
        sig1 = signal(1, :);
        [pks, locs] = findpeaks(smSig, ...
            'MinPeakProminence', 1.5 * median(abs(smSig - median(smSig))), ...
            'MinPeakHeight', median(smSig) + 8 * median(abs(smSig - median(smSig))), ...
            'MinPeakDistance', minDistance);
        figure;
        plot(timeVector, sig1, 'b'); hold on;
        plot(timeVector, smSig, '--r', 'LineWidth', 1);
        scatter(timeVector(locs), pks, 'kx', 'MarkerFaceColor', 'k', 'LineWidth', 1.5);
        xlabel('Time (s)');
        ylabel('Amplitude');
        title(sprintf(['Ch 1 - N Peaks: %d | Mean height: %.2f | Freq: %.2f Hz | Peak RMS Ratio: %.2f\n' ...
            'Mean width: %.2f'], ...
            numPeaks(1), meanPeakHeight(1), peakFreq(1), peakRMSRatio(1), avgPeakWidth(1)));
        legend({'Original signal', 'Smoothed power signal', 'Peaks'});
        grid on;
    end
end
