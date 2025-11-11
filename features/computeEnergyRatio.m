function energyRatio = computeEnergyRatio(signal, fs, spikeBand)
    if nargin < 3
        spikeBand = [300 3000];
    end

    % Always ensure single-channel vector
    signal = signal(:);

    % Compute PSD (periodogram returns single-sided by default)
    [Pxx, F] = periodogram(signal, [], [], fs);

    % Compute power from PSD directly (F and Pxx are aligned)
    spikePower = bandpower(Pxx, F, spikeBand, 'psd');
    totalPower = bandpower(Pxx, F, 'psd');

    % Avoid division by zero
    if totalPower == 0
        energyRatio = 0;
    else
        energyRatio = spikePower / totalPower;
    end
end
