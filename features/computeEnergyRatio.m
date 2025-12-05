function energyRatio = computeEnergyRatio(signal, fs, spikeBand)
    if nargin < 3
        spikeBand = [300 3000];
    end

    signal = signal(:);
    N = length(signal);

    % Compute band power
    spikePower = bandpower(signal, fs, spikeBand);
    totalPower = bandpower(signal, fs, [0 min(fs/2, (N - 1) * fs/(2 * N))]);

    % Avoid division by zero
    if totalPower == 0
        energyRatio = 0;
    else
        energyRatio = spikePower / totalPower;
    end
end
