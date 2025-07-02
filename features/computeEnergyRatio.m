function energyRatio = computeEnergyRatio(signal, fs, spikeBand)
    % Computes the energy ratio in a specified spike band relative to total power
    if nargin < 3
        spikeBand = [300 3000]; % Default spike band if not provided
    end

    % Compute power in the spike band
    spikePower = bandpower(signal, fs, spikeBand);  
    % Compute total power of the signal
    totalPower = bandpower(signal, fs, [0 fs/2]); 

    % Prevent division by zero
    if totalPower == 0
        energyRatio = 0;
    else
        energyRatio = spikePower / totalPower;
    end
end