function irregularity = compIrregularityFactor(segment)
    % Measures how different the signal is from white noise
    complexity_signal = compHjorthComplexity(segment);
    white_noise = randn(size(segment)); % Simulated white noise
    complexity_noise = compHjorthComplexity(white_noise);
    irregularity = complexity_signal / (complexity_noise + eps);
end
