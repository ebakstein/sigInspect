function sparseness = compSparseness(segment)
    % Measures how concentrated the energy is in few large values
    N = length(segment);
    % RMS/Mean absolute value
    sparseness = sqrt(sum(segment.^2) / N) / (sum(abs(segment)) / N);
end
