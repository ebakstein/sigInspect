function mobility = compHjorthMobility(segment)
    % Frequency content and smoothness of the signal
    % High for random noise
    deriv = diff(segment);
    mobility = sqrt(var(deriv) / (var(segment) + eps));
end