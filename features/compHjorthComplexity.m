function complexity = compHjorthComplexity(segment)
    % How much the frequency content changes over time
    % High for irregular signal
    deriv1 = diff(segment);
    deriv2 = diff(deriv1);
    mobility1 = compHjorthMobility(deriv1);
    mobility2 = compHjorthMobility(deriv2);
    complexity = mobility2 / (mobility1 + eps);
end