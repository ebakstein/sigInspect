function entropyVal = compShannonEntropy(segment, numBins)
    % Random signals large entropy
    % H=−∑pi*log2pi
    if nargin < 2
        numBins = 5;
    end
    energy = segment .^ 2;
    probDist = histcounts(energy, numBins, 'Normalization', 'probability');
    probDist(probDist == 0) = []; % Remove zeros to avoid log(0)
    entropyVal = -sum(probDist .* log2(probDist));
end