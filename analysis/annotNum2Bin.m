function ba = annotNum2Bin(numAnnot, maxN)
    % Convert numerical annotations to binary vector
    if nargin < 2
        maxN = 5; % Default number of artifact types
    end
    
    na = numAnnot(:);  % Flatten the input annotations
    ba = false(length(na), maxN);  % Initialize binary annotation matrix
    
    for ii = 1:length(na)
        x = na(ii);
        for ni = maxN:-1:0        
            if x >= 2^ni
                x = x - 2^ni;
                ba(ii, ni + 1) = true;
            end
        end
    end
end
