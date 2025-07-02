function f = compSignalLength(segment)
    f = sum(abs(diff(segment,1,2)), 2); 
end