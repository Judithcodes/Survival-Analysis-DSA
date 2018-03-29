function [ M ] = thresholder2( X )
% Spectrum occupancy thresholding function
%   Takes array of power spectrum density data from RF Explorer or other
%   spectrum analyzer and converts it into array of binary spectrum
%   occupancy values

    L = size(X, 1);
    W = size(X, 2);
    M = zeros(W, L);
    
    noiseFloor = min(X);
    threshold = noiseFloor + 3;

    for i = 1:L
        for j = 1:W
           if X(i, j) >= threshold(j)
              M(j, i) = 1; 
           end
        end
    end
end

