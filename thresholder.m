function [ M ] = thresholder( X, threshold )
% Spectrum occupancy thresholding function
%   Takes array of power spectrum density data from RF Explorer or other
%   spectrum analyzer and converts it into array of binary spectrum
%   occupancy values

    L = size(X, 1);
    W = size(X, 2);
    M = zeros(W, L);

    for i = 1:L
        for j = 1:W
           if X(i, j) >= threshold
              M(j, i) = 1; 
           end
        end
    end
end

