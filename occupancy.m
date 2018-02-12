function [ countsIdle, countsOcc ] = occupancy( M )
%Spectrum occupancy test array processor
%   Takes array of binary spectrum occupancy values and returns two arrays
%   containing the number of each length of occupied and idle period present.
%       [ countsIdle, countsOcc ] = occupancy( M )
    D = size(M);
    Len = D(2);     % number of samples (length)
    Chan = D(1);     % number of channels (height)
    t1 = zeros(Chan, 1);
    t2 = zeros(Chan, 1);
    countsIdle = zeros(Chan, Len);
    countsOcc = zeros(Chan, Len);
    
    for j = 1:Len
       for i = 1:Chan
           if M(i, j) == 0
               t2(i) = 0;
               t1(i) = t1(i) + 1;
               if ((j+1) > Len) || (M(i, (j+1)) == 1)
                   countsIdle(i, t1(i)) = countsIdle(i, t1(i)) + 1;
               end
           elseif M(i, j) == 1
               t1(i) = 0;
               t2(i) = t2(i) + 1;
               if ((j+1) > Len) || (M(i, (j+1)) == 0)
                   countsOcc(i, t2(i)) = countsIdle(i, t2(i)) + 1;
               end
           end
       end
    end

    
end

