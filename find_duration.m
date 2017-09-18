function maxDuration = find_duration( t, p, cdf )
%Request maximum channel availability
%   t = time elapsed since end of last transmission
%   p = probability of successful transmission
%   cdf = cumulative distribution function

length = size(cdf, 2);

for i = t:length
    if cdf(i) < p
        maxDuration = i - 1 - t;
        if maxDuration < 0
           maxDuration = 0; 
        end
        break
    else
        maxDuration = 0;
    end
end

end

