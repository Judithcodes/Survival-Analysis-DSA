function maxDuration = find_duration( t, p, cdf )
%Request maximum channel availability
%   t = time elapsed since end of last transmission
%   p = probability of successful transmission
%   cdf = cumulative distribution function

length = size(cdf, 2);
theta = (-1)*log(p);
q = cdf(t);

% for i = t:length
%     if cdf(i) < p
%         maxDuration = i - 1 - t;
%         if maxDuration < 0
%            maxDuration = 0; 
%         end
%         break
%     elseif cdf(i) >= p
%         maxDuration = 10;  
%     else
%         maxDuration = 0;
%     end
% end
%--------------------------------------------------------------
% maxDuration = 0;
% for i = t:(t+50)
%    if cdf(i) >= p
%        maxDuration = i - t;
%    elseif cdf(i) < p
%        break
%    end
% end
%--------------------------------------------------------------
maxDuration = 0;
for i = t:length
    if (cdf(i)-cdf(t)) <= theta
        maxDuration = i - t;
    elseif (cdf(i)-cdf(t)) > theta
        break
    end
end

end

