% Survival analysis-based dynamic spectrum access algorithm
% Initial framework
%
% Based on 2017 conference paper by T.A. Hall et al.
%--------------------------------------------------------------------------

% Training algorithm with spectrum occupancy data representative of channel
% characteristics
length = 100;
trainer = spectrum_occ(1, length);
counts = zeros(1, length);    % stores number of occurences of each idle period length
t = 0;
tau = 1;         % transmit duration requested

for i = 1:length
    if trainer(i) == 0
        t = t+1;
        if (i + 1) > length
            counts(t) = counts(t) + 1;
        else
            if trainer(i + 1) == 1
                counts(t) = counts(t) + 1;
            end
        end
    elseif trainer(i) == 1        
        t = 0;
    end
end

periodsIdle = sum(counts);
pdf = counts./periodsIdle;

for i = 1:length
    cdf(i) = sum(pdf(1:i));
end
