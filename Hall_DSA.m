% Survival analysis-based dynamic spectrum access algorithm
% Initial framework
%
% Based on 2017 conference paper by T.A. Hall et al.
%--------------------------------------------------------------------------

% Training algorithm with spectrum occupancy data representative of channel
% characteristics

% Simulation variables
length = 100;                    % number of samples in each channel of spectrum occupancy data
channels = 10;                   % number of channels in test occupancy matrix
counts = zeros(1, length);    % stores number of occurences of each idle period length
t = 0;                            % time marker
tau = 1;                          % transmit duration requested

% Occupancy data
trainer = spectrum_occ(1, length);          % training array for DSA algorithm
M = spectrum_occ(channels, length);         % test matrix of occupancy data
occupied = sum(M, 2);
vacant = length - occupied;

% Secondary user transmit request scheduling
duty2nd = 10;                     % duty cycle for secondary user transmit
period2nd = 10;                   % period for secondary user transmit
schedule2nd = [zeros(1, 100-duty2nd), ones(1, duty2nd)];
schedule2nd = repmat(schedule2nd, 1, length/period);     % transmit schedule for secondary user
decision = zeros(channels, length/period2nd);

%--------------------------------------------------------------------------------------------------
% Train DSA algorithm
%--------------------------------------------------------------------------------------------------
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

% Calculate survival/hazard function
periodsIdle = sum(counts);
pdf = counts./periodsIdle;
cdf = zeros(1, length);

for i = 1:length
    cdf(i) = sum(pdf(i:length));
end

% Scan test matrix of occupancy data and grant or deny transmission
% requests
for i = 1:channels
    t = 0;
    for j = 1:samples
        sample = M(i, j);
        if sample == 0
            t = t + 1;
            if schedule2nd == 1
                
            end
        elseif sample == 1
            t = 0;
        end
    end
end
