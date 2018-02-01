% Survival analysis-based dynamic spectrum access algorithm
% Initial framework
%
% Based on 2017 conference paper by T.A. Hall et al.
%--------------------------------------------------------------------------

% Training algorithm with spectrum occupancy data representative of channel
% characteristics

% Simulation variables
length = 10000;                    % number of samples in each channel of spectrum occupancy data
channels = 1;                   % number of channels in test occupancy matrix
counts = zeros(1, length);    % stores number of occurences of each idle period length
t = 0;                            % time marker
tau = 1;                          % transmit duration requested
threshold = 0.50;                  % interference threshold (probability of successful transmission)

% Occupancy data
L1 = 20;             % Occupancy event rate (lambda)
L2 = 20;             % Vacancy event rate
%=============================================================================
% Variant 1: Randomly generated occupancy, exponential
%=============================================================================
% trainer = spectrum_occ(1, length);          % training array for DSA algorithm
% M = spectrum_occ(channels, length);         % test matrix of occupancy data
%=============================================================================
% Variant 2: Randomly generated occupancy, dual poisson processes
%=============================================================================
M = spectrum_occ_poiss(channels, length, L1, L2);
trainer = M(1, :);
%=============================================================================
% Variant 3: Periodic spectrum occupancy
%=============================================================================
% duty1st = 0.3;
% period1st = 10;
% trainer = [ones(1, period1st * duty1st), zeros(1, period1st - (period1st * duty1st))];
% trainer = repmat(trainer, 1, length/period1st);
% M = [ones(channels, period1st * duty1st), zeros(channels, period1st - (period1st * duty1st))];
% M = repmat(M, 1, length/period1st);
%----------------------------------------------------------------------------
occupied = sum(M, 2);
vacant = length - occupied;

% Secondary user transmit request scheduling
duty2nd = 1;                     % duty cycle for secondary user transmit
period2nd = 10;                   % period for secondary user transmit
requests = [zeros(1, period2nd - period2nd*duty2nd), ones(1, period2nd*duty2nd)];
requests = repmat(requests, 1, length/period2nd); 
requests = repmat(requests, channels, 1);       % transmit request schedule for secondary user
schedule = zeros(channels, length + 100);         % transmit grant schedule for secondary user
%decision = zeros(channels, length);
transmit = zeros(channels, length);         % segments where secondary user successfully transmits
interfere = zeros(channels, length);        % segments where secondary user collides with primary user

%=====================================================================================================
% Train DSA algorithm
%=====================================================================================================
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
cdf = cumsum(pdf);
ccdf = cumsum(pdf, 'reverse');
%=============================================================================
% Cumulative Hazard Function
%=============================================================================
H = zeros(1, length);
H(1) = counts(1)*(1/periodsIdle);
for j = 2:length
    if periodsIdle >= j
        H(j) = H(j-1) + counts(j)*(1/(periodsIdle - j + 1));
    else
        H(j:length) = H(j-1);
        break
    end
end
theta = (-1)*log(threshold);
%-----------------------------------------------------------------------------

% Scan test matrix of occupancy data and grant or deny transmission
% requests
for i = 1:channels
    t = 0;
    for j = 1:length
        sample = M(i, j);
        if sample == 0
            t = t + 1;
            if schedule(i, j) == 1
                transmit(i, j) = 1;
            end
            if requests(i, j) == 1
                %=============================================================
                % Algorithm 1
                %=============================================================
%                 T = t + tau;
%                 if T > length
%                     T = length; 
%                 end
%                 if H(T) < threshold
%                     schedule(i, (j + 1) : (j + tau)) = 1;
%                 end
% %                 if cdf(T) >= threshold
% %                     schedule(i, (j + 1) : (j + tau)) = 1;
% %                 end

                %=============================================================
                % Algorithm 2
                %=============================================================
                tau = 1;
                while (H(t + tau)) < threshold
                    tau = tau + 1;
                end
                tau = tau - 1;
                schedule(i, (j + 1) : (j + 1 + tau)) = 1;
                %-------------------------------------------------------------
%                 tau = quantile(cdf, (1 - threshold * ccdf(t)));
%                 schedule(i, (j + 1) : (j + 1 + tau)) = 1;
                %-------------------------------------------------------------
            end
        elseif sample == 1
            t = 0;
            if schedule(i, j) == 1
                interfere(i, j) = 1;
            end
        end
    end
end

% Calculate metrics
transTot = sum(transmit, 2);
util = transTot./vacant;
interfTot = sum(interfere, 2);
interfRate = interfTot./(length);
