% Survival analysis-based dynamic spectrum access algorithm
% Modifying length of training sequence
%
% Based on 2017 journal and conference paper by T.A. Hall et al.
%--------------------------------------------------------------------------

% Training algorithm with spectrum occupancy data representative of channel
% characteristics

% Simulation variables
Length = 1000000;                    % number of samples in each channel of spectrum occupancy data
channels = 1;                   % number of channels in test occupancy matrix
t = 0;                            % time marker
tau = 10;                          % transmit duration requested
threshold = 0.90;                  % interference threshold (probability of successful transmission)
theta = (-1)*log(threshold);

% Occupancy data
P1 = 20;             % Occupancy event rate (lambda)
P2 = 20;             % Vacancy event rate
S1 = 15;            % SU request event rate
S2 = 15;            % SU idle event rate
%=============================================================================
% Variant 1: Randomly generated occupancy, exponential
%=============================================================================
% trainer = spectrum_occ_exp(1, Length, P2, P1);          % training array for DSA algorithm
% M = spectrum_occ_exp(channels, Length, P2, P1);         % test matrix of occupancy data
%=============================================================================
% Variant 2: Randomly generated occupancy, dual poisson processes
%=============================================================================
% M = spectrum_occ_poiss(channels, Length, P1, P2);
% trainer = spectrum_occ_poiss(1, Length, P1, P2);
%=============================================================================
% Variant 3: Periodic spectrum occupancy
%=============================================================================
% duty1st = 0.3;
% period1st = 10;
% trainer = [ones(1, period1st * duty1st), zeros(1, period1st - (period1st * duty1st))];
% trainer = repmat(trainer, 1, Length/period1st);
% M = [ones(channels, period1st * duty1st), zeros(channels, period1st - (period1st * duty1st))];
% M = repmat(M, 1, Length/period1st);
%=============================================================================
% Variant 4: Gaussian distributed spectrum occupancy
%=============================================================================
% M = spectrum_occ(channels, Length);
% trainer = spectrum_occ(1, Length);
%=============================================================================
% Variant 5: Time-varying occupancy, Poisson
%=============================================================================
trainer = spectrum_occ_poiss(1, Length, P1, P2);
Tw = 5000;
nWin = floor(Length/Tw);
m = 1:nWin;

% Sinusoidal channel parameter variations
%----------------------------------------
P1 = 20 - 15*sin(2*pi*(m-1)/100);             % Occupancy event rate
P2 = 20 + 15*sin(2*pi*(m-1)/100);            % Vacancy event rate

% Linear variation
%----------------------
% P2 = 2 + abs(100-m);
% P1(1:100) = 2 + m(1:100);
% P1(101:200) = 102 - abs(100 - m(101:200));

M = [];
for i = 1:nWin
   M =  [M, spectrum_occ_poiss(channels, Tw, P1(i), P2(i))];
end
%=============================================================================
% Variant 6: Time-varying occupancy, exponential
%=============================================================================
% trainer = spectrum_occ_exp(1, Length, P2, P1);
% Tw = 5000;
% nWin = floor(Length/Tw);
% m = 1:nWin;
% 
% % Sinusoidal variation
% %----------------------
% % P2 = 20 + 15*sin(2*pi*m/10);            % Vacancy event rate 
% 
% % Linear variation
% %----------------------
% P2 = abs(100-8*m);
% 
% M = [];
% for i = 1:nWin
%     M =  [M, spectrum_occ_exp(channels, Tw, P2(i), P1)];
% end
%-----------------------------------------------------------------------------

occupied = sum(M);
vacant = Length - occupied;

% Secondary user transmit request scheduling
%==========================================================================================
% Variant 1: Periodic SU transmit request
%==========================================================================================
% duty2nd = 1;                     % duty cycle for secondary user transmit
% period2nd = 10;                   % period for secondary user transmit
% requests = [zeros(1, period2nd - period2nd*duty2nd), ones(1, period2nd*duty2nd)];
% requests = repmat(requests, 1, Length/period2nd); 
% requests = repmat(requests, channels, 1);       % transmit request schedule for secondary user
%==========================================================================================
% Variant 2: Poisson distributed SU transmit request
%==========================================================================================
requests = spectrum_occ_poiss(channels, Length, S1, S2);
%------------------------------------------------------------------------------------------
schedule = zeros(channels, Length + 100);         % transmit grant schedule for secondary user
transmit = zeros(channels, Length);         % segments where secondary user successfully transmits
interfere = zeros(channels, Length);        % segments where secondary user collides with primary user


%=====================================================================================================
% Train DSA algorithm
%=====================================================================================================
% Generate array with number of instances of each length of idle period in
% training vector
counts = occupancy(trainer);

n = length(find(counts));

% Calculate survival/hazard function
periodsIdle = sum(counts);
pdf = counts./periodsIdle;
cdf = cumsum(pdf);
ccdf = cumsum(pdf, 'reverse');
%=============================================================================
% Cumulative Hazard Function
%=============================================================================
Ti = [];
for i = 1:Length
    Ti = [Ti, i*ones(1, counts(i))];
end

h = zeros(1, Length);
H = zeros(1, Length);
j = 1;
for t = 1:n
    temp = 0;    
    while t >= Ti(j)
        temp = temp + 1/(periodsIdle - j + 1);
        j = j + 1;
        if j > periodsIdle
            j = periodsIdle;
            break
        end
    end
    h(t) = temp;
    if t == 1
        H(t) = temp;
    elseif t == n   
        H(t:Length) = H(t-1) + temp;
    else
        H(t) = H(t-1) + temp;
    end    
end
%-----------------------------------------------------------------------------

% Scan test matrix of occupancy data and grant or deny transmission
% requests
for i = 1:channels
    t = 0;
    for j = 1:Length
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
                T = t + tau;
                if T > Length
                    T = Length; 
                end
                if H(T) - H(t) < theta
                    schedule(i, (j + 1) : (j + tau)) = 1;
                end
                %-------------------------------------------------------------
%                 T = t + tau;
%                 if T > Length
%                     T = Length; 
%                 end
%                 if ccdf(T) >= threshold
%                     schedule(i, (j + 1) : (j + tau)) = 1;
%                 end
                %=============================================================
                % Algorithm 2
                %=============================================================
%                 tau = 1;
%                 while (H(t + tau)) < theta
%                     tau = tau + 1;
%                 end
%                 tau = tau - 1;
%                 schedule(i, (j + 1) : (j + 1 + tau)) = 1;
                %-------------------------------------------------------------    
%                 tau = 1;
%                 while (H(t + tau) - H(t)) < theta
%                     tau = tau + 1;
%                     if (t + tau) > Ti(periodsIdle)
%                        break
%                     end
%                 end
%                 tau = tau - 1;
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
util = 100 * transTot./vacant;
interfTot = sum(interfere, 2);
interfRate = 100 * interfTot./(Length);

