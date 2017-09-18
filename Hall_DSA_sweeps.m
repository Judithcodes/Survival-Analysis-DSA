% Survival analysis-based dynamic spectrum access algorithm
% Initial framework, with sweeps
%
% Based on 2017 conference paper by T.A. Hall et al.
%--------------------------------------------------------------------------

% Training algorithm with spectrum occupancy data representative of channel
% characteristics

% Simulation variables
length = 10000;                    % number of samples in each channel of spectrum occupancy data
t = 0;                            % time marker
tau = 1;                          % transmit duration requested
threshold = 0.50;                  % interference threshold (probability of successful transmission)
startP = 0;                     % starting spectrum occupancy density
stopP = 100;                     % ending spectrum occupancy density
sweepsP = 101;                   % number of sweeps for spectrum occupancy
startQ = 1;                     % starting transmit request density
stopQ = 10;                     % ending transmit request density
sweepsQ = 10;                   % number of sweeps for request density

%decision = zeros(channels, length);
transTot = zeros(sweepsP, sweepsQ);         % segments where secondary user successfully transmits
interfTot = zeros(sweepsP, sweepsQ);        % segments where secondary user collides with primary user
util = zeros(sweepsP, sweepsQ);
interfRate = zeros(sweepsP, sweepsQ);

%--------------------------------------------------------------------------------------------------
% Train DSA algorithm
%--------------------------------------------------------------------------------------------------
for p = linspace(startP, stopP, sweepsP)
    x = p+1;
    % Occupancy data
    m = 1.0;
    b = 0.06;
    counts = zeros(1, length);    % stores number of occurences of each idle period length
    %----------------------------------------------------------------------------
    % Variant 1: Generate test array of single channel with random occupancy
    %----------------------------------------------------------------------------
    trainer = zeros(1, length);
    M = zeros(1, length);
    for k = 1:length
        roll1 = stopP * rand;
        roll2 = stopP * rand;
        if p >= roll1
            trainer(k) = 1;
        elseif p < roll1
            trainer(k) = 0;
        end
        if p >= roll2
            M(k) = 1;
        elseif p < roll2
            M(k) = 0;
        end
    end
    %----------------------------------------------------------------------------
    % Variant 2: Randomly generated occupancy, exponential 2
    %----------------------------------------------------------------------------
    % trainer = spectrum_occ_exp(1, length, m, b);
    % M = spectrum_occ_exp(1, length, m, b);
    %----------------------------------------------------------------------------
    % Variant 3: Periodic spectrum occupancy
    %----------------------------------------------------------------------------
%     trainer = [ones(1, p), zeros(1, stopP - p)];
%     trainer = repmat(trainer, 1, length/stopP);
%     M = trainer;
    %----------------------------------------------------------------------------
    occupied = sum(M);
    vacant = length - occupied;
    
    t = 0;
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
    cdf = cumsum(pdf, 'reverse');

    % Scan test matrix of occupancy data and grant or deny transmission
    % requests

    for q = linspace(startQ, stopQ, sweepsQ)
        y = q;
        t = 0;
        
        % Secondary user transmit request scheduling
        %--------------------------------------------------------------------------
        % Variant 1: Periodic secondary transmit request
        %--------------------------------------------------------------------------
%         period2nd = 10;                   % period for secondary user transmit
%         requests = [zeros(1, stopQ - q), ones(1, q)];
%         requests = repmat(requests, 1, length/stopQ); 
        %--------------------------------------------------------------------------
        % Variant 2: Random secondary transmit request
        %--------------------------------------------------------------------------
        requests = zeros(1, length);
        for k = 1:length
            roll = stopQ * rand;
            if q >= roll
                requests(k) = 1;
            elseif q < roll
                requests(k) = 0;
            end
        end
        schedule = zeros(1, length + 10);         % transmit grant schedule for secondary user
        transmit = zeros(1, length);
        interfere = zeros(1, length);
        
        for i = 1:length
            sample = M(i);
            if sample == 0
                t = t + 1;
                if schedule(i) == 1
                    transmit(i) = 1;
                end
                if requests(i) == 1
                    %-------------------------------------------------------------
                    % Algorithm 1
                    %-------------------------------------------------------------
%                     T = t + tau;
%                     if T > length
%                         T = length; 
%                     end
%                     if cdf(T) >= threshold
%                         schedule((i + 1) : (i + tau)) = 1;
%                     end
                    %-------------------------------------------------------------
                    % Algorithm 2
                    %-------------------------------------------------------------
                    tau = find_duration(t, threshold, cdf);
                    schedule((i + 1) : (i + tau)) = 1;
                    %-------------------------------------------------------------
                end
            elseif sample == 1
                t = 0;
                if schedule(i) == 1
                    interfere(i) = 1;
                end
            end
        end
        
        % Calculate metrics
        transTot(x, y) = sum(transmit);
        util(x, y) = transTot(x, y) ./ vacant;
        interfTot(x, y) = sum(interfere);
        interfRate(x, y) = interfTot(x, y) ./ length;
    end    
end
