% Survival analysis-based dynamic spectrum access algorithm
% Initial framework, with sweeps
% Sweeping transmit request period
%
% Based on 2017 conference paper by T.A. Hall et al.
%--------------------------------------------------------------------------

% Training algorithm with spectrum occupancy data representative of channel
% characteristics

% Simulation variables
Length = 10000;                    % number of samples in each channel of spectrum occupancy data
t = 0;                            % time marker
threshold = 0.9;                  % interference threshold (probability of successful transmission)
theta = (-1)*log(threshold);
P1 = 15;             % Occupancy event rate (lambda)
P2 = 15;             % Vacancy event rate
S1 = 15;            % SU request event rate
S2 = 15;            % SU idle event rate
startP = 15;                     % starting spectrum occupancy density
stopP = 15;                     % ending spectrum occupancy density
sweepsP = 1;                   % number of sweeps for spectrum occupancy
startQ = 1;                     % starting transmit request duration
stopQ = 10;                     % ending transmit request duration
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
    x = p;
    % Occupancy data
    P1 = 3;             % Occupancy event rate (lambda)
    P2 = 15;             % Vacancy event rate
    counts = zeros(1, Length);    % stores number of occurences of each idle period length
    %----------------------------------------------------------------------------
    % Variant 1: Generate test array of single channel with random occupancy
    %----------------------------------------------------------------------------
%     trainer = zeros(1, Length);
%     M = zeros(1, Length);
%     for k = 1:Length
%         roll1 = stopP * rand;
%         roll2 = stopP * rand;
%         if p >= roll1
%             trainer(k) = 1;
%         elseif p < roll1
%             trainer(k) = 0;
%         end
%         if p >= roll2
%             M(k) = 1;
%         elseif p < roll2
%             M(k) = 0;
%         end
%     end
    %----------------------------------------------------------------------------
    % Variant 2: Randomly generated occupancy, dual Poisson processes
    %----------------------------------------------------------------------------
    M = spectrum_occ_poiss(1, Length, p, P2);
    trainer = spectrum_occ_poiss(1, Length, p, P2);
    %----------------------------------------------------------------------------
    % Variant 3: Periodic spectrum occupancy
    %----------------------------------------------------------------------------
%     trainer = [ones(1, p), zeros(1, stopP - p)];
%     trainer = repmat(trainer, 1, Length/stopP);
%     M = trainer;
    %----------------------------------------------------------------------------
    occupied = sum(M);
    vacant = Length - occupied;
    
    t = 0;
    for i = 1:Length
        if trainer(i) == 0
            t = t+1;
            if (i + 1) > Length
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
    
    n = length(find(counts));

    % Calculate survival/hazard function
    periodsIdle = sum(counts);
    pdf = counts./periodsIdle;
    cdf = cumsum(pdf, 'reverse');
    %---------------------------------------------------------------------
    % Cumulative Hazard Function
    %---------------------------------------------------------------------
%     H = zeros(1, Length);
%     H(1) = counts(1)*(1/n);
%     for j = 2:Length
%         if periodsIdle >= j
%             H(j) = H(j-1) + counts(j)*(1/(periodsIdle - j + 1));
%         else
%             H(j:Length) = H(j-1);
%             break
%         end
%     end
    %---------------------------------------------------------------------
    Ti = [];
    for i = 1:Length
        Ti = [Ti, i*ones(1, counts(i))];
    end

    H = zeros(1, Length);
    j = 1;
    for t = 1:n
        temp = 0;    
        while t >= Ti(j)
            temp = temp + 1/(periodsIdle - j + 1);
            j = j + 1;
            if j > periodsIdle
                j = periodsIdle;
            end
        end
        if t == 1
            H(t) = temp;
        elseif t == n   
            H(t:Length) = H(t-1) + temp;
        else
            H(t) = H(t-1) + temp;
        end    
    end
%----------------------------------------------------------------------------

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
% %         requests = ones(1, Length);
%         requests = [ 0 0 1 1 ];
%         requests = repmat( requests, 1, Length/4);
        %----------------------------------------------------------------------------
        % Variant 2: Randomly generated occupancy, dual Poisson processes
        %----------------------------------------------------------------------------
        requests = spectrum_occ_poiss(1, Length, S1, S2);
        %--------------------------------------------------------------------------
        % Variant 3: Random secondary transmit request
        %--------------------------------------------------------------------------
%         requests = zeros(1, Length);
%         for k = 1:length
%             roll = stopQ * rand;
%             if q >= roll
%                 requests(k) = 1;
%             elseif q < roll
%                 requests(k) = 0;
%             end
%         end
        %--------------------------------------------------------------------------
        schedule = zeros(1, Length + 100);         % transmit grant schedule for secondary user
        transmit = zeros(1, Length);
        interfere = zeros(1, Length);
        
        for i = 1:Length
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
                    tau = q;
                    T = t + tau;
                    if T > Length
                        T = Length; 
                    end
                    if H(T) - H(t) < theta
                       schedule((i + 1) : (i + tau)) = 1;
                    end
                %-------------------------------------------------------------
                % Algorithm 2
                %-------------------------------------------------------------
%                     tau = 1;
%                     while (H(t + tau) - H(t)) < threshold
%                         tau = tau + 1;
%                         if (t + tau) > Ti(periodsIdle)
%                             break
%                         end
%                     end
%                     tau = tau - 1;
%                     schedule((i + 1) : (i + tau)) = 1;
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
        interfRate(x, y) = interfTot(x, y) ./ Length;
    end    
end
