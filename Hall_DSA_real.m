% Survival analysis-based dynamic spectrum access algorithm
% Initial framework, with sweeps
%
% Based on 2017 conference paper by T.A. Hall et al.
%--------------------------------------------------------------------------

% Training algorithm with spectrum occupancy data representative of channel
% characteristics

% Simulation variables
t = 0;                            % time marker
tau = 1;                          % transmit duration requested
threshold = 0.9;                  % interference threshold (probability of successful transmission)
theta = (-1)*log(threshold);
% P1 = 15;             % Occupancy event rate (lambda)
% P2 = 15;             % Vacancy event rate
S1 = 15;            % SU request event rate
S2 = 15;            % SU idle event rate
startP = -92;                     % starting spectrum occupancy density
stopP = -76;                     % ending spectrum occupancy density
sweepsP = 17;                   % number of sweeps for spectrum occupancy
startQ = 1;                     % starting transmit request density
stopQ = 1;                     % ending transmit request density
sweepsQ = 1;                   % number of sweeps for request density

util = zeros(sweepsP, sweepsQ);
interfRate = zeros(sweepsP, sweepsQ);

fileTraining = 'rfdata2.xlsx';
fileScanning = 'rfdata1.xlsx';
f = xlsread(fileTraining, 'E6:IJ6');
XT = xlsread(fileTraining, 'E7:IJ1006');
XS = xlsread(fileScanning, 'E7:IJ1006');
Length = size(XT, 1);
channels = size(XT, 2);

%--------------------------------------------------------------------------------------------------
% Train DSA algorithm
%--------------------------------------------------------------------------------------------------
for p = linspace(startP, stopP, sweepsP)
    x = (-1)*startP + p + 1;
    
    % Converting sampled spectrum data into binary occupancy arrays
    trainer = thresholder(XT, p);
    M = thresholder(XS, p);
    occupied = sum(M, 2);
    vacant = Length - occupied;
    counts = occupancy(trainer);
    
    % Calculate survival/hazard function
    periodsIdle = sum(counts, 2);
    h = zeros(channels, Length);
    H = zeros(channels, Length);
    Ti = zeros(channels, Length);
    for r = 1:channels
        n(r) = length(find(counts(r, :)));
        Tl = [];
        
        for i = 1:Length
            Tl = [Tl, i*ones(1, counts(r, i))];
        end
        Ti(r, :) = [Tl, zeros(1, Length - size(Tl, 2))];

        j = 1;
        
        for t = 1:n(r)
            temp = 0;    
            while t >= Ti(r, j)
                temp = temp + 1/(periodsIdle(r) - j + 1);
                j = j + 1;
                if j > periodsIdle(r)
                    j = periodsIdle(r);
                    break
                end
            end
            h(r, t) = temp;
            if t == 1
                H(r, t) = temp;
                if n(r) == 1
                    H(r, t:Length) = temp;
                end
            elseif t == n(r)   
                H(r, t:Length) = H(r, t-1) + temp;
            else
                H(r, t) = H(r, t-1) + temp;
            end    
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
%         requests = [zeros(1, stopQ - q), ones(1, q)];
%         requests = repmat(requests, 1, Length/stopQ); 
%         requests = [ 0 1 ];
%         requests = repmat(requests, 1, Length/2);
%--------------------------------------------------------------------------
% Variant 2: Random secondary transmit request
%--------------------------------------------------------------------------
%         requests = zeros(1, Length);
%         for k = 1:Length
%             roll = stopQ * rand;
%             if q >= roll
%                 requests(k) = 1;
%             elseif q < roll
%                 requests(k) = 0;
%             end
%         end
%==========================================================================================
% Variant 2: Poisson distributed SU transmit request
%==========================================================================================
        requests = spectrum_occ_poiss(channels, Length, S1, S2);
%------------------------------------------------------------------------------------------
        schedule = zeros(channels, Length + 100);         % transmit grant schedule for secondary user
        transmit = zeros(channels, Length);
        interfere = zeros(channels, Length);
        
        for i = 1:channels
            t = 0;
            %--------------------------------------------------------------
            if vacant(i) == 0
               vacant(i) = 1;                   % TEMPORARY FIX !! 
            end
            if periodsIdle(i) == 0
               periodsIdle(i) = 1; 
            end
            %--------------------------------------------------------------
            for j = 1:Length
                sample = M(i, j);
                if sample == 0
                    t = t + 1;
                    if schedule(i, j) == 1
                        transmit(i, j) = 1;
                    end
                    if requests(i, j) == 1
                    %-------------------------------------------------------------
                    % Algorithm 1
                    %-------------------------------------------------------------
%                         T = t + tau;
%                         if T > Length
%                             T = Length; 
%                         end
%                         if H(i, T) - H(i, t) < theta
%                             schedule(i, (j + 1) : (j + tau)) = 1;
%                         end
                    %-------------------------------------------------------------
                    % Algorithm 2
                    %-------------------------------------------------------------
%                         tau = 1;
%                         while (H(i, t + tau)) < theta
%                             tau = tau + 1;
%                         end
%                         tau = tau - 1;
%                         schedule(i, (j + 1) : (j + 1 + tau)) = 1;
                    %------------------------------------------------------
                        tau = 1;
                        if (t + tau) > Length
                            tau = 1;
                        else
                            while (H(i, t + tau) - H(i, t)) < theta
                                tau = tau + 1;
                                if (t + tau) > Ti(i, periodsIdle(i))
                                   break
                                end
                            end
                        end
                        tau = tau - 1;
                        schedule(i, (j + 1) : (j + 1 + tau)) = 1;
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
        util(x, y) = 100*mean(transTot./vacant, 1);
        interfTot = sum(interfere, 2);
        interfRate(x, y) = 100*mean(interfTot./(Length), 1);
    end    
end
