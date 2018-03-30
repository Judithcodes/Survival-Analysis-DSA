% Survival analysis-based dynamic spectrum access algorithm
% Modified spectrum segment allocation method
% Channel mean idle time varies over length of sampling period
%
% Based on 2017 journal and conference paper by T.A. Hall et al.
%--------------------------------------------------------------------------

% Training algorithm with spectrum occupancy data representative of channel
% characteristics

% Simulation variables
Length = 1000000;                    % number of samples in each channel of spectrum occupancy data
channels = 1;                   % number of channels in test occupancy matrix
Tw = 5000;                      % window length
nWin = floor(Length/Tw);            % number of windows in sample array
nTrain = 10;                        % number of windows before retraining
t = 0;                            % time marker
% tau = 15;                          % transmit duration requested
threshold = 0.9;                  % interference threshold (probability of successful transmission)
theta = (-1)*log(threshold);

transTot = zeros(nWin, 1);
util = zeros(nWin, 1);
interfTot = zeros(nWin, 1);
interfRate = zeros(nWin, 1);
UTIL = zeros(50, 1);
INTERF = zeros(50, 1);

S1 = 15;            % SU request event rate
S2 = 15;            % SU idle event rate

for g = 1:50
    tau = g;
    
    % Occupancy data
    m = 1:nWin;
    %--------------------------------------------------------------------------
    % P1 = 50;                                % Occupancy event rate (lambda)
    % P2 = 20 + 15*sin(2*pi*m/100);            % Vacancy event rate 
    %--------------------------------------------------------------------------
    P2 = 2 + abs(100-m);
    P1(1:100) = 2 + m(1:100);
    P1(101:200) = 102 - abs(100 - m(101:200));
    %=============================================================================
    % Variant 1: Time-varying mean vacancy, exponential
    %=============================================================================
    % M = [];
    % for i = 1:nWin
    %    M =  [M, spectrum_occ_exp(channels, Tw, P2(i), P1(i))];
    % end
    %=============================================================================
    % Variant 2: Time-varying mean vacancy, dual poisson processes
    %=============================================================================
    M = [];
    for i = 1:nWin
       M =  [M, spectrum_occ_poiss(channels, Tw, P1(i), P2(i))];
    end
    %=============================================================================

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

    for m = 1:nWin
        start = (m-1)*Tw + 1;
        stop = m*Tw;
        window = M( : , start:stop);
        reqLocal = requests( : , start:stop);
        transmit = zeros(channels, Tw);         % segments where secondary user successfully transmits
        interfere = zeros(channels, Tw);        % segments where secondary user collides with primary user
        schedule = zeros(channels, Tw + 100);         % transmit grant schedule for secondary user

        occupied = sum(window);
        vacant = Tw - occupied;

        % Initial training on first window
        if m == 1
            %=====================================================================================================
            % Train DSA algorithm
            %=====================================================================================================
            % Generate array with number of instances of each length of idle period in
            % training vector
            counts = occupancy(window);
            n = length(find(counts));
            H = zeros(1, Tw);

            % Calculate survival/hazard function
            periodsIdle = sum(counts);
            pdf = counts./periodsIdle;
            cdf = cumsum(pdf);
            ccdf = cumsum(pdf, 'reverse');
            %=============================================================================
            % Initial Cumulative Hazard Function
            %=============================================================================
            Ti = [];
            for i = 1:Tw
                Ti = [Ti, i*ones(1, counts(i))];
            end

            h = zeros(1, Tw);
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
                    H(t:Tw) = H(t-1) + temp;
                else
                    H(t) = H(t-1) + temp;
                end    
            end
            %-----------------------------------------------------------------------------
        else        
            n = length(find(counts));

            %==================================================================
            % Retrain every nTrain windows
            %==================================================================
    %         if rem(m, nTrain) == 1
    %             % Calculate survival/hazard function
    %             periodsIdle = sum(counts);
    %             n = length(find(counts));
    %             H = zeros(1, Tw);
    %             
    %             % Periodic Rebuild of Cumulative Hazard Function
    %             %--------------------------------------------------------------
    %             Ti = [];
    %             for i = 1:Tw
    %                 Ti = [Ti, i*ones(1, counts(i))];
    %             end
    % 
    %             h = zeros(1, Tw);
    %             j = 1;
    %             for t = 1:n
    %                 temp = 0;    
    %                 while t >= Ti(j)
    %                     temp = temp + 1/(periodsIdle - j + 1);
    %                     j = j + 1;
    %                     if j > periodsIdle
    %                         j = periodsIdle;
    %                         break
    %                     end
    %                 end
    %                 h(t) = temp;
    %                 if t == 1
    %                     H(t) = temp;
    %                 elseif t == n   
    %                     H(t:Tw) = H(t-1) + temp;
    %                 else
    %                     H(t) = H(t-1) + temp;
    %                 end    
    %             end
    %         end
            %==================================================================
            % Retrain every window
            %==================================================================
    %         % Calculate survival/hazard function
    %         periodsIdle = sum(counts);
    %         n = length(find(counts));
    %         H = zeros(1, Tw);
    % 
    %         % Rebuild Cumulative Hazard Function
    %         %------------------------------------------------------------------
    %         Ti = [];
    %         for i = 1:Tw
    %             Ti = [Ti, i*ones(1, counts(i))];
    %         end
    %         h = zeros(1, Tw);
    %         j = 1;
    %         for t = 1:n
    %             temp = 0;    
    %             while t >= Ti(j)
    %                 temp = temp + 1/(periodsIdle - j + 1);
    %                 j = j + 1;
    %                 if j > periodsIdle
    %                     j = periodsIdle;
    %                     break
    %                 end
    %             end
    %             h(t) = temp;
    %             if t == 1
    %                 H(t) = temp;
    %             elseif t == n   
    %                 H(t:Tw) = H(t-1) + temp;
    %             else
    %                 H(t) = H(t-1) + temp;
    %             end    
    %         end
            %-----------------------------------------------------------------------------

            % Scan test matrix of occupancy data and grant or deny transmission
            % requests
            counts = zeros(channels, Tw);
            for i = 1:channels
                t = 0;
                for j = 1:Tw
                    sample = window(i, j);
                    if sample == 0
                        t = t + 1;
                        if ((j+1) > Tw) || (window(i, (j+1)) == 1)
                            counts(i, t) = counts(i, t) + 1;
                        end
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
                            %=============================================================
                            % Algorithm 2
                            %=============================================================
%                             tau = 1;
%                             while (H(t + tau) - H(t)) < theta
%                                 tau = tau + 1;
%                                 if (t + tau) > Ti(periodsIdle)
%                                    break
%                                 end
%                             end
%                             tau = tau - 1;
%                             schedule(i, (j + 1) : (j + 1 + tau)) = 1;
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
        end

        % Calculate metrics
        transTot(g, m) = sum(transmit, 2);
        util(g, m) = 100*transTot(g, m)./vacant;
        interfTot(g, m) = sum(interfere, 2);
        interfRate(g, m) = 100*interfTot(g, m)./(Tw);

    end

    UTIL(g) = mean(util(g, :));
    INTERF(g) = mean(interfRate(g, :));
    
end


