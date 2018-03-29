function Band = spectrum_occ_exp( channels, samples, L0, L1 )
%Spectrum Occupancy, Exponential Distribution
%   Returns matrix of randomly generated binary values simulating spectrum
%   occupancy.
%       channels = # of rows
%       samples = minimum # of columns

% Number of columns in channel matrix
n = 6*samples;

% Generating single channel occupancy data
T = [];
toggle = 0;
for i = 1:samples
    duration = 0;
    if toggle == 0                  % Generates vacant period
        while duration == 0
            duration = round(exprnd(L0));
        end
        T = [T , zeros(1, duration)];
        toggle = ~toggle;
    elseif toggle == 1              % Generates occupied period
        while duration == 0
            duration = round(exprnd(L1));
        end
        T = [T , ones(1, duration)];
        toggle = ~toggle;
    end
end

% Generating band occupancy data
if channels < 2
    Band = T( 1 , 1:samples );
else    
    G = zeros(channels, n);
    G( 1 , : ) = [ T , zeros(1, n-size(T, 2)) ]; 
    for i = 2:channels
        T = [];
        toggle = 0;
        for j = 1:samples
            if toggle == 0
                duration = round(exprnd(L0));
                T = [T , zeros(1, duration)];
                toggle = ~toggle;
            elseif toggle == 1
                duration = round(exprnd(L1));
                T = [T , ones(1, duration)];
                toggle = ~toggle;
            end
        end
        G( i , : ) = [ T , zeros(1, n - size(T,2)) ];
    end
    Band = G( 1:channels , 1:samples );
end

end

