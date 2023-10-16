function [leak] = interactionMatrix(G,gw,wg,wc,bw)
% interactionMatrix.m
% 
% Compute the interaction matrix of a cascade graphic equalizer containing
% the leak factors to account for the band interaction when assigning
% filter gains. All filters are Orfanidis peak/notch filters with
% adjustable bandwidth gain.
% 
% Input parameters:
% G  = Linear gain at which the leakage is determined
% gw = Gain factor at bandwidth (0.5 refers to db(G)/2)
% wg = Command frequencies i.e. filter center frequencies (in rad/sample)
% wc = Design frequencies (rad/sample) at which leakage is computed
% bw = Bandwidth of filters in radians
% 
% Output:
% leak = N by M matrix showing how the magnitude responses of the
% band filters leak to the design frequencies. N is determined from the
% length of the array wc (number of design frequencies) whereas M is 
% determined from the length of wg (number of filters)
%
% Uses pareq.m
%
% Written by Vesa Valimaki, Espoo, Finland, 12 April 2016
% Modified by Juho Liski, Espoo, Finland, 26 September 2016
% Bug fix by Vesa Valimaki & Juho Liski, Espoo, Finland, 7 November 2018
%
% Aalto University, Dept. of Signal Processing and Acoustics

M = length(wg);  % The number of center frequencies (filters)
N = length(wc);  % The number of design frequencies
leak = zeros(M,N);  % Initialize interaction matrix
Gdb = db(G);  % Convert linear gain factor to dB
Gdbw = gw*Gdb;  % dB gain at bandwidth
Gw = 10.^(Gdbw/20);  % Convert to linear gain factors

% Estimate leak factors of peak/notch filters
if sum(abs(Gdb))  ~= 0
    for m = 1:M;  % Filters
        [num,den] = pareq(G(m), Gw(m), wg(m), bw(m)); % Parametric EQ filter
        H = freqz(num,den,wc);  % Evaluate response at wc frequencies
        Gain = db(H)/Gdb(m);  % Normalized interference (Re 1 dB)
        leak(m,:) = abs(Gain);  % Store gain in a matrix
    end
else
    leak = reshape([zeros(M);eye(M)],M,[]); % To avoid warning caused
    leak = leak(:,2:end);                   % by all 0-dB command gains
end
end
