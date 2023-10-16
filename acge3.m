function [numsopt,densopt, Hopttot] = acge3(Gdb, flag, fs)
% acge3.m
% 
% Design third-octave EQ according to the method presented by J. Liski and
% V. Valimaki in "The Quest for the Best Graphic Equalizer," in Proc. 
% DAFx-17, Edinburgh, UK, Sep. 2017.
% 
% Input parameters:
% Gdb  = command gains in dB, size 31x1
% flag - use median (1) or not (0)
%  fs - sampling rate (Hz)
% 
% Output:
% numsopt = numerator parts of the 31 filters
% densopt = denominator parts of the 31 filters
%
% Uses pareq.m and interactionMatrix.m
%
% Created by Juho Liski and Vesa Valimaki, Otaniemi, Espoo, Finland, 21 June 2017
% Modified by Juho Liski, Otaniemi, Espoo, Finland, 6 May 2019
%
% %% Updated by Karolina Prawda, 13.10.2023
% 
% Aalto University, Dept. of Signal Processing and Acoustics

if flag ==1 
    G0 = median(Gdb);
else 
    G0 = 0;
end

gDB = Gdb;
Gdb = -G0 + Gdb;


fc1 = 10^3 * (2 .^ ([-17:13]/3)); % Log center frequencies for filters
fc2 = zeros(1,61); % Center frequencies and intermediate points between them
fc2(1:2:61) = fc1;
for k = 2:2:61
    fc2(k) = sqrt(fc2(k-1)*fc2(k+1));  % Extra points are at geometric mean frequencies
end
wg = 2*pi*fc1/fs;  % Command gain frequencies in radians
wc = 2*pi*fc2/fs;  % Center frequencies in radians for iterative design with extra points
gw = 0.4; % Gain factor at bandwidth (parameter c)
bw = 2*pi/fs*[9.178 11.56 14.57 18.36 23.13 29.14 36.71 46.25 58.28 73.43 ...
    92.51 116.6 146.9 185.0 233.1 293.7 370.0 466.2 587.4 740.1 932.4 ...
    1175 1480 1865 2350 2846 3502 4253 5038 5689 5570]; % EQ filter bandwidths

leak = interactionMatrix(10^(17/20)*ones(1,31),gw,wg,wc,bw); % Estimate leakage b/w bands
Gdb2 = zeros(61,1);
Gdb2(1:2:61) = Gdb;
for k = 2:2:61
    Gdb2(k) = (Gdb2(k-1)+Gdb2(k+1))/2; % Interpolate target gains linearly b/w command gains
end
Goptdb = leak'\Gdb2;      % Solve first estimate of dB gains based on leakage
Gopt = 10.^(Goptdb/20);    % Convert to linear gain factors

% Iterate once
leak2 = interactionMatrix(Gopt,gw,wg,wc,bw); % Use previous gains
G2optdb = leak2'\Gdb2;     % Solve optimal dB gains based on leakage
G2opt = 10.^(G2optdb/20);   % Convert to linear gain factors
G2woptdb = gw*G2optdb  ;    % Gain at bandwidth wg

G2wopt = 10.^(G2woptdb/20); % Convert to linear gain factor


% Design filters with optimized gains
numsopt = zeros(3,31);  % 3 num coefficients for each 10 filters
densopt = zeros(3,31);  % 3 den coefficients for each 10 filters
for k = 1:31,
    [num,den] = pareq(G2opt(k), G2wopt(k), wg(k), bw(k)); % Design filters
    numsopt(:,k) = num;
    densopt(:,k) = den;
    
end



numsopt = numsopt.*10^(G0/(20*length(fc1)));
% %%% Evaluation and plotting of the frequency response
Nfreq = 2^9;  % Number of frequency points for frequency response evaluation
% Log frequency points
w = [ logspace(log10(1),log10(fs/2-1),Nfreq-1), fs/2];  
% Evaluate frequency responses
Hopt = ones(Nfreq,31);   % Frequency response of individual filters
Hopttot = ones(Nfreq,1); % Frequency response of the cascade EQ
for k = 1:31
    Hopt(:,k) = freqz(numsopt(:,k),densopt(:,k),w,fs);
    Hopttot = Hopt(:,k) .* Hopttot;
end


