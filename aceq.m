function [numsopt,densopt, Hopttot] = aceq(Gdb, flag, fs)
% aceq.m
% 
% Design octave EQ according to the new method presented by V. Valimaki
% and J. Liski in "Accurate Cascade Graphic Equalizer," IEEE Signal 
% Processing Letters, 2017

% 
% Input parameters:
% Gdb  = command gains in dB, size 10x1
% flag - use median (1) or not (0)
%  fs - sampling rate (Hz)
% 
% Output:
% numsopt = numerator parts of the 10 filters
% densopt = denominator parts of the 10 filters
%
% Uses pareq.m and interactionMatrix.m
%
% Created by Juho Liski and Vesa Valimaki, Otaniemi, Espoo, Finland, 29 September 2016
% Modified by Juho Liski, Otaniemi, Espoo, Finland, 23 December 2016
% 
% Updated by Karolina Prawda, 13.10.2023
%
% Aalto University, Dept. of Signal Processing and Acoustics
% 
if flag ==1 
    G0 = median(Gdb);
else 
    G0 = 0;
end

gDB = Gdb;
Gdb = -G0 + Gdb;
numB = length(Gdb);
if numB == 10
    numF = 19;
    fc1 = 16000./(2.^(9:-1:0));

end
% fs  = 44100;%44.1e3;  % Sample rate
% fc1 = 16000./(2.^(9:-1:0)); % Exact log center frequencies for filters
fc2 = zeros(1,numF); % Center frequencies and intermediate points between them
fc2(1:2:numF) = fc1;
for k = 2:2:numF
    fc2(k) = sqrt(fc2(k-1)*fc2(k+1));  % Extra points are at geometric mean frequencies
end
wg = 2*pi*fc1/fs;  % Command gain frequencies in radians
wc = 2*pi*fc2/fs;  % Center frequencies in radians for iterative design with extra points
gw = 0.3; % Gain factor at bandwidth (parameter c)
bw = zeros(1,numB); % EQ filter bandwidths
for m = 1:numB
    bw(m) = 1.5*wg(m); 
end
if numB == 10, bw(8) = 0.93*bw(8); bw(9) = 0.78*bw(9); bw(10) = 0.76*wg(10); end% Additional adjustmenst due to asymmetry

leak = interactionMatrix(10^(numF/20)*ones(1,numB),gw,wg,wc,bw); % Estimate leakage b/w bands
Gdb2 = zeros(numF,1);
Gdb2(1:2:numF) = Gdb;

gDB2 = zeros(numF,1);
gDB2(1:2:numF) = gDB;

for k = 2:2:numF
    Gdb2(k) = (Gdb2(k-1)+Gdb2(k+1))/2; % Interpolate target gains linearly b/w command gains
    gDB2(k) = (gDB2(k-1)+gDB2(k+1))/2;
end

Goptdb = leak'\Gdb2;      % Solve first estimate of dB gains based on leakage
Gopt = 10.^(Goptdb/20) ;   % Convert to linear gain factors

% Iterate once
leak2 = interactionMatrix(Gopt,gw,wg,wc,bw); % Use previous gains
G2optdb = leak2'\Gdb2;     % Solve optimal dB gains based on leakage
G2opt = 10.^(G2optdb/20);   % Convert to linear gain factors
G2woptdb = gw*G2optdb;      % Gain at bandwidth wg
G2wopt = 10.^(G2woptdb/20); % Convert to linear gain factor

% Design filters with optimized gains
numsopt = zeros(3,numB+1);  % 3 num coefficients for each 10 filters
densopt = zeros(3,numB+1);  % 3 den coefficients for each 10 filters
for k = 1:numB,
    [num,den] = pareq(G2opt(k), G2wopt(k), wg(k), bw(k)); % Design filters
    numsopt(:,k) = num;
    densopt(:,k) = den;
end

numsopt = numsopt*(10^(G0/(20*numB)));

%%% Evaluation and plotting of the frequency response
Nfreq = 2^9;  % Number of frequency points for frequency response evaluation
w = logspace(log10(1),log10(fs/2),Nfreq);  % Log frequency points

%%
% Evaluate frequency responses
Hopt = ones(Nfreq,numB);  % Frequency response of individual filters
Hopttot = ones(Nfreq,1);  % Frequency response of the cascade EQ
for k = 1:numB
    Hopt(:,k) = freqz(numsopt(:,k),densopt(:,k),w,fs);
    Hopttot = Hopt(:,k) .* Hopttot;
end

