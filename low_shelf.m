function [num, den] = low_shelf(wc,fs, GL, GH)
% implementation of first-order low-shelf filter
% K. Prawda, 5.10.2023

% INPUT:
% wc - crossover frequency (Hz)
% fs - sampling rate
% GL - gain in the low frequencies (linear)
% GH - gain in the high frequencies (linear)

% OUTPUT
% num - numerator coefficients
% den - denominator coefficients

wH =2*pi*wc/fs; % crossover frequency in [rad]
G = GL/GH;

    
%  filter coefficients
aH0 = tan(wH/2) + sqrt(G);
aH1 = tan(wH/2) - sqrt(G);
bH0 = G*tan(wH/2) + sqrt(G);
bH1 = G*tan(wH/2) - sqrt(G);

num = GH.*[bH0 bH1];
den = [aH0 aH1];

end