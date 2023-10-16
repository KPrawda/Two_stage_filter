function [num, den] = high_shelf(GH, fs)
% implementation of first-order high-shelf filter
% K. Prawda, 5.10.2023

% INPUT:
% GH - gain in the high frequencies (linear)
% fs - sampling rate

% OUTPUT
% num - numerator coefficients
% den - denominator coefficients

wc_ = 20200/fs; % crossover frequency

wH =2*pi*wc_; % frequency in [rad]
    
%  coefficients
aH0 = sqrt(GH)*tan(wH/2) + 1;
aH1 = sqrt(GH)*tan(wH/2) - 1;
bH0 = sqrt(GH)*tan(wH/2) + GH;
bH1 = sqrt(GH)*tan(wH/2) - GH;

num = [bH0 bH1];
den = [aH0 aH1];
end