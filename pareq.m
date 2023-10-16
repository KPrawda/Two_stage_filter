function [num, den] = pareq(G, GB, w0, B)
% [num, den] = pareq(G, GB, w0, B)
% Second-order parametric equalizing filter design 
% with adjustable bandwidth gain 
%
% [num, den] = pareq(G, GB, w0, B)
%
% Parameters
% G = peak gain (linear)
% GB = bandwidth gain (linear)
% w0 = center frequency (rads/sample)
% B = bandwidth (rads/sample)
%
% Output
% num = [b0, b1, b2] = numerator coefficients
% den = [1,  a1, a2] = denominator coefficients
%
% Written by Vesa Valimaki, August 24, 2016
%
% Ref. S. Orfanidis, Introduction to Signal Processing, p. 594
% We have set the dc gain to G0 = 1.

if G==1,
    beta = tan(B/2);  % To avoid division by zero, when G=1
else
    beta = sqrt(abs(GB^2 - 1)/abs(G^2 - GB^2)) * tan(B/2);
end
num = [(1 + G*beta), -2*cos(w0), (1 - G*beta)] / (1+beta);
den = [1, -2*cos(w0)/(1+beta), (1-beta)/(1+beta)];
