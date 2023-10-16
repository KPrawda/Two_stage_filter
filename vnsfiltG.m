function [y,ir,vns,bir] = vnsfiltG(x,lens,RT,dur,d, noseg, fs, method)
%VNSFILT3
% A Velvet-Noise Sequence (VNS) reverberator structure
% of order M (number of parallel paths, defined by lens) is constructed 
% and the input signal is filtered with it. 
% 
% Uses velvetnsmkx & twoFilters and nested functions
%
% Input parameters: 
% x = Input signal (monophonic)
% lens = Lengths of the EVN sequences 
%       (preferably consecutive primes times 20*M)
% g = decay parameter. With g = 1, the reverberation time 
%     at 1 kHz is 1 sec. 
% dur = Duration of the impulse response in seconds (e.g., dur = 2 sec)
%
% Output parameters: 
% y = Output signal (e.g., reverberated signal)
% ir = Impulse response of the VNS reverberator
%      (Note that pseudo-random numbers are used, 
%       so IR will be slightly different every time.)
% vns = velvet noise sequences
% 
% The sampling rate is 44100 kHz.
%
% Created by Vesa Valimaki, Otaniemi, Espoo, 9.3.2014
% Modified by Vesa Valimaki, Otaniemi, Espoo, 20.11.2018


len = round(dur*fs);    % Length is samples
density = 2205;         % A sufficient overall density (2000 impulses/sec) is used

M = length(lens);       % Number of branches
Td =80;% M*fs/density;      % This is the grid size of each branch VN sequence
                        % i.e., average distance of impulses in each sequence
% delta = floor(density/(M+1))/(density);  % Defines the range where impulses will occur 
delta = round(floor(density/M)/(density),2);  % Defines the range where impulses will occur 
                        % (e.g. delta = 0.33, when M = 3)
                        % The FLOOR is used to avoid collisions at the boundaries
N = round(fs/density);  % Delay between branches (for interleaving) is fs/density
                        % which is the same as Td, the grid size of overall system

% Generate M Velvet-Noise sequences (lengths defined by LENS)
vns = zeros(M,max(lens));  % Velvet-Noise sequence matrix (all separately)
% d = 2; %delay parameter - how many multiplications of M*Td will the signal be delayed by
g = -60./RT;
wc = 10000;
for m = 1:M
    RTloss = (RT-(m-1)*Td*d/fs);
    decdB = -60./RTloss ;  % dacay rate in dB/s
    gains = 10.^(decdB/20);  % frequency dependent gain
    gdec(m,:) = gains.^(lens(m)/fs);
    gdB(m,:) = db(gdec(m,:));
    [~,~,~, htot(m,:)] = twoFilters(RTloss, lens(m), fs, method, wc);
    
end

for m = 1:M,
    rng(3)
    [xxx,sm,nnonzero(m)] = velvetnsmkx(lens(m), density/M, delta);
    size(xxx);
    vns(m,1:lens(m)) = sm;
end


cir = zeros(M,len + d*M*(M-1)*Td);  % Reserve 2 sec for IR of each comb filter
bir = zeros(M,len + d*M*(M-1)*Td);  % Reserve 2 sec for IR of each branch (with VNS)



dc_vns = zeros(M,max(lens) );
for m = 1:M % first: delay for making the "steps" in the decay less visible
        if d ==0
            cir(m,(m-1)*Td*d+1:((m-1)*Td*d)+len) = impz( 1,[1 zeros(1,(lens(m)-2)) htot(m,:)],len);
        
        else
            cir(m,(m-1)*Td*d+1:((m-1)*Td*d)+len) = impz( 1,[1 zeros(1,(lens(m)-2)) htot(m,:)],len);

        end

            if nargin >5 && noseg ~= 0
                if noseg > 0
                    d1 = cir(m,(m-1)*Td*d+1);
                    d2 = cir(m,(m-1)*Td*d+lens(m));
                    dc_sc = [1, (d1-(d1-abs(d2))/3)/d1, (d1-2*(d1-abs(d2))/3)/d1];
                    dc_len = [0 round(lens(m)*0.25), round(lens(m)*0.35), round(lens(m)*0.40)];

                    for i = 1:3
    
                        dc_vns(m,sum(dc_len(1:i))+1: sum(dc_len(1:i+1))) = dc_sc(i)*vns(m,sum(dc_len(1:i))+1: sum(dc_len(1:i+1)));  
                    end
                    bir(m,:) = filter(dc_vns(m,1:lens(m)),1,cir(m,:));
                elseif noseg < 0
                    sum(cir(m,(m-1)*Td*d+lens(m)+1:m*Td*d+lens(m)));
                    d1 = 20*log10(abs(sum(cir(m,(m-1)*Td*d+1:m*Td*d))));
                    d2 = 20*log10(abs(sum(cir(m,(m-1)*Td*d+lens(m)+1:m*Td*d+lens(m)))));
                    dLb = d1-d2;
                    alpha = -log(10^(-dLb/20))./nnonzero(m);
                    nnonzero(m)
                    for n = 1: nnonzero(m)
                    dc_vns_(m,(n-1)*Td+1:n*Td) = vns(m,(n-1)*Td+1:n*Td).*exp(-alpha*n);
                    end
                    bir(m,:) = filter(dc_vns_(m,1:lens(m)),1,cir(m,:));
%                      vns = dc_vns_;
                end
            else
                bir(m,:) = filter(vns(m,1:lens(m)),1,cir(m,:)); % Conv with VNS - SFIR filtering
            end
            xxx = [zeros(1,(m-1)*N) bir(m,:)];              % Delay for interleaving z^-D
            bir(m,:) = xxx(1:len + d*M*(M-1)*Td);               % Truncate to DUR*FS samples
end
ir =cir;

y = sum(bir);


 

