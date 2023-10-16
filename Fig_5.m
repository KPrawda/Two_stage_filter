% script to analyze the data from the Arni measurements

% housekeeping
clear variables
close all
clc
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
%%
fs = 44100;

freqRange = [20, 20000];
bpo1 = 1; % full octave
bpo3 = 3; % third octave

% numB3 = 10; % 10 octave bands
numB3 = 30; % 10 octave bands
bandInd = 1:numB3;


%% load RT from Arni
load('two-stage-RT-values.mat')

%%
colors = [240, 149, 161; 201, 109, 121; 161, 82, 92]./255;
colors2 = [0, 117, 196; 161, 205, 244; 126, 168, 190]./255; 
labs = string([0 10 30, 100 300 1000 3000 10000, 20000]./1000);

%% statistical analysis
nTrial = 1000;%numComb;
Nfreq = 2^9; 
w = [0, logspace(log10(1),log10(fs/2-1),Nfreq-2), fs/2];
rng(0)
dls = round(0.3*rand(1, nTrial)*fs);
dls(dls < 0.01*fs) = 0.01*fs;

%% initialize variables
t60_target =zeros(Nfreq, nTrial);

t60_SHE = zeros(Nfreq, nTrial);
t60_GEQ = zeros(Nfreq, nTrial);
t60_MED = zeros(Nfreq, nTrial);
t60_NOT = zeros(Nfreq, nTrial);


RTerror_SHE = zeros(Nfreq, nTrial);
RTerror_SHE_proc = zeros(Nfreq, nTrial);
RTerror_GEQ = zeros(Nfreq, nTrial);
RTerror_GEQ_proc = zeros(Nfreq, nTrial);
RTerror_MED = zeros(Nfreq, nTrial);
RTerror_MED_proc = zeros(Nfreq, nTrial);
RTerror_NOT = zeros(Nfreq, nTrial);
RTerror_NOT_proc = zeros(Nfreq, nTrial);

%% run on Arni dataset

for it  = 1: nTrial
    %% use only GEQ 
    method = 'geq';
    % get the frequency response of the attenuation filter
    [HGEQ] = twoFilters(rt_(:, it), dls(it), fs, method);
    
    % convert to reverberation time values
    t60_GEQ(:, it) = -60*dls(it)./(fs*20*log10(abs(HGEQ)));
    
    %% use median gain and GEQ 
    method = 'median';
    % get the frequency response of the attenuation filter
    [HMED] = twoFilters(rt_(:, it), dls(it), fs, method);
    
    % convert to reverberation time values
    t60_MED(:, it) = -60*dls(it)./(fs*20*log10(abs(HMED)));
    
    %% use median gain, notch filter, and GEQ 
    method = 'notch';
    % get the frequency response of the attenuation filter
    [HNOT] = twoFilters(rt_(:, it), dls(it), fs, method);
    
    % convert to reverberation time values
    t60_NOT(:, it)  = -60*dls(it)./(fs*20*log10(abs(HNOT)));
    
    %% use shelf filter and GEQ 
    method = 'shelf';
    % get the frequency response of the attenuation filter
    [HSHE, w, target_mag] = twoFilters(rt_(:, it), dls(it), fs, method,1000);
    
    % convert to reverberation time values
    t60_SHE(:, it)  = -60*dls(it)./(fs*20*log10(abs(HSHE)));
    
    %% target RT
    t60_target(:, it) = -60*dls(it)./(fs*target_mag);


    %% RT error
    % get the RT error, linear and %
    RTerror_SHE(:,it) = t60_target(:,it)-t60_SHE(:,it); 
    RTerror_SHE_proc(:,it) = 100*RTerror_SHE(:,it)./t60_target(:,it);
       
    RTerror_GEQ(:,it) = t60_target(:,it)-t60_GEQ(:,it);
    RTerror_GEQ_proc(:,it) = 100*RTerror_GEQ(:,it)./t60_target(:,it);
        
    RTerror_MED(:,it) = t60_target(:,it)-t60_MED(:,it);
    RTerror_MED_proc(:,it) = 100*RTerror_MED(:,it)./t60_target(:,it);

    RTerror_NOT(:,it) = t60_target(:,it)-t60_NOT(:,it);
    RTerror_NOT_proc(:,it) = 100*RTerror_NOT(:,it)./t60_target(:,it);
end


%%
lw = 3;
f = figure(3); clf; hold on
edges  = -1000:1000;
binCenters = (edges(1:end-1) + edges(2:end))/2;
[RT_0] = histcounts((RTerror_GEQ_proc), edges, 'Normalization','probability');
[RT_1] = histcounts((RTerror_MED_proc), edges, 'Normalization','probability');
[RT_sh] = histcounts((RTerror_NOT_proc), edges, 'Normalization','probability');
[RT_s] = histcounts((RTerror_SHE_proc), edges, 'Normalization','probability');
RT_0(RT_0==0) = 10^-9;
RT_1(RT_1==0) = 10^-9;
RT_sh(RT_sh==0) = 10^-9;
RT_s(RT_s==0) = 10^-9;

xlim([-110 110])
ylim([10^-6 1])
set(gca, 'yscale', 'log')
ylabel('Probability', 'interpreter', 'latex')
xlabel('$T_{60}$ error ($\%$)', 'interpreter', 'latex')
set(gca, 'yTick',[0.0001   0.01 1], 'yticklabels',[0.0001  0.01 1], 'xtick', -100:25:100, 'Fontsize',12)

plot(binCenters, (RT_0), 'color',colors(1, :), 'linewidth', lw)
plot(binCenters, (RT_1), 'color',colors2(1, :), 'linewidth', lw)
plot(binCenters, (RT_sh),'-.', 'color',colors2(2, :), 'linewidth', lw*0.5)
plot(binCenters, (RT_s), 'r', 'linewidth', lw)

box on

legend('GEQ only', 'Median gain','Notch filter','Proposed', 'Location','northwest', 'interpreter', 'latex', 'numcolumns',1)
f.Position(4)=250
pos = f.Position;
