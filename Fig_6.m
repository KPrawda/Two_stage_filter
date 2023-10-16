%% synthesize RIRs with Interleaved Velvet Noise and different filter designs
% K. Prawda, 13.10.2023
% uses vnsfiltG and nested functions
%% IVN examples
% housekeeping
clear variables
close all
clc
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
%% load the sound example and the RT values
filename = 'pori.wav';
[rir, fs] = audioread(filename);
rir(1:1699) = []; % delete the zeros from before the direct sound
rir = rir./max(abs(rir)); % normalize 

load('RT_pori.mat');
%% synthesize reverb with IVN
% for reference, see: 
% V. Välimäki and K. Prawda. “Late-Reverberation Synthesis Using Interleaved Velvet-Noise Sequences”. In: IEEE/ACM Trans.
% Audio, Speech, and Lang. Process. 29 (Feb. 2021), pp. 1149–1160. doi: 10.1109/TASLP.2021.3060165.

lens = 20*4*[79 83 89 97];  % lenghts of EVN sequences (samples)
dur = 5;                    % duration of the final signal (s)
t = 0:1/fs:dur-1/fs;        % time vector 
d = 1;                      % delay parameter for the EVN sequences
noseg = 3;                  % number of segments in one loop iteration
x = [1, zeros(1, 100)];     % input is an impulse


% create signals using different filters 
y_shelf = vnsfiltG(x,lens,RT_,dur,d, noseg, fs, 'shelf'); % proposed
y_geq = vnsfiltG(x,lens,RT_,dur,d, noseg, fs,'GEQ');
y_median = vnsfiltG(x,lens,RT_,dur,d, noseg, fs,'median');
y_notch = vnsfiltG(x,lens,RT_,dur,d, noseg, fs,'notch');

%% early reflections, fadein and fadeout of the impulse response
earlyRef = rir(1:0.12*fs);  % early refl 120 ms
cf_length = 0.05*fs;        % fadein and fadeout length
fadeout = linspace(1,0, cf_length);
fadeout = [ones(1, length(earlyRef)-cf_length), fadeout];
fadein = linspace(0,1, cf_length);
fadein = [fadein, ones(1, length(y_shelf)-round(0.07*fs) - cf_length+1)];


ER_fadeout = earlyRef.*(fadeout');
REV_fadein_shelf = [zeros(1,round(0.07*fs)-1), y_shelf(round(0.07*fs):end).*fadein];
REV_fadein_geq = [zeros(1,round(0.07*fs)-1), y_geq(round(0.07*fs):end).*fadein];
REV_fadein_median = [zeros(1,round(0.07*fs)-1), y_median(round(0.07*fs):end).*fadein];
REV_fadein_notch = [zeros(1,round(0.07*fs)-1), y_notch(round(0.07*fs):end).*fadein];

ER_fadeout= [ER_fadeout; zeros(length(REV_fadein_shelf) - length(ER_fadeout), 1)];


% calculate signal powers for power adjustment between the early and late
% part
P_rir = sum(abs(rir(0.07*fs:0.12*fs)).^2);
P_er = sum(abs(ER_fadeout(0.07*fs:0.12*fs)).^2);

P_rev_shelf = sum(abs(REV_fadein_shelf(0.07*fs:0.12*fs)).^2);
P_rev_geq = sum(abs(REV_fadein_geq(0.07*fs:0.12*fs)).^2);
P_rev_median = sum(abs(REV_fadein_median(0.07*fs:0.12*fs)).^2);
P_rev_notch = sum(abs(REV_fadein_notch(0.07*fs:0.12*fs)).^2);

gain_shelf = (sqrt(P_rir) - sqrt(P_er)).^2/P_rev_shelf;
gain_geq = (sqrt(P_rir) - sqrt(P_er)).^2/P_rev_geq;
gain_median = (sqrt(P_rir) - sqrt(P_er)).^2/P_rev_median;
gain_notch = (sqrt(P_rir) - sqrt(P_er)).^2/P_rev_notch;

% fade in and fade out for the created reverbs
ivn_rir = ER_fadeout' + 2*gain_shelf*REV_fadein_shelf;
ivn_rir_geq = ER_fadeout' + 2*gain_geq*REV_fadein_geq;
ivn_rir_median = ER_fadeout' + 2*gain_median*REV_fadein_median;
ivn_rir_notch = ER_fadeout' + 2*gain_notch*REV_fadein_notch;
%% calculate spectrograms
NFFT = 2.^9;

overlap = 512;
window_length = 2*overlap;
FF = logspace(log10(1),log10(fs/2),NFFT);

[out_sh,f_sh,tt_sh] = spectrogram(ivn_rir,window_length, overlap, FF, fs, 'yaxis');
psn_sh = db(out_sh);
psn_sh = psn_sh-max(max(psn_sh));



[out_0,f_0,tt_0] = spectrogram(ivn_rir_geq,window_length, overlap, FF, fs, 'yaxis');
psn_0 = db(out_0);
psn_0 = psn_0-max(max(psn_0));



[outR,fR,ttR, psR] = spectrogram(rir,window_length, overlap, FF, fs, 'yaxis');
psnR = db(outR);
psnR = psnR-max(max(psnR));
[~, locMR] = max(max(psnR));
psnR = psnR(:, locMR:end);
ttR = ttR(1:end-locMR+1);
%% plot figure 6 from the paper
f = figure(2);clf

subplot(3,1,1) % reference RIR

pcol =  pcolor(ttR,fR,psnR);
pcol.EdgeColor = 'none';
set(gca, 'YDir', 'normal','Yscale', 'log')
colormap(flipud(gray))
caxis([-65 0])
ylim([5 fs/2+6000])
xlim([0 3])

clb = colorbar;
set(get(clb,'Label'),'String','Energy (dB)', 'Fontsize',12,'interpreter', 'latex')
clb.TickLabelInterpreter = 'latex';
clb.Ticks = [-60 -30 0];


xlabel('(a)', 'interpreter', 'latex')
set(gca,'Fontsize',12,'YTick', [ 10 100  1000  10000], 'YTicklabel',{ '0.01', '0.1'  '1'  '10'})
box on
set(gca, 'layer','top')



subplot(3,1,2) % RIR synthesized with IVN and GEQ only

pcol =  pcolor(tt_0,f_0,psn_0);
pcol.EdgeColor = 'none';
set(gca, 'YDir', 'normal','Yscale', 'log')
colormap(flipud(gray))
caxis([-65 0])
ylim([5 fs/2+6000])
xlim([0 3])

clb = colorbar;
set(get(clb,'Label'),'String','Energy (dB)', 'Fontsize',12,'interpreter', 'latex')
clb.TickLabelInterpreter = 'latex';
clb.Ticks = [-60 -30 0];

ylabel('Frequency (kHz)', 'interpreter', 'latex')
hold on

xlabel('(b)', 'interpreter', 'latex')
set(gca,'Fontsize',12,'YTick', [ 10 100  1000  10000], 'YTicklabel',{ '0.01', '0.1'  '1'  '10'})
box on
set(gca, 'layer','top')



subplot(3,1,3) % RIR synthesized with IVN and Proposed method 

pcol =  pcolor(tt_sh,f_sh,psn_sh);
pcol.EdgeColor = 'none';
set(gca, 'YDir', 'normal','Yscale', 'log')
colormap(flipud(gray))
caxis([-65 0])
ylim([5 fs/2+6000])
xlim([0 3])

clb = colorbar;
set(get(clb,'Label'),'String','Energy (dB)', 'Fontsize',12,'interpreter', 'latex')
clb.TickLabelInterpreter = 'latex';
clb.Ticks = [-60 -30 0];




xlabel({'Time (s)'; '(c)'}, 'interpreter', 'latex'); %ylabel('Frequency (kHz)', 'interpreter', 'latex')
set(gca,'Fontsize',12,'YTick', [ 10 100  1000  10000], 'YTicklabel',{ '0.01', '0.1'  '1'  '10'})
box on
set(gca, 'layer','top')


f.Position(end) =550;
