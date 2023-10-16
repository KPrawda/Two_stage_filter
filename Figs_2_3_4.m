%% implement two stage loss filter for reverberation synthesis
% K. Prawda, 5.10.2023
% uses twoFilters and nested functions
%% housekeeping
clear all 
close all
clc


%% load RT from Arni
load('two-stage-RT-values.mat')

%% setup some params
fs = 44100;
freqRange = [20, 20000];
bpo1 = 1; % full octave
bpo3 = 3; % third octave

numB1 = 10; % 10 octave bands
numB3 = 31; % 10 octave bands

%% running example of a reverbeartion time and target attenuation for FIgs 2--5
RT = rt_(:, 1); 
dL = 0.1*fs; %   delay line length


%% use only GEQ 
method = 'geq';
% get the frequency response of the attenuation filter
[HGEQ] = twoFilters(RT, dL, fs, method);

% convert to reverberation time values
t60_GEQ = -60*dL./(fs*20*log10(abs(HGEQ)));

%% use median gain and GEQ 
method = 'median';
% get the frequency response of the attenuation filter
[HMED] = twoFilters(RT, dL, fs, method);

% convert to reverberation time values
t60_MED = -60*dL./(fs*20*log10(abs(HMED)));

%% use median gain, notch filter, and GEQ 
method = 'notch';
% get the frequency response of the attenuation filter
[HNOT] = twoFilters(RT, dL, fs, method);

% convert to reverberation time values
t60_NOT  = -60*dL./(fs*20*log10(abs(HNOT)));

%% use shelf filter and GEQ 
method = 'shelf';
% get the frequency response of the attenuation filter
[HSHE, w, target_mag] = twoFilters(RT, dL, fs, method,300);

% convert to reverberation time values
t60_SHE  = -60*dL./(fs*20*log10(abs(HSHE)));

%% target RT
t60 = -60*dL./(fs*target_mag);

%%
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
%%
pinks = [240, 149, 161; 201, 109, 121; 161, 82, 92]./255;
blues = [0, 117, 196; 161, 205, 244]./255;
labs = ["1", "30", "300", "3k", "10k", "20k"];
%2, 103, 193 green blue
%% plot magnitude for GEQ only
f = figure(1); clf

%% subplot(2,2,1)
[ax11, ax12] = subplotResponse(w, target_mag, 20*log10(abs(HGEQ)), [3055, 2900, 3000, fs/2],  [-20.5 3],  pinks(1, :), [2,4,1]);

yline(0, ':', 'color', pinks(1,:), 'LineWidth', 1.5, 'Parent',ax11)

yline(0, ':', 'color', pinks(1,:), 'LineWidth', 1.5, 'Parent',ax12)

set(ax11, 'XTick',[1 30, 300 ],'xticklabels', labs(1:3), 'YTick',[-30 :10 :0], 'Fontsize',12, 'XTickLabelRotation',45)

ax11.YLabel.String = 'Magnitude (dB)';
ax11.YLabel.Interpreter = 'latex';

ax11.XLabel.String = '(a)';
ax11.XLabel.Interpreter = 'latex';

ax11.XLabel.Position(1) = ax11.XLabel.Position(1) +1800;
ax11.XLabel.Position(2) = ax11.XLabel.Position(2) -0.3;

ax11.Title.String = 'GEQ only';
ax11.Title.Interpreter = 'latex';
ax11.Title.Position(1) = 1800;


set(ax12, 'XTick',[3000 10000, 20000],'xticklabels', labs(end-2:end), 'YTick',[], 'Fontsize',12, 'YColor', 'none', 'XTickLabelRotation',45)


%% subplot(2,2,2)
[ax21, ax22] = subplotResponse(w, target_mag, 20*log10(abs(HMED)), [3055, 2900, 3000, fs/2], [-20.5 3],  blues(2, :), [2,4,3]);

yline(20*log10(abs(HMED(1))), ':', 'color', blues(2,:), 'LineWidth', 1.5, 'Parent',ax21)

yline(20*log10(abs(HMED(1))), ':', 'color', blues(2,:), 'LineWidth', 1.5, 'Parent',ax22)


set(ax21, 'XTick',[1 30, 300 ],'xticklabels', labs(1:3), 'YTick',[-30 :10 :0], 'Fontsize',12, 'XTickLabelRotation',45)


ax21.XLabel.String = '(b)';
ax21.XLabel.Interpreter = 'latex';

ax21.XLabel.Position(1) = ax21.XLabel.Position(1) +1800;
ax21.XLabel.Position(2) = ax21.XLabel.Position(2) -0.3;

ax21.Title.String = 'Median gain';
ax21.Title.Interpreter = 'latex';
ax21.Title.Position(1) = 1800;


set(ax22, 'XTick',[3000 10000, 20000],'xticklabels', labs(end-2:end), 'YTick',[], 'Fontsize',12, 'YColor', 'none', 'XTickLabelRotation',45)



%% subplot(2,2,3)
[ax31, ax32] = subplotResponse(w, target_mag, 20*log10(abs(HNOT)), [3055, 2900, 3000, fs/2], [-20.5 3],  blues(1, :), [2,4,5]);

%
GH =10^((target_mag(end) - median(target_mag))/20); % shelf filter gain
[Hnum, Hden] = high_shelf(GH, fs);
Hshelf_h = freqz(Hnum, Hden, w, fs);
%
semilogx(w,20*log10(abs(HNOT(1)))+20*log10(abs(Hshelf_h)),':', 'color', blues(1, :), 'LineWidth',1.5, 'Parent',ax31)
semilogx(w,20*log10(abs(HNOT(1)))+20*log10(abs(Hshelf_h)),':', 'color', blues(1, :), 'LineWidth',1.5, 'Parent',ax32)

set(ax31, 'XTick',[1 30, 300 ],'xticklabels', labs(1:3), 'YTick',[-30 :10 :0], 'Fontsize',12, 'XTickLabelRotation',45)

ax31.YLabel.String = 'Magnitude (dB)';
ax31.YLabel.Interpreter = 'latex';

ax31.XLabel.String = {'Frequency (Hz)', '(c)'};
ax31.XLabel.Interpreter = 'latex';

ax31.XLabel.Position(1) = ax31.XLabel.Position(1) +1800;
ax31.XLabel.Position(2) = ax31.XLabel.Position(2) -0.3;

ax31.Title.String = 'Notch filter';
ax31.Title.Interpreter = 'latex';
ax31.Title.Position(1) = 1800;


set(ax32, 'XTick',[3000 10000, 20000],'xticklabels', labs(end-2:end), 'YTick',[], 'Fontsize',12, 'YColor', 'none', 'XTickLabelRotation',45)

%% subplot(2,2,4)
[ax41, ax42] = subplotResponse(w, target_mag, 20*log10(abs(HSHE)), [3055, 2900, 3000, fs/2], [-20.5 3],  [1 0 0], [2,4,7]);

%
GL = max(10.^(target_mag./20)); % the gain for low frequencies
GH =10.^(target_mag(end)./20); % the gain for high frequencies

[Hnum, Hden] = low_shelf(300,fs, GL, GH);

% evaluate filter's frequency response
Hshelf_L = freqz(Hnum, Hden, w, fs);
%

semilogx(w,20*log10(abs(Hshelf_L)),'r:', 'LineWidth',1.5, 'Parent',ax41)
semilogx(w,20*log10(abs(Hshelf_L)),'r:', 'LineWidth',1.5, 'Parent',ax42)

set(ax41, 'XTick',[1 30, 300 ],'xticklabels', labs(1:3), 'YTick',[-30 :10 :0], 'Fontsize',12, 'XTickLabelRotation',45)

ax41.XLabel.String = {'Frequency (Hz)', '(d)'};
ax41.XLabel.Interpreter = 'latex';

ax41.XLabel.Position(1) = ax41.XLabel.Position(1) +1800;
ax41.XLabel.Position(2) = ax41.XLabel.Position(2) -0.3;

ax41.Title.String = 'Proposed';
ax41.Title.Interpreter = 'latex';
ax41.Title.Position(1) = 1800;


set(ax42, 'XTick',[3000 10000, 20000],'xticklabels', labs(end-2:end), 'YTick',[], 'Fontsize',12, 'YColor', 'none', 'XTickLabelRotation',45)

%% print the magnitude figure
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f,'all_filters_October2023','-dpdf','-r0')
%% plot RT

f = figure(2); clf

%% subplot(2,2,1)
[ax11, ax12] = subplotResponse(w, t60, t60_GEQ, [3055, 2900, 3000, fs/2],  [0 11],  pinks(1, :), [2,4,1]);


set(ax11, 'XTick',[1 30, 300 ],'xticklabels', labs(1:3), 'YTick',[0:2 :10], 'Fontsize',12, 'XTickLabelRotation',45)

ax11.YLabel.String = '$T_{60}$(s)';
ax11.YLabel.Interpreter = 'latex';

ax11.XLabel.String = '(a)';
ax11.XLabel.Interpreter = 'latex';

ax11.XLabel.Position(1) = ax11.XLabel.Position(1) +1800;
ax11.XLabel.Position(2) = ax11.XLabel.Position(2) -0.3;

ax11.Title.String = 'GEQ only';
ax11.Title.Interpreter = 'latex';
ax11.Title.Position(1) = 1800;


set(ax12, 'XTick',[3000 10000, 20000],'xticklabels', labs(end-2:end), 'YTick',[], 'Fontsize',12, 'YColor', 'none', 'XTickLabelRotation',45)


%% subplot(2,2,2)
[ax21, ax22] = subplotResponse(w, t60, t60_MED, [3055, 2900, 3000, fs/2],  [0 11],  blues(2, :), [2,4,3]);


set(ax21, 'XTick',[1 30, 300 ],'xticklabels', labs(1:3), 'YTick',[0:2 :10], 'Fontsize',12, 'XTickLabelRotation',45)

ax21.XLabel.String = '(b)';
ax21.XLabel.Interpreter = 'latex';

ax21.XLabel.Position(1) = ax21.XLabel.Position(1) +1800;
ax21.XLabel.Position(2) = ax21.XLabel.Position(2) -0.3;

ax21.Title.String = 'Median gain';
ax21.Title.Interpreter = 'latex';
ax21.Title.Position(1) = 1800;


set(ax22, 'XTick',[3000 10000, 20000],'xticklabels', labs(end-2:end), 'YTick',[], 'Fontsize',12, 'YColor', 'none', 'XTickLabelRotation',45)



%% subplot(2,2,3)
[ax31, ax32] = subplotResponse(w, t60, t60_NOT, [3055, 2900, 3000, fs/2],  [0 11],  blues(1, :), [2,4,5]);


set(ax31, 'XTick',[1 30, 300 ],'xticklabels', labs(1:3), 'YTick',[0:2 :10], 'Fontsize',12, 'XTickLabelRotation',45)

ax31.YLabel.String = '$T_{60}$(s)';
ax31.YLabel.Interpreter = 'latex';

ax31.XLabel.String = {'Frequency (Hz)', '(c)'};
ax31.XLabel.Interpreter = 'latex';

ax31.XLabel.Position(1) = ax31.XLabel.Position(1) +1800;
ax31.XLabel.Position(2) = ax31.XLabel.Position(2) -0.3;

ax31.Title.String = 'Notch filter';
ax31.Title.Interpreter = 'latex';
ax31.Title.Position(1) = 1800;


set(ax32, 'XTick',[3000 10000, 20000],'xticklabels', labs(end-2:end), 'YTick',[], 'Fontsize',12, 'YColor', 'none', 'XTickLabelRotation',45)

%% subplot(2,2,4)
[ax41, ax42] = subplotResponse(w, t60, t60_SHE, [3055, 2900, 3000, fs/2],  [0 11],  [1 0 0], [2,4,7]);

set(ax41, 'XTick',[1 30, 300 ],'xticklabels', labs(1:3), 'YTick',[0:2 :10], 'Fontsize',12, 'XTickLabelRotation',45)

ax41.YLabel.String = 'RT(s)';
ax41.YLabel.Interpreter = 'latex';

ax41.XLabel.String = {'Frequency (Hz)', '(d)'};
ax41.XLabel.Interpreter = 'latex';

ax41.XLabel.Position(1) = ax41.XLabel.Position(1) +1800;
ax41.XLabel.Position(2) = ax41.XLabel.Position(2) -0.3;

ax41.Title.String = 'Proposed';
ax41.Title.Interpreter = 'latex';
ax41.Title.Position(1) = 1800;


set(ax42, 'XTick',[3000 10000, 20000],'xticklabels', labs(end-2:end), 'YTick',[], 'Fontsize',12, 'YColor', 'none', 'XTickLabelRotation',45)


%% print the RT figure
%%%%%%%%%%%%%%%%%%%%%%
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f,'all_T60s_October2023','-dpdf','-r0')
%% biquad filter
labs_ = ["1", "10","30","100", "300","1k", "3k", "10k", "20k"];
diff_hs_ = target_mag - 20*log10(abs(Hshelf_L));
f = figure(3); clf 
%%%%%%%%%%%%%%%
subplot(2,1,1); 
semilogx(w,20*log10(abs(Hshelf_L)), 'r-', 'LineWidth',1.5); hold on
semilogx(w,diff_hs_, '-','color', blues(1, :), 'LineWidth',1.5)

xlabel('Frequency (Hz)', 'interpreter', 'latex')
ylabel('Magnitude (dB)', 'interpreter', 'latex')


set(gca, 'XTick',[1 10 30, 100 300 1000 3000 10000, 20000],'xticklabels', labs_, 'YTick',[-30:10: 10], 'Fontsize',12)
box on
grid off
% grid on
xlim([0 fs/2])
ylim([-22, 8])
legend( '$H_{\textrm{I}}  (\omega)$', '$H_{\textrm{II}}  (\omega)$', 'location','southwest', 'interpreter', 'latex', 'numcolumns', 2)


%%%%%%%%%%%%%%%%%%%
subplot(2,1,2)
semilogx(w,target_mag, 'k', 'LineWidth', 2); hold on
semilogx(w,20*log10(abs(HSHE)), '--','color', pinks(1, :),  'LineWidth',2)

xlabel('Frequency (Hz)', 'interpreter', 'latex')
ylabel('Magnitude (dB)', 'interpreter', 'latex')


set(gca, 'XTick',[1 10 30, 100 300 1000 3000 10000, 20000],'xticklabels', labs_, 'YTick',[-30:10: 10], 'Fontsize',12)
box on
grid off
% grid on
xlim([0 fs/2])
ylim([-22, 8])
legend('Target',  '$H_{\textrm{att}} (\omega)$',  'location','southwest', 'interpreter', 'latex', 'numcolumns', 2)
f.Position(end) = 320;
pos = f.Position;

%% 
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f,'H_I_H_II_H_att_Sept2023','-dpdf','-r0')

%% function for figure plotting


function [ax1, ax2] = subplotResponse(w, target, response, wlim,  ylims,  colors, plotxyz)

    ax1 = subplot(plotxyz(1), plotxyz(2), plotxyz(3)); % left side of the axis - logarithmic

    semilogx(w(w < wlim(1)),target(w <  wlim(1)), 'k--', 'LineWidth',2); hold on
    semilogx(w(w < wlim(1)),response(w < wlim(1)),'-', 'color', colors, 'LineWidth',2)

    xlim([1 wlim(3)])
    yline(ylims(2), '-k')

    ylim(ylims)
    
    pos = ax1.Position;
    pos(end) = 0.25;
    pos(2) = pos(2) +0.09;
    ax1.Position = pos;
    
    grid off
    box off
    
    
    
    ax2 = subplot(plotxyz(1), plotxyz(2), plotxyz(3)+1); % right side of the axis - linear


    plot(w(w > wlim(2)),target(w > wlim(2)), 'k--', 'LineWidth',2); hold on
    plot(w(w >wlim(2)),response(w > wlim(2)),'-', 'color', colors, 'LineWidth',2)
    set(ax2,'units','normalized','position',[pos(1)+pos(3) pos(2) pos(3) pos(4)])
    
    xline(wlim(4), 'k')
    yline(ylims(2), '-k')
    
    xlim([wlim(3) wlim(4)])
    ylim(ylims)
    
    
    grid off
    box off
    
end