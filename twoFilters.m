function [Hatt, w, target_mag, h] = twoFilters(RT, dL, fs, method, wc)
% implementation on the two-stage filter for recursive systems in
% artificial reverberation synthesis

% uses functions:
%  low_shelf
%  high_shelf
%  aceq
%  acge3

% INPUT: 
% RT - reverberation time values for octave/third-octave bands [s]
% dL - delay-line length [samples]
% fs - sampling rate
% method - which filter to use in the first stage:
%           'GEQ' - no prefiltering
%           'median' - median gain
%           'notch' - median  + shelving filter in high frequencies 
%           'shelf' - shelving filter
% wc - crossover frequency for the shelf filter [Hz]

% OUTPUT:
% Hatt - frequency response of the attenuation filterfilter

%% step one - prepare the target filter gains

gdB = -60./(RT);            % reverberation time converted to target gains in dB
gains = 10.^(gdB/20);       % frequency-dependent linear gain
gLin_dl = gains.^(dL/fs);   % linear gains adjusted for the delay-line length
gdB_dl = 20*log10(gLin_dl); % delay-adjusted gains in dB

%% step two - interpolate the shape of the target attenuation
nBand  = length(RT);        % number of octave bands
Nfreq = 2^9;                % Number of frequency points for frequency response evaluation
w = [ logspace(log10(1),log10(fs/2-1),Nfreq-1), fs/2];  % frequency points for frequency response evaluation

if nBand == 10 % octave
f =  16000./(2.^(9:-1:0));
elseif nBand == 31 % third octave
 f =  10^3 * (2 .^ ([-17:13]/3)); 
elseif nBand == 30 % third octave minus the highest band
 f =  10^3 * (2 .^ ([-17:13]/3)); 
 f(end) = [];
end 

% 

ind = zeros(nBand, 1); % locations of the band frequencies
for i = 1:length(f)
    [~, ind(i)] = min(abs(w-f(i)));
end

target_mag =zeros(Nfreq,1); % this will be the target magnitude response
target_mag(ind(1):ind(end)) = interp1(f, gdB_dl, w(ind(1):ind(end)), "linear"); % interpolate linearly between known points
target_mag = interp1(w(ind(1):ind(end)), target_mag(ind(1):ind(end)), w, "nearest", "extrap"); % extrapolate to DC and Nyquist, nearest neighbor to avoid weird extrapolations
    


gLin_ = 10.^(target_mag./20);   % convert the full target frequency response into linear gains

%% use only GEQ
if contains(method, 'geq') || contains(method, 'GEQ') || contains(method, 'Geq')
    Gdb_d = target_mag(ind);
    flag = 0; % do not use broadband gain
    Hshelf = ones(Nfreq,1); % do not use shelf filter
end

%% use median gain and GEQ
if contains(method, 'median') || contains(method, 'Median') 
    Gdb_d = target_mag(ind);
    flag = 1; % use broadband gain
    Hshelf = ones(Nfreq,1);
end

%% use median gain and notch filter and GEQ
%  for reference, see : K. Prawda, S. J. Schlecht, and V. Välimäki. “Improved reverberation time control for feedback delay networks”. In: Proc. Int.
% Conf. Digital Audio Effects (DAFx). Birmingham, UK, Sept. 2019
if contains(method, 'notch') || contains(method, 'Notch') 
    Gdb_d = target_mag(ind);
    flag = 1;

    GH =10^((Gdb_d(end) - median(Gdb_d))/20); % shelf filter gain
    [Hnum, Hden] = high_shelf(GH, fs);
    Hshelf = freqz(Hnum, Hden, w, fs);
   
end

%% pre-filter using first-order shelf filter

% get the filter coefficients
if contains(method, 'shelf') || contains(method, 'Shelf')
    GL = gLin_(1);% the gain for low frequencies
    GH = gLin_(end); % the gain for high frequencies
    

   [Hnum, Hden] = low_shelf(wc,fs, GL, GH);

    % evaluate filter's frequency response
    Hshelf = freqz(Hnum, Hden, w, fs);
   
    % get the difference between the target magnitude and the shelf filter
    % magnitude
    diff_mag= target_mag-20*log10(abs(Hshelf));

    % get the target gains for the GEQ
    Gdb_d = diff_mag(ind);
    flag = 0;
end


%%



if nBand == 10
    [nums,dens, Hatt] = aceq(Gdb_d, flag, fs);
elseif nBand == 31 || 30
    [nums,dens, Hatt] = acge3(Gdb_d, flag, fs);
end 

h = impz(nums(:,1), dens(:, 1), Nfreq, fs);
for i = 2:nBand
    h = filter(nums(:, i), dens(:, i), h);
end
H_geq = Hatt; % save in another variable for plotting

if exist('Hnum') % if using the pre- or post-filter in the form of shelf filter
    h = filter(Hnum, Hden, h);
    Hatt = Hatt.*(Hshelf.');
end

%% plot the filter responses

figure(10); clf

semilogx(w, target_mag, 'k', 'LineWidth',2); hold on
semilogx(w, db(Hshelf), 'r:',  'LineWidth',1);

if contains(method, 'shelf') || contains(method, 'Shelf')
    semilogx(w, diff_mag, 'k--',  'LineWidth',1); hold on
end

semilogx(w, db(Hatt), 'r--',  'LineWidth',2); hold on
semilogx(w, db(H_geq), 'b--',  'LineWidth',2); hold on

xlim([1 fs/2])
set(gca, 'XTick',[10 30, 100 300 1000 3000 10000], 'Fontsize',12,'fontname','Times')
box on
grid on
end