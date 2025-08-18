tic; close all; clc;

%% -------------------- PARAMETERS --------------------
saveFigs   = true;                            % Set false to skip saving JPEGs
outFolder  = fileparts(mfilename('fullpath'));% Save to this script's folder

% Filter settings (match your analysis)
fs   = 1000;          % Sampling rate (Hz)  [aka fa]
fc   = 0.5;           % High-pass cutoff (Hz)
ord  = 3;             % Butterworth order

% Synthetic-test settings
testDur   = 10;                         % seconds
testFreqs = [0.5 1 2 5 10 15];          % Hz (typical oscillometry band)
t_syn     = (0:1/fs:testDur-1/fs)';     % time vector

% Real-data file and a single known segment (re-using your indices)
dataFile  = ;
dataVar   = 'data';
seg_start = 25015;                       
seg_end   = 49100;                       

%% -------------------- FILTER DESIGN ----------------------------
% Design Butterworth HP and get b,a for generic filter()/filtfilt() use.
[b_hp, a_hp] = butter(ord, fc/(fs/2), 'high');
[z,p,k] = butter(ord, fc/(fs/2), 'high');
[sos,g] = zp2sos(z,p,k);

%% -------------------- THEORETICAL RESPONSE ---------------------
% Frequency response
nFFT = 4096;
[H, w]   = freqz(b_hp, a_hp, nFFT, fs);
mag_db    = 20*log10(abs(H));
ph_unwrap = unwrap(angle(H));                 % radians

% Numerical group delay 
f_vec   = w;                                  % freq vector 
omega   = 2*pi*f_vec;
dphi    = gradient(ph_unwrap);
domega  = gradient(omega);
tau_g_s = -dphi ./ domega;                    % seconds
tau_g_samples = tau_g_s * fs;                 % samples

% Plot magnitude, phase, and group delay
figure;
tiledlayout(3,1);
nexttile;
plot(f_vec, mag_db); grid on; xlim([0 20]);
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)'); title('HP Butterworth Magnitude');

nexttile;
plot(f_vec, ph_unwrap); grid on; xlim([0 20]);
xlabel('Frequency (Hz)'); ylabel('Phase (rad)'); title('Unwrapped Phase (single-pass)');

nexttile;
plot(f_vec, tau_g_samples); grid on; xlim([0 20]);
xlabel('Frequency (Hz)'); ylabel('Group delay (samples)'); title('Group Delay (single-pass)');
if saveFigs
    saveas(gcf, fullfile(outFolder, 'theoretical_bode_groupdelay.jpeg'));
end

% Measure empirical delay via cross-correlation.
results_syn = table('Size',[numel(testFreqs) 4], ...
    'VariableTypes',{'double','double','double','double'}, ...
    'VariableNames',{'freq_Hz','theory_delay_samples','emp_filter_delay_samp','emp_filtfilt_delay_samp'});

for i = 1:numel(testFreqs)
    f0 = testFreqs(i);
    x  = sin(2*pi*f0*t_syn);                   

    % Single-pass (phase-shifting)
    y1 = filter(b_hp, a_hp, x);

    % Zero-phase
    y2 = filtfilt(b_hp, a_hp, x);

    % Empirical delay via cross-correlation peak (y vs x)
    [c1, lags1] = xcorr(y1, x, 'normalized');
    [~, idx1]   = max(c1);
    lag1_samp   = lags1(idx1);

    [c2, lags2] = xcorr(y2, x, 'normalized');
    [~, idx2]   = max(c2);
    lag2_samp   = lags2(idx2);

    
    [~, iF] = min(abs(f_vec - f0));
    gd_theory = tau_g_samples(iF);

    results_syn.freq_Hz(i)                = f0;
    results_syn.theory_delay_samples(i)   = gd_theory;
    results_syn.emp_filter_delay_samp(i)  = lag1_samp;
    results_syn.emp_filtfilt_delay_samp(i)= lag2_samp;
end

disp('Delay results (samples):');
disp(results_syn);

% Quick overlay for one representative frequency 
repF = 5;
[~, idxRep] = min(abs(testFreqs - repF));
x_rep  = sin(2*pi*testFreqs(idxRep)*t_syn);
y1_rep = filter(b_hp, a_hp, x_rep);
y2_rep = filtfilt(b_hp, a_hp, x_rep);

figure;
plot(t_syn, x_rep,  'DisplayName','Raw'); hold on;
plot(t_syn, y1_rep, 'DisplayName','filter (single-pass)');
plot(t_syn, y2_rep, 'DisplayName','filtfilt (zero-phase)');
grid on; xlim([0 1]); % 1-second zoom for clarity
xlabel('Time (s)'); ylabel('Amplitude'); title(sprintf('Overlay at %.1f Hz', testFreqs(idxRep)));
legend('Location','best');
if saveFigs
    saveas(gcf, fullfile(outFolder, sprintf('synthetic_overlay_%gHz.jpeg', testFreqs(idxRep))));
end

%% -------------------- REAL DATA ------------------
% Load and slice one of your known segments to visualize edge effects to
S = load(dataFile, dataVar);
x_real_full = S.(dataVar);
x_seg = x_real_full(seg_start:seg_end);



t_seg = (0:numel(x_seg)-1)'/fs;

y1_seg = filter(b_hp, a_hp, x_seg);       % single-pass
y2_seg = filtfilt(b_hp, a_hp, x_seg);     % zero-phase

% Empirical delay via xcorr
[c1, l1] = xcorr(y1_seg, x_seg, 'normalized'); [~,i1]=max(c1); lag1 = l1(i1);
[c2, l2] = xcorr(y2_seg, x_seg, 'normalized'); [~,i2]=max(c2); lag2 = l2(i2);

fprintf('\nReal segment delays (samples): single-pass = %d, filtfilt = %d\n', lag1, lag2);

figure;
subplot(2,1,1);
plot(t_seg, x_seg - 100, 'DisplayName','Raw'); hold on;
plot(t_seg, y1_seg, 'DisplayName','filter'); grid on;
xlabel('Time (s)'); ylabel('Pressure'); title('Real Segment: Raw vs filter');
legend;

subplot(2,1,2);
plot(t_seg, x_seg - 100, 'DisplayName','Raw'); hold on;
plot(t_seg, y2_seg, 'DisplayName','filtfilt'); grid on;
xlabel('Time (s)'); ylabel('Pressure'); title('Real Segment: Raw vs filtfilt');
legend;

if saveFigs
    saveas(gcf, fullfile(outFolder, 'real_segment_overlays.jpeg'));
end
toc;
