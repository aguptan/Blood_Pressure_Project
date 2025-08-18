%% Algorithm Testing: BP Oscillometry + Simple Agreement Analysis (no plots)
tic; clc; close all; clearvars;
set(0, 'DefaultFigureWindowStyle', 'docked'); % For docked plots

%% -------------------- PARAMETERS --------------------
saveFigs = true; % Save the pictures or not

% Detect this script's folder
figuresFolder = fileparts(mfilename('fullpath'));

% File & data
dataFile   = ;
dataVar    = 'data';
fa         = 1000;  % Sampling rate (Hz)

% Filter parameters
Fcp        = 0.5;   % High-pass cutoff frequency (Hz)
Fsp        = 500;   % Filter design frequency (Hz)
filterOrder = 3;

% Envelope window sizes 
envWinDown = 400;

% Device readings (CVS) and selection labels you used when pairing
diastole   = [92; 94; 86; 85; 84; 85; 83; 90; 91; 86];
systole    = [122; 127; 127; 121; 116; 117; 129; 121; 132; 126];
trial_id   = [0,1,1,0,1,1,1,1,0,0,0,0]; 
MAP_device = round((2*diastole + systole) / 3);

% Manual start & end indices 
start_ncpl_down = [25015,  123907, 544616, 605601, 686961, 748412];
end_ncpl_down   = [49100,  146100, 573973, 632877, 711598, 776347];
[data, meta] = bp_osc.load_bp_data(dataFile, dataVar);
ncpldata = data(1:782170); 
[sos, g] = bp_osc.design_hp_filter(Fsp, Fcp, filterOrder);
segments = bp_osc.segment_cuff_cycles(ncpldata, fa, start_ncpl_down, end_ncpl_down, repmat({'down'}, 1, numel(start_ncpl_down)));

alg_map_down = zeros(numel(segments),1);
for block = 1:numel(segments)
    [analysis_down{block}, alg_map_down(block)] = bp_osc.analyze_oscillogram(segments(block), fa, sos, g, envWinDown);
end
figure
for block = 1:numel(segments)
    subplot(3,2,block)
    plot(segments(block).t, segments(block).y)
    sgtitle("Downstroke")
    xlabel('Time (s)'); ylabel('Pressure (mmHg)')
end
if saveFigs, saveas(gcf, fullfile(figuresFolder, 'non_coupled_down_raw.jpeg')); end

figure
for block = 1:numel(analysis_down)
    subplot(3,2,block)
    plot(analysis_down{block}.t2p, analysis_down{block}.envDiff)
    hold on
    plot(analysis_down{block}.t2p, analysis_down{block}.filtered)
    legend([num2str(round(analysis_down{block}.peakPressure)), ' mmHg'])
    sgtitle("Downstroke Oscillogram in time")
    xlabel('Pressure (mmHg)'); ylabel('\Delta Pressure (mmHg)')
    set(gca, 'xdir', 'reverse')
end
if saveFigs, saveas(gcf, fullfile(figuresFolder, 'non_coupled_down_filtered.jpeg')); end

%% -------------------- ANALYSIS --------------------
device_subset = MAP_device;
n_pairs = min(numel(alg_map_down), numel(device_subset));

% Use the numeric vector produced by analyze_oscillogram
alg_vec = double(alg_map_down(1:n_pairs));
dev_vec = double(device_subset(1:n_pairs));

% Regression / agreement analysis
err   = alg_vec - dev_vec;
bias  = mean(err);
sdE   = std(err, 0);
loaLo = bias - 1.96*sdE;
loaHi = bias + 1.96*sdE;
mae   = mean(abs(err));
rmse  = sqrt(mean(err.^2));

% Prepare for correlation and fit
X = dev_vec(:);   % Device MAP (x-axis)
y = alg_vec(:);   % Algorithm MAP (y-axis)
R  = corr(X, y);
R2 = R.^2;

% Summary
fprintf('Pairs used: %d\n', n_pairs);
fprintf('Bias (Alg - Device): %.2f mmHg\n', bias);
fprintf('SD of error: %.2f mmHg\n', sdE);
fprintf('95%% LoA: [%.2f, %.2f] mmHg\n', loaLo, loaHi);
fprintf('MAE: %.2f mmHg\n', mae);
fprintf('RMSE: %.2f mmHg\n', rmse);
fprintf('Pearson r: %.3f\n', R);

% Table
results_table = table((1:n_pairs)', X, y, err(:), 'VariableNames', {'Pair','MAP_device','MAP_alg','Error'});


% Least-squares fit
p    = polyfit(X, y, 1);
xfit = linspace(min(X), max(X), 100);
yfit = polyval(p, xfit);

% Plot
figure;
scatter(X, y, 40, 'filled'); hold on; grid on;
plot(xfit, yfit, 'LineWidth', 2);               % regression line
% plot(xfit, xfit, '--', 'LineWidth', 1);         % identity line 
xlabel('Device MAP (mmHg)');
ylabel('Algorithm MAP (mmHg)');
title('Air Downstroke: Algorithm vs Device MAP');
legend('Data','Linear fit','Identity (y = x)','Location','best');
xlim([90 110])
ylim([90 125])
% Text box with fit and R^2
eqn = sprintf('y = %.2fx + %.2f\nR^2 = %.3f', p(1), p(2), R2);
text(min(X), max(y), eqn, 'FontSize',10, 'Color','r', 'VerticalAlignment','top')

if saveFigs
    saveas(gcf, fullfile(figuresFolder, 'air_down_linear_regression.jpeg'));
end

toc

