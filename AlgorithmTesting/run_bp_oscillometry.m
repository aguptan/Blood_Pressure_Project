%% Algorithm Testing: BP Oscillometry
% load raw cuff data, define segments, design a high-pass filter, compute oscillogram envelopes, pick peaks,
% and plot/save figures for upstroke/downstroke blocks

tic; clc; close all; clearvars;
set(0, 'DefaultFigureWindowStyle', 'docked'); % For docked plots

%% -------------------- PARAMETERS --------------------
saveFigs = true; % Save the pictures or not

% Detect this script's folder
figuresFolder = fileparts(mfilename('fullpath'));

% File & data
dataFile   = DataFileLocation; 
dataVar    = 'data';
fa         = 1000;  % Sampling rate (Hz)

% Filter parameters
Fcp        = 0.5;   % High-pass cutoff frequency (Hz)
Fsp        = 500;   % Filter design frequency (Hz)
filterOrder = 3;

% Envelope window sizes (samples)
envWinDown = 400;
envWinUp   = 500;

% Device readings
diastole   = [92; 94; 86; 85; 84; 85; 83; 90; 91; 86];
systole    = [122; 127; 127; 121; 116; 117; 129; 121; 132; 126];
MAP_device = round((2*diastole + systole) / 3);

% Manual start/end indices
start_ncpl_down = [25015,  123907, 544616, 605601, 686961, 748412];
end_ncpl_down   = [49100,  146100, 573973, 632877, 711598, 776347];
start_dev_up    = [2348800, 2462880, 2555810, 2659740];
end_dev_up      = [2373640, 2487660, 2579280, 2683200];
start_dev_down  = [2393550, 2512030, 2603540, 2706130];
end_dev_down    = [2413810, 2530890, 2618940, 2724360];

%% -------------------- LOAD DATA --------------------
[data, meta] = bp_osc.load_bp_data(dataFile, dataVar);
ncpldata = data(1:782170); % Non-coupled subset

%% -------------------- DESIGN FILTER --------------------
[sos, g] = bp_osc.design_hp_filter(Fsp, Fcp, filterOrder);

%% -------------------- NON-COUPLED DOWNSLOPE --------------------
segments = bp_osc.segment_cuff_cycles(ncpldata, fa, start_ncpl_down, end_ncpl_down, repmat({'down'}, 1, numel(start_ncpl_down)));
for block = 1:numel(segments)
    analysis_down{block} = bp_osc.analyze_oscillogram(segments(block), fa, sos, g, envWinDown);
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

%% -------------------- DEVICE UPSTROKE --------------------
segments_up = bp_osc.segment_cuff_cycles(data, fa, start_dev_up, end_dev_up, repmat({'up'}, 1, numel(start_dev_up)));
for block = 1:numel(segments_up)
    analysis_up{block} = bp_osc.analyze_oscillogram(segments_up(block), fa, sos, g, envWinUp);
end

figure
for block = 1:numel(segments_up)
    subplot(2,2,block)
    plot(segments_up(block).t, segments_up(block).y)
    sgtitle("Upstroke")
    xlabel('Time (s)'); ylabel('Pressure (mmHg)')
end
if saveFigs, saveas(gcf, fullfile(figuresFolder, 'device_up_raw.jpeg')); end

figure
for block = 1:numel(analysis_up)
    subplot(2,2,block)
    plot(segments_up(block).t, analysis_up{block}.filtered)
    hold on
    plot(segments_up(block).t, analysis_up{block}.envLo, segments_up(block).t, analysis_up{block}.envUp)
    sgtitle("Upstroke")
    xlabel('Time (s)'); ylabel('Pressure (mmHg)')
end
if saveFigs, saveas(gcf, fullfile(figuresFolder, 'device_up_envelopes.jpeg')); end

figure
for block = 1:numel(analysis_up)
    subplot(2,2,block)
    plot(analysis_up{block}.t2p, analysis_up{block}.envDiff)
    hold on
    plot(analysis_up{block}.t2p, analysis_up{block}.filtered)
    plot(analysis_up{block}.t2p(analysis_up{block}.peakIdx), analysis_up{block}.peakVal, 'ro')
    legend([num2str(round(analysis_up{block}.peakPressure)), ' mmHg'])
    sgtitle("Upstroke Oscillogram in time")
    xlabel('Pressure (mmHg)'); ylabel('\Delta Pressure (mmHg)')
    set(gca, 'xdir', 'reverse')
end
if saveFigs, saveas(gcf, fullfile(figuresFolder, 'device_up_filtered.jpeg')); end

%% -------------------- DEVICE DOWNSLOPE --------------------
segments_dev_down = bp_osc.segment_cuff_cycles(data, fa, start_dev_down, end_dev_down, repmat({'down'}, 1, numel(start_dev_down)));
for block = 1:numel(segments_dev_down)
    analysis_dev_down{block} = bp_osc.analyze_oscillogram(segments_dev_down(block), fa, sos, g, envWinUp);
end

figure
for block = 1:numel(segments_dev_down)
    subplot(2,2,block)
    plot(segments_dev_down(block).t, segments_dev_down(block).y)
    sgtitle("Downstroke Device")
    xlabel('Time (s)'); ylabel('Pressure (mmHg)')
end
if saveFigs, saveas(gcf, fullfile(figuresFolder, 'device_down_raw.jpeg')); end

figure
for block = 1:numel(analysis_dev_down)
    subplot(2,2,block)
    plot(segments_dev_down(block).t, analysis_dev_down{block}.filtered)
    hold on
    plot(segments_dev_down(block).t, analysis_dev_down{block}.envLo, segments_dev_down(block).t, analysis_dev_down{block}.envUp)
    sgtitle("Downstroke Device")
    xlabel('Time (s)'); ylabel('Pressure (mmHg)')
end
if saveFigs, saveas(gcf, fullfile(figuresFolder, 'device_down_envelopes.jpeg')); end

figure
for block = 1:numel(analysis_dev_down)
    subplot(2,2,block)
    plot(analysis_dev_down{block}.t2p, analysis_dev_down{block}.envDiff)
    hold on
    plot(analysis_dev_down{block}.t2p(analysis_dev_down{block}.peakIdx), analysis_dev_down{block}.peakVal, 'ro')
    legend(num2str(analysis_dev_down{block}.peakPressure))
    sgtitle("Downstroke Device Oscillogram")
    xlabel('Pressure (mmHg)'); ylabel('Pressure (mmHg)')
end
if saveFigs, saveas(gcf, fullfile(figuresFolder, 'device_down_peaks.jpeg')); end

figure
for block = 1:numel(analysis_dev_down)
    subplot(2,2,block)
    plot(analysis_dev_down{block}.t2p, analysis_dev_down{block}.envDiff)
    hold on
    plot(analysis_dev_down{block}.t2p, analysis_dev_down{block}.filtered)
    plot(analysis_dev_down{block}.t2p(analysis_dev_down{block}.peakIdx), analysis_dev_down{block}.peakVal, 'ro')
    legend([num2str(round(analysis_dev_down{block}.peakPressure)), ' mmHg'])
    sgtitle("Downstroke Oscillogram in time")
    xlabel('Pressure (mmHg)'); ylabel('\Delta Pressure (mmHg)')
    set(gca, 'xdir', 'reverse')
end
if saveFigs, saveas(gcf, fullfile(figuresFolder, 'device_down_filtered.jpeg')); end

toc
