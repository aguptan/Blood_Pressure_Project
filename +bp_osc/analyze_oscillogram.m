function [analysis, MAP_est] = analyze_oscillogram(segment, fs, sos, g, envWin)
% ANALYZE_OSCILLOGRAM  Filters and computes envelope for a BP segment,
% and returns a numeric MAP estimate from the envelope peak.
%
% Inputs:
%   segment - Struct with fields t (time) and y (pressure signal)
%   fs      - Sampling frequency (Hz)
%   sos, g  - High-pass filter parameters
%   envWin  - Envelope window length (samples)
%
% Outputs:
%   analysis - Struct with filtered signal, envelopes, and peak info
%   MAP_est  - Estimated Mean Arterial Pressure (mmHg)

% --- Filter ---
filtSig = filtfilt(sos, g, segment.y);

% --- Envelope ---
[envUp, envLo] = envelope(filtSig, envWin, 'peak');
envDiff = envUp - envLo;

% --- Map index to pseudo-pressure (linear mapping start→end) ---
t2p = linspace(segment.y(1), segment.y(end), numel(envDiff));

% --- Peak ---
[peakVal, peakIdx] = max(envDiff);

% --- MAP estimate ---
MAP_est = round(t2p(peakIdx));

% --- Store results ---
analysis = struct();
analysis.t            = segment.t;
analysis.raw          = segment.y;
analysis.filtered     = filtSig;
analysis.envUp        = envUp;
analysis.envLo        = envLo;
analysis.envDiff      = envDiff;
analysis.t2p          = t2p;
analysis.peakVal      = peakVal;
analysis.peakIdx      = peakIdx;
analysis.peakPressure = MAP_est;   % same value
analysis.MAP_est      = MAP_est;   % convenience duplicate
end
