function [sos, g] = design_hp_filter(fs, fc, order)
% DESIGN_HP_FILTER  Designs a high-pass Butterworth filter.
%
% Inputs:
%   fs    - Sampling frequency (Hz)
%   fc    - Cutoff frequency (Hz)
%   order - Filter order
%
% Outputs:
%   sos   - Second-order sections matrix
%   g     - Gain factor

[z, p, k] = butter(order, fc / (fs/2), 'high');
[sos, g]  = zp2sos(z, p, k);

end