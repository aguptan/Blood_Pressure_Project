function summary = summarize_results(segments, analyses, MAP_device)
% SUMMARIZE_RESULTS  Creates a summary table of all segments and results.
%
% Inputs:
%   segments   - Struct array of segments
%   analyses   - Struct array of analyses
%   MAP_device - Device MAP values
%
% Outputs:
%   summary - Table with trial info and MAP comparison

nSeg = numel(segments);
MAP_est = arrayfun(@(a) round(a.peakPressure), analyses);

summary = table((1:nSeg)', {segments.type}', MAP_device(:), MAP_est(:), ...
    MAP_est(:) - MAP_device(:), ...
    'VariableNames', {'Trial', 'Type', 'MAP_device', 'MAP_est', 'Error'});

end