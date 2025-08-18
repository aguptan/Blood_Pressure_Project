function segments = segment_cuff_cycles(data, fs, start_idx, end_idx, types)
% SEGMENT_CUFF_CYCLES  Creates segment definitions from manual indices.
%
% Inputs:
%   data      - Full BP signal
%   fs        - Sampling rate (Hz)
%   start_idx - Vector of start indices
%   end_idx   - Vector of end indices
%   types     - Cell array of type labels ('up', 'down')
%
% Outputs:
%   segments  - Struct array with segment info

nSeg = numel(start_idx);
segments = struct('id', [], 'type', [], 't', [], 'y', []);

for k = 1:nSeg
    idxRange = start_idx(k):end_idx(k);
    segments(k).id   = k;
    segments(k).type = types{k};
    segments(k).t    = (0:numel(idxRange)-1) / fs;
    segments(k).y    = data(idxRange);
end

end