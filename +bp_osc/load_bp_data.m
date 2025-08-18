function [data, meta] = load_bp_data(filepath, varname)
% LOAD_BP_DATA  Loads blood pressure data from a .mat file.
%
% Inputs:
%   filepath - Path to .mat file
%   varname  - Variable name in .mat file
%
% Outputs:
%   data - Data vector
%   meta - Struct with metadata

S = load(filepath);
if isfield(S, varname)
    data = S.(varname);
else
    error('Variable "%s" not found in file %s', varname, filepath);
end

meta.filepath = filepath;
meta.length   = numel(data);
meta.fs       = []; % Set manually later if not stored

end