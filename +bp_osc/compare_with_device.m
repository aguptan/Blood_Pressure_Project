function results = compare_with_device(MAP_est, device_systole, device_diastole)
% COMPARE_WITH_DEVICE  Compares estimated MAP to device MAP.
%
% Inputs:
%   MAP_est         - Vector of estimated MAP values (mmHg)
%   device_systole  - Vector of device systolic values (mmHg)
%   device_diastole - Vector of device diastolic values (mmHg)
%
% Outputs:
%   results - Table with comparison and error stats

MAP_device = round((2*device_diastole + device_systole) / 3);
err = MAP_est(:) - MAP_device(:);

results = table(MAP_device(:), MAP_est(:), err, ...
    'VariableNames', {'MAP_device', 'MAP_est', 'Error'});

end