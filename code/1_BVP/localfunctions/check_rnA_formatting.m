function check_rnA_formatting(sig,fs,rnA_idx,rnA_prv,fixed_length,filename)

% sig           PPG signal @fs
% fs            sampling frequency
% rnA_idx       annotation indices in [samples]
% rnA_prv       rough PRV signal [ms]
% fixed_length  [samples]
% filename      string


%% Default settings
edge_allowance  = 2;    % [s]
prv_lower_bound = 500;  % [ms]
prv_upper_bound = 1500; % [ms]


%% Run assertions, raise warning on catch

try
    assert(rnA_idx(end)<=fixed_length);
catch
    fprintf('Warning - %s only has %d samples \n',filename,length(sig));
end
    
try
    assert(rnA_idx(end)>fixed_length-fs*edge_allowance);
catch
    fprintf('Warning - %s only has %d samples \n',filename,length(sig));
end

try
    assert(~any(isnan(rnA_idx)));
catch
    fprintf('Warning - %s contains NaN values \n',filename);
end

try
    assert(isempty(find(rnA_prv>prv_upper_bound,1)));
catch
    fprintf('Warning - %s has PRV values over %d \n',filename,prv_upper_bound);
end

try
    assert(isempty(find(rnA_prv<prv_lower_bound,1)));
catch
    fprintf('Warning - %s has PRV values under %d \n',filename,prv_lower_bound);
end


end