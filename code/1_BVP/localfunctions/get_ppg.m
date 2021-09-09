function ppg = get_ppg(ppg_dir,filename,fs_org,fs_int,fixed_size)

% ppg_dir       directory of target data
% filename      string specifier for target signal
% fs_org        sample frequency of original measurements
% fs_int        target fs to interpolate to (skip if nan)
% fixed_size    standardize sample size of measurements

% ppg           PPG signal of fixed length @fs_int

% WARNING
% check all format assumptions! Customize for your use case.
% .txt suffix is added to filename by default


%% Apply steps

tmp = readtable([ppg_dir filename '.txt']);
sig_org = tmp.Var3;
t_org = 0:1/fs_org:(length(sig_org)-1)/fs_org;

if isnan(fs_int)
    % continue as is
    sig = sig_org;
else
    % resample signal to fs_int
    t_int = t_org(1) : 1/fs_int : t_org(end);
    sig = spline(t_org, sig_org, t_int);
end

try
    ppg = sig(1:fixed_size);
catch
    fprintf('Warning - %s only has %d samples \n',filename,length(sig))
    ppg = sig;
end


end