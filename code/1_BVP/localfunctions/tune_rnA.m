function [nA_new,RMSXD] = tune_rnA(nA,sig,windowsize,fixed_length)

%% adjust annotations to ensure they're on the primary peak
% get nA_org from ex51_validate_nA.m
% tune them by finding local max around each annotation

% nA            annotation point indices [samples]
% sig           PPG signal @fs
% fs            sampling freq [Hz]
% windowsize    size of search window centered around nA [samples]

% nA_new        adjusted nA annotations [samples]
% RMSXD         root mean square difference on X-axis [scalar]


%% Apply shifting window to adjust each nA point

nA_new = zeros(size(nA));

RMSXD = 0; % RMS on x-axis deviation, to quantify adjustment

for i = 1:length(nA)
    
    if nA(i) - 0.5*windowsize <= 1
        start = 1;
    else
        start = nA(i) - 0.5*windowsize;
    end
    
    if nA(i) + 0.5*windowsize >= fixed_length
        stop = fixed_length;
    else
        stop = nA(i) + 0.5*windowsize;
    end

    local = sig(round(start:stop));
    [~,ref_idx] = max(local);
    nA_new(i) = nA(i)+ref_idx-(windowsize/2+1);
    RMSXD = RMSXD + (ref_idx-(windowsize/2+1))^2;
    
end

end