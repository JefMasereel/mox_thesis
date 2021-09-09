function mask = artf_mask(rsp,method,windowsize)

%% args

% mask   logical array
%               1 where artefacts were detected by specified method
%               0 where signal shows no significant deviations

% rsp           respiratory signal
% windowsize    size of sliding window, [samples]

% for clipping method, use
%   method.name = 'clip';
%   method.CF = clipfactor;
% for stdev method, use
%   method.name = 'stdev';
% with no further fields required. 


%% Collect local metrics as specified

localmetrics = zeros(1,length(rsp)-windowsize);

% apply clipping method
if isequal(method.name,'clip')
    
    % prepare clipping metric
    [~,clipflag,~] = clip_rsp(rsp,2);  %% CF=2, finetuning required!
    
    % collect local clipping ratios
    for i = 1:length(rsp)-windowsize
        window = i:i+windowsize;
        clips_in_segment = length(find(clipflag(window)==1));
        localratio = clips_in_segment/windowsize;
        localmetrics(i) = localratio;
    end

% Apply stdev method
elseif isequal(method.name,'stdev')
    
    % collect local stdev values
    for i = 1:length(rsp)-windowsize
        window = i:i+windowsize;
        localmetrics(i) = std(rsp(window));
    end
    
% Unknown method specifier
else
    error('unknown method specifier, check method struct');
end


%% Detect outliers and build mask

% estimate outliers
% detection method specified by method.rmoutliers
[~,rmbool] = rmoutliers(localmetrics,method.rmoutliers); 

% derive mask
mask = zeros(size(rsp));
for i = 1:length(rmbool)-windowsize
    window = i:i+windowsize;
    if any(rmbool(window))
        mask(window) = 1;
    end
end

end