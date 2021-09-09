function [rsp,flag,ratio] = clip_rsp(rsp,factor)

% IN 
% rsp     provided respiration signal
% factor  specified factor for clipping threshold

% OUT
% rsp     returns clipped respiration signal
% ratio   ratio (0..1) of clipped samples in signal


%% basic amplitude clipping

CF = factor;
MN = mean(rsp);
SD = std(rsp);
idx_upper = rsp>MN+CF*SD;
idx_lower = rsp<MN-CF*SD;
if ~isempty(find(idx_upper | idx_lower,1))
    ratio = length(find(idx_upper | idx_lower))/length(rsp);
    rsp(idx_upper) = MN+CF*SD;
    rsp(idx_lower) = MN-CF*SD;
else
    % no clipping necessary at this threshold
    % keep rsp as is, ratio is zero
    ratio = 0;
end

% 1 where sample was clipped, zero elsewhere
flag = idx_upper | idx_lower;

end