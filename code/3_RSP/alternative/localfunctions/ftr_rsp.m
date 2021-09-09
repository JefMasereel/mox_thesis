function ftr = ftr_rsp(bmObj)

% read baseline corrected respiration signal from bmObj
% estimate basic features based on ref015 etc

% bmObj     objoect class defined by BreathMetrics
%           WARNING: assumed that simplify==1


%% Estimate breath-per-breath biomarker scores

idx_inh = bmObj.inhaleOnsets;
idx_exh = bmObj.exhaleOnsets;
time    = bmObj.time;
ampl    = bmObj.baselineCorrectedRespiration;

% due to effects of simplify, ignore last exhale onset
loopsize = length(idx_inh)-1;
Ti = zeros(1,loopsize);
Te = zeros(1,loopsize);
RR = zeros(1,loopsize);
Vt = zeros(1,loopsize);
MV = zeros(1,loopsize);

for br = 1:loopsize
    Ti(br) = time(idx_exh(br)) - time(idx_inh(br));
    Te(br) = time(idx_inh(br+1)) - time(idx_exh(br));
    RR(br) = 60/(Ti(br)+Te(br));
    Vt(br) = ampl(idx_exh(br)) - ampl(idx_inh(br));
    MV(br) = RR(br)*Vt(br);
end


%% Compute cumulative scores

ftr.scorenames = {'MN','SD','MAD','CV','AR'};
ftr.biomarkers = {'Ti','Te','RR','Vt','MV'};
ftr.full.Ti = Ti;
ftr.full.Te = Te;
ftr.full.RR = RR;
ftr.full.Vt = Vt;
ftr.full.MV = MV;

for i=1:length(ftr.biomarkers)
    tmp = ftr.full.(ftr.biomarkers{i});
    
    if length(tmp) <= 1
        [MN SD MAD CV AR] = deal(NaN);
        warning = true;
        note = 'insufficient breaths detected in segment';
        fprintf('Warning - %d breath(s) detected in segment, skipping ftrs \n',length(tmp));
    else
        MN = mean(tmp);
        SD = std(tmp);
        MAD = mad(tmp);
        CV = SD/MN;
        AR = corr(tmp(1:end-1)',tmp(2:end)');
        warning = false;
        note = '';
    end
    
    ftr.MN.(ftr.biomarkers{i})  = MN;
    ftr.SD.(ftr.biomarkers{i})  = SD;
    ftr.MAD.(ftr.biomarkers{i}) = MAD;
    ftr.CV.(ftr.biomarkers{i})  = CV;
    ftr.AR.(ftr.biomarkers{i})  = AR;
    ftr.warning.(ftr.biomarkers{i}) = {warning;note};
end


end