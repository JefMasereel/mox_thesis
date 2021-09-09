function [idx,pids] = define_hypotheses(filenames,placebo)

% Definition of hypotheses ported to external script for better overview
% Builds sample indices for all defined hypotheses as needed in stattests


% easier indexing formats
pids.oxytocin = placebo.labels.pid(contains(placebo.labels.treatment,'oxytocin'));
pids.placebo  = placebo.labels.pid(contains(placebo.labels.treatment,'placebo'));
idx.oxytocin  = contains(filenames,pids.oxytocin);
idx.placebo   = contains(filenames,pids.placebo);
idx.RS   = contains(filenames,'_RS_');
idx.EYE  = contains(filenames,'_EYE_');
idx.PRE  = contains(filenames,'PRE');
idx.POST = contains(filenames,'POST');

% % remove specific patients from entire analysis if wanted
% ignore_pids = [1 2 3 ...];
% idx_names = {'RS','EYE','PRE','POST','oxytocin','placebo'};
% for i = 1:length(idx_names)
%     tmp = idx.(idx_names{i});
%     for p = 1:length(ignore_pids)
%         pid_idx = ignore_pids(p);
%         pid_str = placebo.labels.pid{pid_idx};
%         rmv_idx = find(contains(filenames,pid_str));
%         tmp(rmv_idx) = 0;
%     end
%     idx.(idx_names{i}) = tmp;
% end

% % To define hypothesis, use these formats
% [State,Time,Treatment]
% paired (P) or non-paired (NP)

% Hyp1: for each ftr_value, test [RS,PRE-POST,oxytocin] (P)
idx.Hyp1.PRE  = (idx.RS & idx.PRE  & idx.oxytocin);
idx.Hyp1.POST = (idx.RS & idx.POST & idx.oxytocin);

% Hyp2: for each ftr_value, test [RS,PRE-POST,placebo]  (P)
idx.Hyp2.PRE  = (idx.RS & idx.PRE  & idx.placebo);
idx.Hyp2.POST = (idx.RS & idx.POST & idx.placebo);

% Hyp3: for each ftr_value, test [RS,PRE,oxytocin-placebo] (NP)
idx.Hyp3.OT = (idx.RS & idx.PRE & idx.oxytocin);
idx.Hyp3.PL = (idx.RS & idx.PRE & idx.placebo);

% Hyp4: for each ftr_value, test [EYE,PRE,oxytocin-placebo] (NP)
idx.Hyp4.OT = (idx.EYE & idx.PRE & idx.oxytocin);
idx.Hyp4.PL = (idx.EYE & idx.PRE & idx.placebo);

% add more if needed
% ...


end


