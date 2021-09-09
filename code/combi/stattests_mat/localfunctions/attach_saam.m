function test = attach_saam(test,saam,pids)
% given a single sample from ftr or dftr, attach a SAAM score sample
% saam.var specifies which score to use (anx,sec,avd)
% saam.varnames links saam.var to .scores matrix
% saam.scores provides the scores for each subject
% --> read specified scores from data

% cf. build_samples.m:
% test.included_pids{ftr} tells us which pids are used in the given sample
% test.sample_rmv{ftr} provides info on missing data due to NaN elements
% --> reconstruct SAAM score sample that matches given sample
    
    %% basic prep

    test.nb_samples = test.nb_samples + 1;
    
    % read SAAM scores
    column = find(contains(saam.varnames,saam.var));
    scores_full = saam.scores(:,column);
    
    
    %% match given test.sample
    
    % check for missing pids (independent of features)
    include_rows = contains(pids.full,test.included_pids);
    scores_pidmatch = scores_full(include_rows);
    
    % adjust for invalid datapoints in each feature
    % then append the matched SAAM score sample
    scores_ftrmatch = cell(test.nb_features,1);
    for ftr = 1:test.nb_features
        scores_ftrmatch{ftr} = scores_pidmatch(~test.sample_rmv{ftr});
        test.sample{ftr} = [test.sample{ftr} , scores_ftrmatch{ftr}];
    end
end