function test = build_samples(test,pids,filenames)
% note: also compatible with dftr, but requires adjusted filenames var
% ftr  -> full filenames cell array
% dftr -> filenames = rownames_ext from dftr_ext() (line 41)
    
    %% basic prep
    
    % logical indexing arrays required for correct sample construction!
    assert(islogical(test.idx));
    
    test.nb_samples    = size(test.idx,2);
    test.nb_features   = size(test.data,2);
    test.nb_datapoints = size(test.idx,1);
    test.sample        = cell(test.nb_features,1);

    %% pairwise removal of incomplete subject data
    % "removal" = overwrite idx(i,j) with false (no ID mismatch caused)
    test.included_pids = {};
    test.excluded_pids = {};
    if test.paired==true
        for i = 1:length(pids.full)
            
            % pid = subject/patient ID (cfr filenames)
            % pid set = all datapoints (idx) related to this pid
            pid = pids.full{i};
            pid_set = contains(filenames,pid);
            
            if ~all(any(test.idx(pid_set,:),1))
                % remove broken pairs (no idx in 1+ cols of pid set)
                test.excluded_pids{end+1} = pid;
                test.idx(pid_set,:) = false;
            else
                % include idx for complete/valid pairs
                test.included_pids{end+1} = pid;
            end
        end
    end
    
    %% build & match samples
    
    test.sampsize = nan(test.nb_features,test.nb_samples);
    
    for ftr = 1:test.nb_features
        
        % basic prep, test.sample empty array
        % add minimal nan padding for unpaired samples
        nb_rows = max(sum(test.idx,1));
        test.sample{ftr} = nan(nb_rows,test.nb_samples);
        
        for i = 1:test.nb_samples
            sample = test.data(test.idx(:,i),ftr);
            samplesize = length(sample);
            test.sample{ftr}(1:samplesize,i) = sample;
            test.sampsize(ftr,i) = samplesize;
        end

        % pairwise removal of NaN values (if paired)
        % "removal" = drop rows that contain nan values
        % no new ID mismatches, be aware of this!
        % rmv is logged to enable later matching if needed
        if test.paired==true
            [test.sample{ftr},rmv] = rmmissing(test.sample{ftr});
            test.sample_rmv{ftr} = rmv; % 1 if removed
        end
    end
end