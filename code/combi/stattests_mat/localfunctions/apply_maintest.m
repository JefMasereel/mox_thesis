function test = apply_maintest(test)

    assert(test.nb_samples==2);         % expand if/when needed
    nb_features = test.nb_features; 
    test.p  = zeros(nb_features,1);
    test.h  = zeros(nb_features,1);
    
    % for correlation tests, save rho (corr coefficient)
    test.rho = zeros(nb_features,1);

    for ftr = 1:nb_features
        
        A = test.sample{ftr}(:,1);
        B = test.sample{ftr}(:,2);
        
        if strcmp(test.method,'signrank')           % paired
            [p,h] = signrank(A,B);
        elseif strcmp(test.method,'ranksum')        % independent
            [p,h] = ranksum(A,B);
        elseif strcmp(test.method,'ttest')          % paired
            [h,p] = ttest(A,B);
        elseif strcmp(test.method,'ttest2')         % independent
            [h,p] = ttest2(A,B);
        elseif strcmp(test.method,'pearson_corr')
            [rho,p] = corr(A,B,'type','Pearson','rows','pairwise');
            h = (p<0.05); 
        elseif strcmp(test.method,'spearman_corr')
            [rho,p] = corr(A,B,'type','Spearman','rows','pairwise');
            h = (p<0.05); 
        else
            fprintf('test.method = "%s" unknown, review test specs\n',test.method);
        end
        
        test.p(ftr) = p;
        test.h(ftr) = h;
        
        try test.rho(ftr) = rho; end
        
    end
    
    nb_rejects = sum(test.h);
    percentage = nb_rejects/nb_features*100;
    fprintf('Test as specified rejects H0 for %d of %d features (%.2f%%)\n',...
            nb_rejects,nb_features,percentage);
end