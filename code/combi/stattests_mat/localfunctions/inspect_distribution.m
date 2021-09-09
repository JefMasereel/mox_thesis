function inspect = inspect_distribution(test,inspect,verbose)
    
    nb_samples  = test.nb_samples;
    nb_features = test.nb_features; 

    inspect.vartestn.p  = zeros(nb_features,1);
    inspect.vartestn.h  = zeros(nb_features,1);
    inspect.normality.p = zeros(nb_features,nb_samples);
    inspect.normality.h = zeros(nb_features,nb_samples);
    
    % prepare arrays for missingness counters
    inspect.sampsize = zeros(nb_features,nb_samples);
    inspect.nb_nans  = zeros(nb_features,nb_samples);
    
    % suppress warnings for small P values
    % note to self: use `warning('query','last')` to identify warnings
    warning('off','stats:jbtest:PTooSmall');        % jbtest, high Pval
    warning('off','stats:jbtest:PTooBig');          % jbtest, low Pval
    warning('off','stats:adtest:OutOfRangePLow');   % adtest, low Pval
    warning('off','stats:adtest:OutOfRangePHigh');  % adtest, high Pval

    for ftr = 1:nb_features

        % homogenous variance
        p = vartestn(test.sample{ftr},'TestType',inspect.vartestn.method,'Display','off');
        inspect.vartestn.p(ftr) = p;
        inspect.vartestn.h(ftr) = (p < inspect.vartestn.alpha);
        
        % normal distribution
        for i = 1:nb_samples
            if strcmp(inspect.normality.method,'kstest')
                [h,p] = kstest(test.sample{ftr}(:,i),'Alpha',inspect.normality.alpha);
            elseif strcmp(inspect.normality.method,'adtest')
                [h,p] = adtest(test.sample{ftr}(:,i),'Alpha',inspect.normality.alpha);
            elseif strcmp(inspect.normality.method,'jbtest')
                [h,p] = jbtest(test.sample{ftr}(:,i),inspect.normality.alpha);
            end
            inspect.normality.p(ftr,i) = p;
            inspect.normality.h(ftr,i) = h;
        end
        
        % count nan values to quantify missingness
        for i = 1:nb_samples
            inspect.sampsize(ftr,i) = length(test.sample{ftr}(:,i));
            inspect.nb_nans(ftr,i)  = sum(isnan(test.sample{ftr}(:,i)));
        end
    end
    
    % turn warnings for small P values on again
    warning('on','stats:jbtest:PTooSmall');        % jbtest, high Pval
    warning('on','stats:jbtest:PTooBig');          % jbtest, low Pval
    warning('on','stats:adtest:OutOfRangePLow');   % adtest, low Pval
    warning('on','stats:adtest:OutOfRangePHigh');  % adtest, high Pval

    % display rejection ratio of tests
    if verbose
        
        % Normality tests (average rejection for all samples)
        nb_candidates = (nb_features*(nb_samples));
        nb_rejections = sum(inspect.normality.h,'all');
        ratio = nb_rejections/nb_candidates;
        fprintf('%.2f%% of normality tests (%s) is rejected.\n',ratio*100,inspect.normality.method);
        
        % Homoskedasticity test (avg pval for 2 samples)
        nb_candidates = nb_features;
        nb_rejections = sum(inspect.vartestn.h,'all');
        ratio = nb_rejections/nb_candidates;
        fprintf('%.2f%% of homogenous variance tests (%s) is rejected.\n',ratio*100,inspect.vartestn.method);
    end
end