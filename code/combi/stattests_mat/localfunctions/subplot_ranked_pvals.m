function subplot_ranked_pvals(p,FDR,alpha,ftr_names,ytext,topN)
    
    [B,I] = sort(p);
    ranked_pFDR = FDR(I);
    ranked_pval = B;
    ranked_rows = ftr_names(I);
    
    if isa(topN,'double') % only include top N features
        ranked_pFDR = ranked_pFDR(1:topN);
        ranked_pval = ranked_pval(1:topN);
        ranked_rows = ranked_rows(1:topN);
    elseif strcmp(topN,'sign_p')
        stop = find(ranked_pval>alpha,1,'first');
        if ~isempty(stop)
            ranked_pFDR = ranked_pFDR(1:stop);
            ranked_pval = ranked_pval(1:stop);
            ranked_rows = ranked_rows(1:stop);
        end
    elseif strcmp(topN,'sign_q')
        stop = find(ranked_pFDR>alpha,1,'first');
        if ~isempty(stop)
            ranked_pFDR = ranked_pFDR(1:stop);
            ranked_pval = ranked_pval(1:stop);
            ranked_rows = ranked_rows(1:stop);
        end
    end
    
    X = categorical(ranked_rows);
    X = reordercats(X,ranked_rows);
    
    hold on
    bar(X,[ranked_pval,ranked_pFDR]);
    yline(alpha)
%     legend('P values','pFDR','alpha') % call legend elsewhere
    ylabel(ytext)
    hold off
end