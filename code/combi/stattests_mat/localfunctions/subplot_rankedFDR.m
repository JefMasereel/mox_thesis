function subplot_rankedFDR(p,FDR,alpha,ftr_names,ytext,topN)
    
    [B,I] = sort(FDR); % adjust here if you want to sort by p-values
    
    ranked_pFDR = B;
    ranked_pval = p(I);
    ranked_rows = ftr_names(I);
    
    if topN ~= 'all' % only include top N features
        ranked_pFDR = ranked_pFDR(1:topN);
        ranked_pval = ranked_pval(1:topN);
        ranked_rows = ranked_rows(1:topN);
    end
    
    X = categorical(ranked_rows);
    X = reordercats(X,ranked_rows);
    
    hold on
    bar(X,[ranked_pval,ranked_pFDR]);
    yline(alpha)
    legend('P values','pFDR','alpha')
    ylabel(ytext)
    hold off
end