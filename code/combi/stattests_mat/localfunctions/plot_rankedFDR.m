function plot_rankedFDR(test,log)
    [B,I] = sort(test.FDR);
    printable_ftrs = strrep(log.ftr_names,'_','\_');
    sorted_ftrs = printable_ftrs(I);
    X = categorical(sorted_ftrs);
    X = reordercats(X,sorted_ftrs);

    figure
    hold on
    bar(X,[test.p(I),B]);
    yline(0.05)
    legend('P values','FDR','risk tolerance')
    title('FDR estimations for pvals in ascending order')
end