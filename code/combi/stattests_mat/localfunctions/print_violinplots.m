function print_violinplots(log,idx,savedir)

% builds an overview figure per feature
% 4 subplots, with 2 elements each

% define idx for violinplot data samples
vd.idx=cell(4,2);
vd.idx{1,1} = (idx.RS  & idx.PRE  & idx.oxytocin);    
vd.idx{1,2} = (idx.RS  & idx.POST & idx.oxytocin);   
vd.idx{2,1} = (idx.RS  & idx.PRE  & idx.placebo);     
vd.idx{2,2} = (idx.RS  & idx.POST & idx.placebo);    
vd.idx{3,1} = (idx.EYE & idx.PRE  & idx.oxytocin);   
vd.idx{3,2} = (idx.EYE & idx.POST & idx.oxytocin);  
vd.idx{4,1} = (idx.EYE & idx.PRE  & idx.placebo);    
vd.idx{4,2} = (idx.EYE & idx.POST & idx.placebo);   
vd.ctg = {'PRE','POST'};

% build figure for each feature
for ftr = 1:length(log.ftr_names)
    
    disp(log.ftr_names{ftr})
    
    % pad with nan to get equal length arrays, read specified feature
    vd.mat = nan(length(log.filenames),8);
    for i=1:4
        for j = 1:2
            idx_ij = vd.idx{i,j};
            vd.mat(idx_ij,2*(i-1)+j) = log.ftr_values(idx_ij,ftr);
        end
    end

    % get ampl window for figure
    ylim_upp = ceil(max(vd.mat,[],'all'));
    ylim_low = floor(min(vd.mat,[],'all'));
    
    % build figure
    figure('visible','off','name',['Violinplot for ' log.ftr_names{ftr}]); 
    set(gcf, 'Units', 'centimeters', 'Position', [0 0 12 5], 'PaperUnits', 'centimeters', 'PaperSize', [12 5]);

    subplot(1,4,1)
    violinplot(vd.mat(:,1:2),vd.ctg,'ShowNotch',true,'ShowMean',true);
    ylim([ylim_low,ylim_upp])
    title('RS oxytocin')
    subplot(1,4,2)
    violinplot(vd.mat(:,3:4),vd.ctg,'ShowNotch',true,'ShowMean',true);
    ylim([ylim_low,ylim_upp])
    title('RS placebo')
    subplot(1,4,3)
    violinplot(vd.mat(:,5:6),vd.ctg,'ShowNotch',true,'ShowMean',true);
    ylim([ylim_low,ylim_upp])
    title('EYE oxytocin')
    subplot(1,4,4)
    violinplot(vd.mat(:,7:8),vd.ctg,'ShowNotch',true,'ShowMean',true);
    ylim([ylim_low,ylim_upp])
    title('EYE placebo')
    
    % save figure as .tif
    filepath = [savedir log.ftr_names{ftr}];
    print(gcf,filepath,'-dtiff','-r300');
end

end


