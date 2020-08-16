% unsupervised analysis of metabolon data
%% load data
breast = readtable('breastnorm.xlsx');

%% heatmap
breastmat=table2array(breast(10:end,15:end));
breastmat=cellfun(@str2num,breastmat);
clustergram(breastmat,'Symmetric','false','Standardize','row','Colormap','redbluecmap');

%glycolysis pathway only
keepmets=find(strcmp(breast.Var4,'Glycolysis, Gluconeogenesis, and Pyruvate Metabolism'));
breastmatglyc=breast(keepmets,:);
breastmatglyc=table2array(breastmatglyc(:,15:end));
breastmatglyc=cellfun(@str2num,breastmatglyc);

breastparentalmean=mean(breastmatglyc(:,1:5),2);
parentalFC=breastmatglyc(:,1:5)./breastparentalmean;
brainFC=breastmatglyc(:,6:10)./breastparentalmean;
lungFC=breastmatglyc(:,11:15)./breastparentalmean;
FCglyc=[parentalFC,brainFC,lungFC];
cgo=clustergram(FCglyc(:,6:end),'Symmetric','false','Colormap','redbluecmap');
glycmets=breast.Var2(keepmets);
set(cgo,'RowLabels',glycmets)

% get bar plot of glyc metabs
breastmatglyc_mean=[mean(breastmatglyc(:,1:5),2),mean(breastmatglyc(:,6:10),2),mean(breastmatglyc(:,11:15),2)];
breastmatglyc_mean(5,:)=[];
bar(breastmatglyc_mean);
Bp=[];
Lp=[];
BLp=[];
for i=1:length(glycmets)
    [h,p] = ttest2(breastmatglyc(i,1:5),breastmatglyc(i,6:10));
    Bp(i)=p;
    [h,p] = ttest2(breastmatglyc(i,1:5),breastmatglyc(i,11:15));
    Lp(i)=p;
    [h,p] = ttest2(breastmatglyc(i,6:10),breastmatglyc(i,11:15));
    BLp(i)=p;
end

%% pca
breastmat=table2array(breast(10:end,15:end))';
breastmat=cellfun(@str2num,breastmat);
[coeff,score,latent,tsquared,explained,mu] = pca(breastmat);
group=repmat({'p'},1,5);
group=[group,repmat({'b'},1,5)];
group=[group,repmat({'l'},1,5)];

selectlabels=breast.Var2;
selectlabels(1:9)=[];
keepfew=[319];
keepfew2=zeros(645,1);
keepfew2(keepfew)=1;
selectlabels2=selectlabels;
selectlabels2(~keepfew2)={' '};
selectlabels3=selectlabels(keepfew);

 figure
 h=biplot(coeff(:,1:2),'Scores',score(:,1:2),'Varlabels',selectlabels2,'Color',[.7,.7,.7]);
hID = get(h, 'tag');
hPt = h(strcmp(hID,'obsmarker'));
set(hPt(1:5), 'Color', [0.1804    0.1922    0.5725],'MarkerSize', 5)
set(hPt(6:10), 'Color', [0.5765    0.1529    0.5608],'MarkerSize', 5)
set(hPt(11:15), 'Color', [0    0.6627    0.6157],'MarkerSize', 5)
set(hPt([1,6,11]), 'Marker', '.')
set(hPt([2,7,12]), 'Marker', '*')
set(hPt([3,8,13]), 'Marker', 'd')
set(hPt([4,9,14]), 'Marker', '+')
set(hPt([5,10,15]), 'Marker', 'p')

 
 % glucose inset
 breastmatglyc_mean=[mean(breastmatglyc(:,1:5),2),mean(breastmatglyc(:,6:10),2),mean(breastmatglyc(:,11:15),2)];
figure,bar(breastmatglyc_mean(5,:));

%% plsr to compare brain and lung breast mets
breastmat=table2array(breast(10:end,15:end));
breastmat=cellfun(@str2num,breastmat);
X=log([breastmat(:,6:10)';breastmat(:,11:15)']);
Y=[repmat(0,5,1);repmat(1,5,1)];
[Xloadings,Yloadings,Xscores,Yscores,betaPLS10,PLSPctVar,mse,stats] = plsregress(X,Y,9);

%select top 20% of metabolites contributing to comp1
[B,I] = sort(stats.W(:,1));
bottom10p=breast.Var11(I(1:65));
top10p=breast.Var11(I(end-64:end));
bottom10p=bottom10p(~strcmp(bottom10p,{''}));
top10p=top10p(~strcmp(top10p,{''}));

%convert to metab names
bottom10p_ind=ismember(breast.Var11,bottom10p);
top10p_ind=ismember(breast.Var11,top10p);
breastlungmets_downoverbrain=breast.Var2(bottom10p_ind);
breastlungmets_upoverbrain=breast.Var2(top10p_ind);

% put selected metabolites into pathways
lungmets_downpathways=table(breast.Var4(bottom10p_ind));
lungmets_downpathways=grpstats(lungmets_downpathways,'Var1');
lungmets_uppathways=table(breast.Var4(top10p_ind));
lungmets_uppathways=grpstats(lungmets_uppathways,'Var1');

%heatmaps
% heatmap of top 20% metabs only
breast_top20=breast(bottom10p_ind,:);
breast_top20=[breast_top20;breast(top10p_ind,:)];
top20mets=breast_top20.Var4;
breastmat_top20=table2array(breast_top20(:,20:end));
breastmat_top20=cellfun(@str2num,breastmat_top20);
cgo=clustergram(breastmat_top20,'Symmetric','false','Standardize','row','Colormap','redbluecmap');
set(cgo,'RowLabels',top20mets)
