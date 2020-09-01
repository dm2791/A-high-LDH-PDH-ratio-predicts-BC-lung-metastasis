% unsupervised analysis of metabolon data
%% load the metabolomics data from file breastnorm.xlsx
breast = readtable('breast.csv');

sampleName = {'P1' 'P2' 'P3' 'P4' 'P5'...
    'B1' 'B2' 'B3' 'B4' 'B5'...
    'L1' 'L2' 'L3' 'L4' 'L5'};

%% heatmap
breastmat=table2array(breast(10:end,15:end));
% Windows PCs and other systems may load the xlsx differently
% if the data is loaded as cell convert to numeric here
if ~isnumeric(breastmat)
    breastmat=cellfun(@str2num,breastmat);
end
clustergram(breastmat,'Symmetric','false','Standardize','row',...
    'Colormap','redbluecmap', 'ColumnLabels', sampleName);

%% get the data for the glycolysis pathway only
glucolysisIdx=find(strcmp(breast.Var4,'Glycolysis, Gluconeogenesis, and Pyruvate Metabolism'));
breastmatglyc=breast(glucolysisIdx,:);
breastmatglyc=table2array(breastmatglyc(:,15:end));

% Windows PCs and other systems may load the xlsx differently
% if the data is loaded as cell convert to numeric here
if ~isnumeric(breastmatglyc)
    breastmatglyc=cellfun(@str2num,breastmatglyc);
end

breastparentalmean=mean(breastmatglyc(:,1:5),2);
parentalFC=breastmatglyc(:,1:5)./breastparentalmean;
brainFC=breastmatglyc(:,6:10)./breastparentalmean;
lungFC=breastmatglyc(:,11:15)./breastparentalmean;
FCglyc=[parentalFC,brainFC,lungFC];
cgo=clustergram(FCglyc(:,6:end),'Symmetric','false',...
    'Colormap','redbluecmap', 'ColumnLabels', sampleName(6:end));
glycmets=breast.Var2(glucolysisIdx);
set(cgo,'RowLabels',glycmets)

%% get bar plot of glyc metabs
breastmatglyc_mean=[mean(breastmatglyc(:,1:5),2),mean(breastmatglyc(:,6:10),2),mean(breastmatglyc(:,11:15),2)];

% remove glucose from this plot
breastmatglyc_mean(5,:)=[];
glycmets(5) = [];

figure
bar(categorical(glycmets),breastmatglyc_mean);
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
h = gca;
h.XTickLabels = glycmets;
h.XTickLabelRotation = 90;

%% pca
[coeff,score,latent,tsquared,explained,mu] = pca(breastmat');
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
set(hPt(1:5), 'MarkerFaceColor', [0.1804    0.1922    0.5725],...
    'MarkerSize', 5, 'MarkerEdgeColor', [0 0 0 ])
set(hPt(6:10), 'MarkerFaceColor', [0.5765    0.1529    0.5608],...
    'MarkerSize', 5, 'MarkerEdgeColor', [0 0 0 ])
set(hPt(11:15), 'MarkerFaceColor', [0    0.6627    0.6157],...
    'MarkerSize', 5, 'MarkerEdgeColor', [0 0 0 ])
set(hPt([1,6,11]), 'Marker', 'o', 'MarkerSize', 10)
set(hPt([2,7,12]), 'Marker', '^', 'MarkerSize', 10)
set(hPt([3,8,13]), 'Marker', 'd', 'MarkerSize', 10)
set(hPt([4,9,14]), 'Marker', 's', 'MarkerSize', 10)
set(hPt([5,10,15]), 'Marker', 'v', 'MarkerSize', 10)

xlabel(sprintf('PC1 (%d%% variance explained)', round(explained(1))))
ylabel(sprintf('PC2 (%d%% variance explained)', round(explained(2))))

 
%% glucose inset
breastmatglyc_mean=[mean(breastmatglyc(:,1:5),2),mean(breastmatglyc(:,6:10),2),mean(breastmatglyc(:,11:15),2)];
figure,bar(breastmatglyc_mean(5,:));

%% plsr to compare brain and lung breast mets
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

% Windows PCs and other systems may load the xlsx differently
% if the data is loaded as cell convert to numeric here
if ~isnumeric(breastmat_top20)
    breastmat_top20=cellfun(@str2num,breastmat_top20);
end

cgo=clustergram(breastmat_top20,'Symmetric','false','Standardize','row','Colormap','redbluecmap');
set(cgo,'RowLabels',top20mets)
