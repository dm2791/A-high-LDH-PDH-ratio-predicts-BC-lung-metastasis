%% analyze mithril output
%% waterfall plots
lungout=readtable('output_mithril_lung.txt');
ranked_lungout=sortrows(lungout,3);
figure
b=bar(ranked_lungout.RawAccumulator);

%select non significant and grey out
nonsig=find(ranked_lungout.pValue>.05);
b.FaceColor = 'flat';
b.CData(nonsig,:) = repmat([.7 .7 .7],length(nonsig),1);

%highlight glycolysis
glyc=find(contains(ranked_lungout.PathwayName,'Glycolysis'));
b.CData(glyc,:) = [0 0 0];

brainout=readtable('mithriloutputbrainh.xls');
ranked_lungout=sortrows(lungout,3);
figure
b=bar(ranked_lungout.RawAccumulator);
nonsig=find(ranked_lungout.pValue>.05);
b.FaceColor = 'flat';
b.CData(nonsig,:) = repmat([.7 .7 .7],length(nonsig),1);
glyc=find(contains(ranked_lungout.PathwayName,'Glycolysis'));
b.CData(glyc,:) = [0 0 0];
%% heatmap
pathway_both=lungout.PathwayName;
[loca locb]=ismember(pathway_both,brainout.PathwayName);
heatmaptable=table(pathway_both,lungout.CorrectedAccumulator,brainout.CorrectedAccumulator(locb));
cgo=clustergram(table2array(heatmaptable(1:61,2:3)),'Colormap','redbluecmap','Standardize','column');
set(cgo,'RowLabels',pathway_both(1:61))
labelnew=flip(cgo.RowLabels); %so I can rotate fig in illustrator

%% heatmap brain vs lung
bvslout=readtable('mithriloutputbrm2lm2.txt');
pathway_bl=bvslout.PathwayName;
heatmaptable_bl=table(pathway_bl,bvslout.CorrectedAccumulator,zeros(249,1),bvslout.AdjustedPValue);
cgo=clustergram(table2array(heatmaptable_bl(1:61,2:3)),'Colormap','redbluecmap','Standardize','column');
set(cgo,'RowLabels',pathway_bl(1:61))
labelnewbl=flip(cgo.RowLabels); %so I can rotate fig in illustrator

