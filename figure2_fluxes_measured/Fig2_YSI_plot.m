% YSI calculations
%% load data
YSIdata=readmatrix('YSIdatamat.csv');
glucosedata=YSIdata(:,1);
lactatedata=YSIdata(:,2);
glutaminedata=YSIdata(:,3);
glutamatedata=YSIdata(:,4);

%% fold changes
FCB=mean(YSIdata(4:6,:))./mean(YSIdata(1:3,:));
FCL=mean(YSIdata(7:9,:))./mean(YSIdata(1:3,:));
[h,pb]=ttest2(YSIdata(1:3,:),YSIdata(4:6,:));
[h,pl]=ttest2(YSIdata(1:3,:),YSIdata(7:9,:));
[h,bl]=ttest2(YSIdata(7:9,:),YSIdata(4:6,:));

%% plot
figure,bar([mean(glucosedata(1:3)),mean(glucosedata(4:6)),mean(glucosedata(5:9))])
figure,bar([mean(lactatedata(1:3)),mean(lactatedata(4:6)),mean(lactatedata(5:9))])
figure,bar([mean(glutaminedata(1:3)),mean(glutaminedata(4:6)),mean(glutaminedata(5:9))])
figure,bar([mean(glutamatedata(1:3)),mean(glutamatedata(4:6)),mean(glutamatedata(5:9))])
