% Analyze Seahorse Data Mito Stress Assay

InPath = '/Volumes/xavierlab/Deepti/MetastasisMetabolismPaper/MatlabFig_files';
cd(InPath)
%% make table
time=[1.31
7.76
14.22
20.77
27.24
33.70
40.28
46.75
53.23
59.80
66.28
72.75];
time=repmat(time,8,1);
time=sortrows(time);
replicate=[1:8]';
replicate=repmat(replicate,12,1);

parentalO2=[];
brm2O2=[];
lm2O2=[];
data=xlsread('231_mito.xlsx','Normalized Rate (Columns)');
for i=1:12
    parentalO2=[parentalO2;data(i,29:36)'];
    brm2O2=[brm2O2;data(i,38:45)'];
    lm2O2=[lm2O2;data(i,47:54)'];
end

treatment=repmat([0:3]',24,1);
treatment=sort(treatment);
treatment=categorical(treatment);

O2data = table(time,treatment,replicate,parentalO2,brm2O2,lm2O2);
O2stacked=stack(O2data,{'parentalO2','brm2O2','lm2O2'},...
    'NewDataVariableName','O2',...
    'IndexVariableName','cellline');

%% pca

O2array=table2array(O2stacked(:,5));
k=1;
for i=1:8
    ind=[k:k+2,k+24:k+26,k+48:k+50,k+72:k+74,k+96:k+98,k+120:k+122,k+144:k+146,k+168:k+170,k+192:k+194,k+216:k+218,k+240:k+242,k+264:k+266];
    O2pca(i,:)=O2array(ind);
    k=k+3;
end
[coeff,score,latent,tsquared] = pca(O2pca');

color=[];
color=[0,0,1;1,0,1;0,1,0];
color=repmat(color,12,1);
figure,scatter(score(1:9,1),score(1:9,2),50,color(1:9,:),'o')
hold on
scatter(score(10:18,1),score(10:18,2),50,color(10:18,:),'d')
scatter(score(19:27,1),score(19:27,2),50,color(19:27,:),'*')
scatter(score(28:36,1),score(28:36,2),50,color(28:36,:),'+')
hold off
%%
O2datamat=table2array(O2data(:,4:6));

k=1;
for i=1:32
    if rem(k,8)==1
        k=i*3-2;
    end
    ind=[k,k+8,k+16];
    averageO2m(i,:)=mean(O2datamat(ind,:),1);
    stdevO2m(i,:)=std(O2datamat(ind,:),1);
    k=k+1;
end
nonmitom=mean(averageO2m(25:32,:),1);
nonmitodevm=std(averageO2m(25:32,:),1);
basalrespm=mean(averageO2m(1:8,:)-averageO2m(25:32,:),1);
basalrespdevm=std(averageO2m(1:8,:)-averageO2m(25:32,:),1);
protonleakm=mean(averageO2m(9:16,:)-averageO2m(25:32,:),1);
protonleakdevm=std(averageO2m(9:16,:)-averageO2m(25:32,:),1);
ATPm=mean(averageO2m(1:8,:)-averageO2m(9:16,:),1);
ATPdevm=std(averageO2m(1:8,:)-averageO2m(9:16,:),1);
maxrespm=mean(averageO2m(17:24,:)-averageO2m(25:32,:));
maxrespdevm=std(averageO2m(17:24,:)-averageO2m(25:32,:));
sparerespm=mean(averageO2m(17:24,:)-averageO2m(1:8,:),1);
sparerespdevm=std(averageO2m(17:24,:)-averageO2m(1:8,:),1);

x=[1:6];
datam=[basalrespm,ATPm];
errm=[basalrespdevm,ATPdevm,];
color=[];
color=[0,0,.3;0,0,.6;0,0,1];
color=repmat(color,2,1);
figure,b=bar(x,datam)
b.FaceColor = 'flat';
b.CData(:,:) = color;
hold on
er = errorbar(x,datam,errm,errm);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
hold off

treatment=repmat({'nonmito','basal','leak','atp','max','spare'}',8,1);
treatment=sort(treatment);
replicate=[1:8]';
replicate=repmat(replicate,6,1);
O2results=table(treatment,replicate,[averageO2m(1:8,:)-averageO2m(9:16,:);averageO2m(1:8,:)-averageO2m(25:32,:);averageO2m(9:16,:)-averageO2m(25:32,:);averageO2m(17:24,:)-averageO2m(25:32,:);averageO2m(25:32,:);averageO2m(17:24,:)-averageO2m(1:8,:)]);
O2results = splitvars(O2results,'Var3','NewVariableNames',{'parental','brm2','lm2'});
O2stackedresults=stack(O2results,{'parental','brm2','lm2'},...
    'NewDataVariableName','O2',...
    'IndexVariableName','cellline');
mdl=fitlme(O2stackedresults,'O2~treatment+cellline+treatment:cellline');
mdlatp=fitlme(O2stackedresults(1:24,:),'O2~cellline');
mdlbasal=fitlme(O2stackedresults(25:48,:),'O2~cellline');
mdlleak=fitlme(O2stackedresults(49:72,:),'O2~cellline');
mdlmax=fitlme(O2stackedresults(73:96,:),'O2~cellline');
mdlnonmito=fitlme(O2stackedresults(97:120,:),'O2~cellline');
mdlspare=fitlme(O2stackedresults(121:144,:),'O2~cellline');

%Fold changes
FCbasal=[(mdlbasal.Coefficients.Estimate(2)+mdlbasal.Coefficients.Estimate(1))/mdlbasal.Coefficients.Estimate(1);(mdlbasal.Coefficients.Estimate(3)+mdlbasal.Coefficients.Estimate(1))/mdlbasal.Coefficients.Estimate(1);(mdlbasal.Coefficients.Estimate(3)+mdlbasal.Coefficients.Estimate(1))/(mdlbasal.Coefficients.Estimate(2)+mdlbasal.Coefficients.Estimate(1))];
FCatp=[(mdlatp.Coefficients.Estimate(2)+mdlatp.Coefficients.Estimate(1))/mdlatp.Coefficients.Estimate(1);(mdlatp.Coefficients.Estimate(3)+mdlatp.Coefficients.Estimate(1))/mdlatp.Coefficients.Estimate(1);(mdlatp.Coefficients.Estimate(3)+mdlatp.Coefficients.Estimate(1))/(mdlatp.Coefficients.Estimate(2)+mdlatp.Coefficients.Estimate(1))];
[h,pbbasal]=ttest2(O2results.parental(9:16),O2results.brm2(9:16));
[h,plbasal]=ttest2(O2results.lm2(9:16),O2results.parental(9:16));
[h,blbasal]=ttest2(O2results.lm2(9:16),O2results.brm2(9:16));
[h,pbatp]=ttest2(O2results.parental(1:8),O2results.brm2(1:8));
[h,platp]=ttest2(O2results.lm2(1:8),O2results.parental(1:8));
[h,blatp]=ttest2(O2results.lm2(1:8),O2results.brm2(1:8));
%% Do the same with ECAR data
%make table
parentalECAR=[];
brm2ECAR=[];
lm2ECAR=[];
for i=18:29
    parentalECAR=[parentalECAR;data(i,29:36)'];
    brm2ECAR=[brm2ECAR;data(i,38:45)'];
    lm2ECAR=[lm2ECAR;data(i,47:54)'];
end

treatment=repmat([0:3]',24,1);
treatment=sort(treatment);
treatment=categorical(treatment);
replicate=[1:8]';
replicate=repmat(replicate,12,1);
ECARdata = table(time,treatment,replicate,parentalECAR,brm2ECAR,lm2ECAR);
ECARstacked=stack(ECARdata,{'parentalECAR','brm2ECAR','lm2ECAR'},...
    'NewDataVariableName','ECAR',...
    'IndexVariableName','cellline');

% compare groups
ECARdatamat=table2array(ECARdata(:,4:6));
k=1;
for i=1:32
    if rem(k,8)==1
        k=i*3-2;
    end
    ind=[k,k+8,k+16];
    averageECARm(i,:)=mean(ECARdatamat(ind,:),1);
    stdevECARm(i,:)=std(ECARdatamat(ind,:),1);
    k=k+1;
end
ecar0=mean(averageECARm(1:8,:),1);
ecar0dev=std(averageECARm(1:8,:),1);
ecar1=mean(averageECARm(9:16,:),1);
ecar1dev=std(averageECARm(9:16,:),1);
ecar2=mean(averageECARm(17:24,:),1);
ecar2dev=std(averageECARm(17:24,:),1);
ecar3=mean(averageECARm(25:32,:),1);
ecar3dev=std(averageECARm(25:32,:),1);


x=[1:6];
datam=[ecar0,ecar3];
errm=[ecar0dev,ecar3dev];
color=[];
color=[0,0,.3;0,0,.6;0,0,1];
color=repmat(color,2,1);
figure,b=bar(x,datam)
b.FaceColor = 'flat';
b.CData(:,:) = color;
hold on
er = errorbar(x,datam,errm,errm);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
hold off

treatment=repmat([0:3]',8,1);
treatment=sort(treatment);
treatment=categorical(treatment);
replicate=[1:8]';
replicate=repmat(replicate,4,1);
ECARresults=table(treatment,replicate,[averageECARm(1:8,:);averageECARm(9:16,:);averageECARm(17:24,:);averageECARm(25:32,:)]);
ECARresults = splitvars(ECARresults,'Var3','NewVariableNames',{'parental','brm2','lm2'});
ECARstackedresults=stack(ECARresults,{'parental','brm2','lm2'},...
    'NewDataVariableName','ECAR',...
    'IndexVariableName','cellline');
mdl=fitlme(ECARstackedresults,'ECAR~treatment + cellline+treatment:cellline');
mdl0=fitlme(ECARstackedresults(1:24,:),'ECAR~cellline');
mdl1=fitlme(ECARstackedresults(25:48,:),'ECAR~cellline');
mdl2=fitlme(ECARstackedresults(49:72,:),'ECAR~cellline');
mdl3=fitlme(ECARstackedresults(73:96,:),'ECAR~cellline');

%Fold changes
FC0=[(mdl0.Coefficients.Estimate(2)+mdl0.Coefficients.Estimate(1))/mdl0.Coefficients.Estimate(1);(mdl0.Coefficients.Estimate(3)+mdl0.Coefficients.Estimate(1))/mdl0.Coefficients.Estimate(1);(mdl0.Coefficients.Estimate(3)+mdl0.Coefficients.Estimate(1))/(mdl0.Coefficients.Estimate(2)+mdl0.Coefficients.Estimate(1))];
FC3=[(mdl3.Coefficients.Estimate(2)+mdl3.Coefficients.Estimate(1))/mdl3.Coefficients.Estimate(1);(mdl3.Coefficients.Estimate(3)+mdl3.Coefficients.Estimate(1))/mdl3.Coefficients.Estimate(1);(mdl3.Coefficients.Estimate(3)+mdl3.Coefficients.Estimate(1))/(mdl3.Coefficients.Estimate(2)+mdl3.Coefficients.Estimate(1))];
[h,pb0]=ttest2(ECARresults.parental(1:8),ECARresults.brm2(1:8));
[h,pl0]=ttest2(ECARresults.lm2(1:8),ECARresults.parental(1:8));
[h,bl0]=ttest2(ECARresults.lm2(1:8),ECARresults.brm2(1:8));
[h,pb3]=ttest2(ECARresults.parental(25:32),ECARresults.brm2(25:32));
[h,pl3]=ttest2(ECARresults.lm2(25:32),ECARresults.parental(25:32));
[h,bl3]=ttest2(ECARresults.lm2(25:32),ECARresults.brm2(25:32));