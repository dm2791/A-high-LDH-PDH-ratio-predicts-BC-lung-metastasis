% Analyze Seahorse Data Glycolytic Rate assay

%% make tables
time=[1.3
7.8
14.2
20.8
27.3
33.7
40.3
46.7
53.2
59.7
66.2];
time=repmat(time,8,1);
time=sortrows(time);
replicate=[1:8]';
replicate=repmat(replicate,11,1);

parentalO2=[];
brm2O2=[];
lm2O2=[];
data=readtable('Deepti_231 test_GlycolyticRate_050219_norm.xlsx','Sheet','Normalized Rate (Columns)');
for i=9:19
    parentalO2=[parentalO2;cellfun(@str2num,table2array(data(i,2:9)))'];
    brm2O2=[brm2O2;cellfun(@str2num,table2array(data(i,11:18)))'];
    lm2O2=[lm2O2;cellfun(@str2num,table2array(data(i,20:27)))'];
end
treatment=repmat([0:1]',24,1);
treatment=[treatment;repmat([2],40,1)];
treatment=sort(treatment);
treatment=categorical(treatment);

O2data = table(time,treatment,replicate,parentalO2,brm2O2,lm2O2);
O2stacked=stack(O2data,{'parentalO2','brm2O2','lm2O2'},...
    'NewDataVariableName','O2',...
    'IndexVariableName','cellline');

parentalPER=[];
brm2PER=[];
lm2PER=[];
for i=41:51
    parentalPER=[parentalPER;cellfun(@str2num,table2array(data(i,2:9)))'];
    brm2PER=[brm2PER;cellfun(@str2num,table2array(data(i,11:18)))'];
    lm2PER=[lm2PER;cellfun(@str2num,table2array(data(i,20:27)))'];
end
PERdata = table(time,treatment,replicate,parentalPER,brm2PER,lm2PER);
PERstacked=stack(PERdata,{'parentalPER','brm2PER','lm2PER'},...
    'NewDataVariableName','PER',...
    'IndexVariableName','cellline');

%% analyze
O2datamat=table2array(O2data(:,4:6));
k=25;
for i=1:8
        ind=[k,k+8,k+16];
    averageOCR_rotAA(i,:)=mean(O2datamat(ind,:),1);
    k=k+1;
end


for i=1:24
    ind=rem(i,8);
    if ind==0
    ind=8;
    end
    mitoOCR(i,:)=O2datamat(i,:)-averageOCR_rotAA(ind,:);
end
mitoPER=mitoOCR*.6;

PERdatamat=table2array(PERdata(:,4:6));
basalglycolysis=mean(PERdatamat(1:24,:)-mitoPER(:,:),1);
basalglycolysisdev=std(PERdatamat(1:24,:)-mitoPER(:,:),1);
compensatoryglycolysis2=mean(PERdatamat(25:48,:),1);
compensatoryglycolysis2dev=std(PERdatamat(25:48,:),1);

x=[1:6];
datam=[basalglycolysis,compensatoryglycolysis2];
errm=[basalglycolysisdev,compensatoryglycolysis2dev];
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

treatment=repmat({'basal','compensatory2'}',24,1);
treatment=sort(treatment);
replicate=[1:24]';
replicate=repmat(replicate,2,1);
glycresults=table(treatment,replicate,[PERdatamat(1:24,:)-mitoPER(:,:);PERdatamat(25:48,:)]);
glycresults = splitvars(glycresults,'Var3','NewVariableNames',{'parental','brm2','lm2'});
glycstackedresults=stack(glycresults,{'parental','brm2','lm2'},...
    'NewDataVariableName','glyc',...
    'IndexVariableName','cellline');
mdl=fitlm(glycstackedresults,'glyc~treatment+cellline+treatment:cellline');
mdlbasal=fitlme(glycstackedresults(1:72,:),'glyc~cellline');
mdlcompensatory=fitlme(glycstackedresults(73:144,:),'glyc~cellline');

%% fold changes
%b/p, l/p, l/b
FCbasal=[(mdlbasal.Coefficients.Estimate(2)+mdlbasal.Coefficients.Estimate(1))/mdlbasal.Coefficients.Estimate(1);(mdlbasal.Coefficients.Estimate(3)+mdlbasal.Coefficients.Estimate(1))/mdlbasal.Coefficients.Estimate(1);(mdlbasal.Coefficients.Estimate(3)+mdlbasal.Coefficients.Estimate(1))/(mdlbasal.Coefficients.Estimate(2)+mdlbasal.Coefficients.Estimate(1))];
FCcompensatory=[(mdlcompensatory.Coefficients.Estimate(2)+mdlcompensatory.Coefficients.Estimate(1))/mdlcompensatory.Coefficients.Estimate(1);(mdlcompensatory.Coefficients.Estimate(3)+mdlcompensatory.Coefficients.Estimate(1))/mdlcompensatory.Coefficients.Estimate(1);(mdlcompensatory.Coefficients.Estimate(3)+mdlcompensatory.Coefficients.Estimate(1))/(mdlcompensatory.Coefficients.Estimate(2)+mdlcompensatory.Coefficients.Estimate(1))];

%ttests
[h,pbbasal]=ttest2(glycresults.parental(1:24),glycresults.brm2(1:24));
[h,plbasal]=ttest2(glycresults.parental(1:24),glycresults.lm2(1:24));
[h,blbasal]=ttest2(glycresults.lm2(1:24),glycresults.brm2(1:24));

[h,pbcompens]=ttest2(glycresults.parental(25:48),glycresults.brm2(25:48));
[h,plcompens]=ttest2(glycresults.parental(25:48),glycresults.lm2(25:48));
[h,blcompens]=ttest2(glycresults.lm2(25:48),glycresults.brm2(25:48));