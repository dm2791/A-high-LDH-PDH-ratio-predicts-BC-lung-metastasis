% analyze project GENIE breast cancer metastasis data from AACR

%% load data
alldata=readtable('mbc_genie_2020_clinical_data.csv');
metastatic=strcmp(alldata.SampleType,'Metastatic Recurrence');
metasdata=alldata(metastatic,:);
%% where do cancers metastasize?
tblTopMets = groupcounts(metasdata, 'SiteOfSampleTested');
tblTopMets = sortrows(tblTopMets, 'GroupCount', 'descend');
%combine pleural effusions
tblTopMets.GroupCount(8)=11;
tblTopMets(19,:)=[];
figure
bar(tblTopMets.GroupCount)
xlabel('Metastasis site')
ylabel('Cases')
h = gca;
h.XTick = 1:(height(tblTopMets));
h.XTickLabel = tblTopMets.SiteOfSampleTested;
h.XTickLabelRotation = 90;
title('BRCA')
%% top mutations and copy number changes
mutations=readtable('Mutated_Genes.txt');
mutations=sortrows(mutations,'x_','descend');

mutations.Freq=extractBefore(mutations.Freq,"%");
mutations.Freq=cellfun(@str2num,mutations.Freq);
cna.Freq=extractBefore(cna.Freq,"%");
cna.Freq=cellfun(@str2num,cna.Freq);

% bar plot with top most frequent mutations
ngenes = 50;
figure
bar(mutations.Freq(1:ngenes))
xlabel('Gene mutated')
ylabel('%samples with one or more mutations')
h = gca;
h.XTick = 1:ngenes;
h.XTickLabel = mutations.Gene(1:ngenes);
h.XTickLabelRotation = 90;
title('BRCA')
%% mutation vs met site
%order: liver, LN, brain,bone,lung
%data from cbioportal interface: 
pik3ca=[14 8 0 13 4];
tp53=[20 24 17 3 9];
both=[4 5 10 3 9];
neither=[ 14 6 3 4 0]; % approximate, didn't include all possible mutations (just top 5 after pik3ca and p53)

y = [pik3ca; tp53; both; neither];
y=y(:,[4,1,2,5,3]); %change order to be same as MBCP graph
figure,bar(y,'stacked')
legend({'bone','liver','regional lymph node','lung','brain/cns'},'Location','northeast')
xlabel('Mutation Status')
ylabel('Cases')
h = gca;
h.XTickLabel = {'PIK3CA','TP53','BOTH','NEITHER'};
h.XTickLabelRotation = 90;
title('BRCA')
%% hormone status vs met site
hormone=[];
hormone=strcmp(metasdata.SubtypeOfInitialDiagnosis,'TNBC'); %tripleneg
hormone=[hormone,strcmp(metasdata.SubtypeOfInitialDiagnosis,'HER2+/HR-')]; %her2+
hormone=[hormone,strcmp(metasdata.SubtypeOfInitialDiagnosis,'HER2-/HR+')]; %er/pr+
hormone= [hormone,strcmp(metasdata.SubtypeOfInitialDiagnosis,'HER2+/HR+')];%er/pr/her2+

tripleneg=groupcounts(metasdata(hormone(:,1),:),'SiteOfSampleTested');
her2=groupcounts(metasdata(hormone(:,2),:),'SiteOfSampleTested');
erpr=groupcounts(metasdata(hormone(:,3),:),'SiteOfSampleTested');
erprher2=groupcounts(metasdata(hormone(:,3),:),'SiteOfSampleTested');

%data from cbioportal interface:
triplenegprop=tripleneg.GroupCount([1 4 6 5 2]); %bone, liver, rln, lung, brain
her2prop=her2.GroupCount([1 4 6 5 2]);
erprprop=erpr.GroupCount([1 4 6 5 2]);
erprher2prop=erprher2.GroupCount([1 4 6 5 2]);

y = [triplenegprop'; her2prop'; erprprop'; erprher2prop'];
figure,bar(y,'stacked')
legend({'bone','liver','regional lymph node','lung','brain/cns'},'Location','northeast')
xlabel('Hormone Receptor Status')
ylabel('Cases')
h = gca;
h.XTickLabel = {'TripleNegative','Her2+','ER/PR+','ER/PR/Her2+'};
h.XTickLabelRotation = 90;
title('BRCA')
%% feature selection
stage=zeros(height(metasdata),1);
stage(strcmp(metasdata.StageAtInitialBreastCancerDiagnosis,'I'))=1;
stage(strcmp(metasdata.StageAtInitialBreastCancerDiagnosis,'II'))=2;
stage(strcmp(metasdata.StageAtInitialBreastCancerDiagnosis,'III'))=3;
stage(strcmp(metasdata.StageAtInitialBreastCancerDiagnosis,'IV'))=4;
stage=table(stage,'VariableNames',{'stage'});
predictiontable=[metasdata(:,18),stage,metasdata(:,[32,12,15,5,6,19,20])];
predictiontable.stage=categorical(predictiontable.stage);

responsevar=metasdata.SiteOfSampleTested;
%% model
predictiontable.response=responsevar;
[idx score] = fscchi2(predictiontable,"response");
pvals=exp(-score); %stage,subtype and traztuzumab significant
figure,bar(score(idx))
set(gca,'xtick',1:9,'xticklabel',predictiontable.Properties.VariableNames(idx));
xtickangle(90)
xlabel('Predictor rank')
ylabel('Predictor importance score')

figure,bar(log(pvals(idx)))
hold on
plot([0 10],[log(.05) log(.05)],'LineStyle','--','Color','k')
set(gca,'xtick',1:9,'xticklabel',predictiontable.Properties.VariableNames(idx));
xtickangle(90)
xlabel('Predictor rank')
ylabel('pvalue')
%% regression lung
predictiontable.response=responsevar;
predictiontable.lung=zeros(height(predictiontable),1);
predictiontable.lung(strcmp(predictiontable.response,'Lung'))=1;
predictiontable.lung=categorical(predictiontable.lung);
predictiontable.response=[];
mdl=fitglm(predictiontable); 
figure, hold on
 scatter(log(mdl.Coefficients.pValue([3:10,26])),[1:9],'c','filled')
%% regression brain
predictiontable.response=responsevar;
predictiontable.brain=zeros(height(predictiontable),1);
predictiontable.brain(strcmp(predictiontable.response,'Brain/CNS'))=1;
predictiontable.brain=categorical(predictiontable.brain);
predictiontable.response=[];
predictiontable.lung=[];
mdl=fitglm(predictiontable); 

 scatter(log(mdl.Coefficients.pValue([3:10,26])),[1:9],'m','filled')
%% regression liver
predictiontable.response=responsevar;
predictiontable.liver=zeros(height(predictiontable),1);
predictiontable.liver(strcmp(predictiontable.response,'Liver'))=1;
predictiontable.liver=categorical(predictiontable.liver);
predictiontable.response=[];
predictiontable.brain=[];
mdl=fitglm(predictiontable); 
 scatter(log(mdl.Coefficients.pValue([3:10,26])),[1:9],'k','filled')
%% regression bone
predictiontable.response=responsevar;
predictiontable.bone=zeros(height(predictiontable),1);
predictiontable.bone(strcmp(predictiontable.response,'Bone'))=1;
predictiontable.bone=categorical(predictiontable.bone);
predictiontable.response=[];
predictiontable.liver=[];
mdl=fitglm(predictiontable); 
 scatter(log(mdl.Coefficients.pValue([3:10,26])),[1:9],'g','filled')
%% regression lymph node
predictiontable.response=responsevar;
predictiontable.ln=zeros(height(predictiontable),1);
predictiontable.ln(strcmp(predictiontable.response,'Lymph Node'))=1;
predictiontable.ln=categorical(predictiontable.ln);
predictiontable.response=[];
predictiontable.bone=[];
mdl=fitglm(predictiontable); 
 scatter(log(mdl.Coefficients.pValue([3:10,26])),[1:9],'b','filled')
 %%
 plot([log(.05),log(.05)],[0,10],'k','LineStyle','--')
 hold off
 title('Predictors p val')
ylim([0 10])
yticks(1:9)
set(gca,'YTickLabel',mdl.CoefficientNames([3:10,26]))
xlabel('logp')
legend('lung', 'brain/cns' ,'liver' ,'bone', 'lymph node', 'p=.05','Location','bestoutside')