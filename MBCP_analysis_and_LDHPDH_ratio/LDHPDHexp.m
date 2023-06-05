% analyze RNA expression data 231 cells

%% read files
ldhbraindata=readmatrix('ldhbrain.csv');
ldhlungdata=readmatrix('ldhlung.csv');
pdhbraindata=readmatrix('pdhbrain.csv');
pdhlungdata=readmatrix('pdhlung.csv');
%% 
%FC
brainfc=mean(ldhbraindata(1:2,4:5),2)./mean(ldhbraindata(1:2,1:3),2);
lungfc=mean(ldhlungdata(1:2,3:6),2)./mean(ldhlungdata(1:2,1:2),2);
[h,pb]=ttest2(ldhbraindata(1:2,4:5)',ldhbraindata(1:2,1:3)');
[h,pl]=ttest2(ldhlungdata(1:2,3:6)',ldhlungdata(1:2,1:2)');
labels={'LDHA';'LDHB'};
labels=categorical(labels);
figure
bar(labels,[mean(ldhbraindata(1:2,1:3),2),mean(ldhbraindata(1:2,4:5),2)])
figure
bar(labels,[mean(ldhlungdata(1:2,1:2),2),mean(ldhlungdata(1:2,3:6),2)])

%% ratio of ldhb to pdh
%FC
brainfc_l=mean(ldhbraindata(:,4:5),2)./mean(ldhbraindata(:,1:3),2);
lungfc_l=mean(ldhlungdata(:,3:6),2)./mean(ldhlungdata(:,1:2),2);
brainfc_p=mean(pdhbraindata(:,4:5),2)./mean(pdhbraindata(:,1:3),2);
lungfc_p=mean(pdhlungdata(:,3:6),2)./mean(pdhlungdata(:,1:2),2);

%FC ldh/pdh
brainfcsum=sum(mean(ldhbraindata(:,4:5),2))/sum(mean(pdhbraindata(:,4:5),2));
lungfcsum=sum(mean(ldhlungdata(:,3:6),2))/sum(mean(pdhlungdata(:,3:6),2));
parentalfcbsum=sum(mean(ldhbraindata(:,1:3),2))/sum(mean(pdhbraindata(:,1:3),2));
parentalfclsum=sum(mean(ldhlungdata(:,1:2),2))/sum(mean(pdhlungdata(:,1:2),2));

brainfc=sum(ldhbraindata(:,4:5),1)./sum(pdhbraindata(:,4:5),1);
lungfc=sum(ldhlungdata(:,3:6),1)./sum(pdhlungdata(:,3:6),1);
parentalfcb=sum(ldhbraindata(:,1:3),1)./sum(pdhbraindata(:,1:3),1);
parentalfcl=sum(ldhlungdata(:,1:2),1)./sum(pdhlungdata(:,1:2),1);
[h,pbl]=ttest2(brainfc,lungfc);
[h,ppb]=ttest2(brainfc,parentalfcb);
[h,ppl]=ttest2(parentalfcl,lungfc);
grp = [zeros(1,length(parentalfcb)),ones(1,length(brainfc)),2*ones(1,length(lungfc))];
figure, boxplot([parentalfcb,brainfc,lungfc],grp)

%% analyze MBCP LDH and PDH levels
metdata_clinical=readtable('data_clinical_patient.txt');
metdata_rna=readtable('data_RNA_Seq_v2_expression_median.txt');

patients=metdata_rna.Properties.VariableNames;
patients=erase(patients,'MBC_');
patients=extractBefore(patients,'_Tumor');

metdata_clinical.SampleID=strrep(metdata_clinical.SampleID,'-','_');
patient_select=ismember(metdata_clinical.SampleID,metdata_rna.Properties.VariableNames);
patient_select=metdata_clinical(patient_select,:);

%% make table
combined_data=[patient_select.SampleID,patient_select.MedREverMetastaticSites,patient_select.MedREverLungMets,patient_select.MedREverBrain_CNSMets,patient_select.PATHProcedureLocation,num2cell(patient_select.PATHSampleCollectionTime)];
br=strcmp(combined_data(:,5),'BREAST') | strcmp(combined_data(:,5),'CHEST WALL') ;
combined_data_keep=combined_data(br,:);
%% restrict to earliest timepoint
patients_keep=extractBefore(combined_data_keep(:,1),'SM');
patients_keep=unique(patients_keep);
sel_min=[];
for i=1:length(patients_keep)
    findpat=strfind(combined_data_keep(:,1),patients_keep(i));
    findpat=~cellfun(@isempty,findpat);
    sel=cell2mat(combined_data_keep(findpat,6));
    add=zeros(length(sel),1);
    [m k]=min(sel);
    add(k)=1;
    sel_min=[sel_min;add];
end
combined_data_keep=combined_data_keep(logical(sel_min),:);

lldhb=[];
ldha=[];
pdha1=[];
pdha2=[];
pdhb=[];

for i=1:length(combined_data_keep)
    ind=find(strcmp(metdata_rna.Properties.VariableNames,combined_data_keep(i,1)));
    ldhb(i)=table2array(metdata_rna(4124,ind(1)));
    ldha(i)=table2array(metdata_rna(7217,ind(1)));
    pdha1(i)=table2array(metdata_rna(6812,ind(1)));
    pdha2(i)=table2array(metdata_rna(11518,ind(1)));
    pdhb(i)=table2array(metdata_rna(13008,ind(1)));
end
ldhb=ldhb';
ldha=ldha';
pdha1=pdha1';
pdha2=pdha2';
pdhb=pdhb';
ratio=(ldhb+ldha)./(pdha1+pdha2+pdhb);
%% plot
lung=strcmp(combined_data(:,3),'YES');
brain=strcmp(combined_data(:,4),'YES');
other2=strcmp(combined_data(:,3),'NO') & strcmp(combined_data(:,4),'NO'); %not lung or brain

lungdata=[ldhb(lung),ldha(lung),pdha1(lung),pdha2(lung)];
braindata=[ldhb(brain),ldha(brain),pdha1(brain),pdha2(brain)];
otherdata2=[ldhb(other2),ldha(other2),pdha1(other2),pdha2(other2)];
lung_ratio=(lungdata(:,1)+lungdata(:,2))./(lungdata(:,3)+lungdata(:,4));
brain_ratio=(braindata(:,1)+braindata(:,2))./(braindata(:,3)+braindata(:,4));
other_ratio2=(otherdata2(:,1)+otherdata2(:,2))./(otherdata2(:,3)+otherdata2(:,4));

grp = [ones(1,length(other_ratio2)),ones(1,length(brain_ratio))*2,ones(1,length(lung_ratio))*3];
figure, boxplot([other_ratio2;brain_ratio;lung_ratio],grp)
%scatterplot
lungsum=[(lungdata(:,1)+lungdata(:,2)),(lungdata(:,3)+lungdata(:,4)+lungdata(:,5))];
brainsum=[(braindata(:,1)+braindata(:,2)),(braindata(:,3)+braindata(:,4)+braindata(:,5))];
other2sum=[(otherdata2(:,1)+otherdata2(:,2)),(otherdata2(:,3)+otherdata2(:,4)+otherdata2(:,5))];
figure, gscatter([other2sum(:,1);brainsum(:,1);lungsum(:,1)],[other2sum(:,2);brainsum(:,2);lungsum(:,2)],grp)

%% make archetypes: nnmf
sel=ismember(metdata_rna.Properties.VariableNames,combined_data_keep(:,1));
mExpression=table2array(metdata_rna(:,sel));
[W2,H2] = nnmf(mExpression',6);
[M,arch]=max(W2,[],2);
archetype=zeros(99,6);
for i=1:99
    archetype(i,arch(i))=1;
end
%%
[B,I] = sort(ratio);
[coeff,score,latent,tsquared,explained,mu] = pca(mExpression');
[score, umap, clusterIdentifiers] = run_umap(mExpression');

figure
h = gscatter(score(:, 1), score(:, 2), archetype);
hold on
scatter(score(find(lung),1),score(find(lung),2),100,'ko'); %lung
hold off
legend({'arch1','arch2','arch3','arch4','arch5','arch6','lung ratios'})
xlabel('UMAP1')
ylabel('UMAP2')

%plot color range
figure,scatter(score(:, 1), score(:, 2), [], ratio, 'filled')
set(gca, 'CLim', [.5, 5])
%circle lung
hold on
scatter(score(find(lung),1),score(find(lung),2),100,'ko'); %lung
hold off
colorbar
legend({'LDH/PDH','lung'})
xlabel('UMAP1')
ylabel('UMAP2')

%% feature selection: use features also available in AACR GENIE dataset
%make table
[ia ic]=ismember(combined_data_keep(:,1),metdata_clinical.SampleID);
metdata_clinical_keep=metdata_clinical(ic,:);
mutnum=metdata_clinical_keep.MutationCount;
hist=metdata_clinical_keep.PATHSampleHistology;

stage=metdata_clinical_keep.MedRStageAtDiagnosis;
stage=cellfun (@(x) x(1),stage,'UniformOutput',false);
stage(strcmp(stage,'U')) = {'0'};
stage(strcmp(stage,'N')) = {'0'};
stage=str2num(cell2mat(stage));
stage(stage==0)=NaN;

hormone=[];
hormone=[hormone,strcmp(metdata_clinical_keep.MedRDiagnosticERStatus,'NEGATIVE')&strcmp(metdata_clinical_keep.MedRDiagnosticHER2Status,'NEGATIVE');]; %tripleneg
hormone=[hormone,~strcmp(metdata_clinical_keep.MedRDiagnosticERStatus,'POSITIVE')&strcmp(metdata_clinical_keep.MedRDiagnosticHER2Status,'POSITIVE');]; %her2+
hormone=[hormone,strcmp(metdata_clinical_keep.MedRDiagnosticERStatus,'POSITIVE')&~strcmp(metdata_clinical_keep.MedRDiagnosticHER2Status,'POSITIVE');]; %er/pr+
hormone=[hormone,strcmp(metdata_clinical_keep.MedRDiagnosticERStatus,'POSITIVE')&strcmp(metdata_clinical_keep.MedRDiagnosticHER2Status,'POSITIVE');]; %er/pr/her2+

subtype=[];
for i=1:99
    sub=find(hormone(i,:),1);
    if isempty(sub)
        subtype(i)=NaN;
    else
        subtype(i)=find(hormone(i,:),1);
    end
end

bone=contains(combined_data_keep(:,2),'BONE');
liver=contains(combined_data_keep(:,2),'LIVER');
lung=contains(combined_data_keep(:,2),'LUNG');
brain=contains(combined_data_keep(:,2),'BRAIN');
lymph=contains(combined_data_keep(:,2),'LYMPH NODE');
metmatrix_top5=[bone,liver,lung,brain,lymph];
metmatrix_top5=[metmatrix_top5,sum(metmatrix_top5(:,1:5),2)];
predictors=[mutnum,stage,subtype',arch,ratio];
predictormatrix=[];
histadd=[];
for i=1:99
    numrep=metmatrix_top5(i,6);
    if numrep==0
        numrep=1;
    end
    predictormatrix=[predictormatrix;repmat(predictors(i,:),[numrep,1])];
    histadd=[histadd;repmat(hist(i),[numrep,1])];
end
predictortable=array2table(predictormatrix);
predictortable.Properties.VariableNames={'mutation num' 'stage' 'subtype' 'archetype','ldhpdh'};

predictortable.hist=histadd;
predictortable = movevars(predictortable, 'ldhpdh', 'After', 'hist');

predictortable.subtype=categorical(predictortable.subtype);
predictortable.subtype(predictortable.subtype=='1')={'tripleneg'};
predictortable.subtype(predictortable.subtype=='2')={'her2'};
predictortable.subtype(predictortable.subtype=='3')={'erpr'};
predictortable.subtype(predictortable.subtype=='4')={'erprher2'};
predictortable.archetype=categorical(predictortable.archetype);
predictortable.stage=categorical(predictortable.stage);
%%
responsevar=[];
for i=1:99
    response=find(metmatrix_top5(i,1:5));
    if isempty(response)
        responsevar=[responsevar;7];
    else
        responsevar=[responsevar;response'];
    end
end

responsevar_names={};
responsevar_names(responsevar==1)={'bone'};
responsevar_names(responsevar==2)={'liver'};
responsevar_names(responsevar==3)={'lung'};
responsevar_names(responsevar==4)={'brain'};
responsevar_names(responsevar==5)={'lymph'};
responsevar_names(responsevar==7)={'other'};
predictortable.response=responsevar_names';

%% add ldh and pdh individually
ldhaval=[];
ldhbval=[];
pdha1val=[];
pdha2val=[];
pdhbval=[];
for i=1:99
    numrep=metmatrix_top5(i,6);
    if numrep==0
        numrep=1;
    end
    ldhaval=[ldhaval;repmat(ldha(i),[numrep,1])];
    ldhbval=[ldhbval;repmat(ldhb(i),[numrep,1])];
    pdha1val=[pdha1val;repmat(pdha1(i),[numrep,1])];
    pdha2val=[pdha2val;repmat(pdha2(i),[numrep,1])];
    pdhbval=[pdhbval;repmat(pdhb(i),[numrep,1])];
end
predictortable_full=array2table([predictormatrix,ldhaval,ldhbval,pdha1val,pdha2val,pdhbval]);
predictortable_full.Properties.VariableNames={'mutation num' 'stage' 'subtype' 'archetype','ldhpdh','ldha','ldhb','pdha1','pdha2','pdhb'};
predictortable_full.subtype=categorical(predictortable_full.subtype);
predictortable_full.subtype(predictortable_full.subtype=='1')={'tripleneg'};
predictortable_full.subtype = reordercats(predictortable_full.subtype,{'tripleneg','1','2','3','4'});
predictortable_full.archetype=categorical(predictortable_full.archetype);
predictortable_full.stage=categorical(predictortable_full.stage);

predictortable_full.response=responsevar_names';
[idx score] = fscchi2(predictortable_full,"response");
[idx score] = fscchi2(predictortable_full,"response",'NumBins',5);
pvals=exp(-score); 
figure,bar(score(idx))
set(gca,'xtick',1:10,'xticklabel',predictortable_full.Properties.VariableNames(idx));
xtickangle(90)
xlabel('Predictor rank')
ylabel('Predictor importance score')

figure,bar(log(pvals(idx)))
hold on
plot([0 10],[log(.05) log(.05)],'LineStyle','--','Color','k')
set(gca,'xtick',1:10,'xticklabel',predictortable_full.Properties.VariableNames(idx));
xtickangle(90)
xlabel('Predictor rank')
ylabel('logpvalue')

%% regression 
lungrep=[];
brainrep=[];
bonerep=[];
liverrep=[];
lnrep=[];
for i=1:99
    numrep=metmatrix_top5(i,6);
    if numrep==0
        numrep=1;
    end
    if metmatrix_top5(i,3)==1
        lungrep=[lungrep;repmat(1,[numrep,1])];
    else
        lungrep=[lungrep;repmat(0,[numrep,1])];
    end
    if metmatrix_top5(i,4)==1
        brainrep=[brainrep;repmat(1,[numrep,1])];
    else
        brainrep=[brainrep;repmat(0,[numrep,1])];
    end
    if metmatrix_top5(i,1)==1
        bonerep=[bonerep;repmat(1,[numrep,1])];
    else
        bonerep=[bonerep;repmat(0,[numrep,1])];
    end
    if metmatrix_top5(i,2)==1
        liverrep=[liverrep;repmat(1,[numrep,1])];
    else
        liverrep=[liverrep;repmat(0,[numrep,1])];
    end
    if metmatrix_top5(i,5)==1
        lnrep=[lnrep;repmat(1,[numrep,1])];
    else
        lnrep=[lnrep;repmat(0,[numrep,1])];
    end
end
predictortable_full.lung=categorical(lungrep);
predictortable_full.brain=categorical(brainrep);
predictortable_full.bone=categorical(bonerep);
predictortable_full.liver=categorical(liverrep);
predictortable_full.ln=categorical(lnrep);
predictortable_full.response=[];
predictortable_full_keep=unique(predictortable_full);
mdl=fitglm(predictortable_full_keep(:,[3,2,1,5,11]),'Intercept',false);
figure, hold on
scatter(log(mdl.Coefficients.pValue([1:end])),[1:9],'c','filled')
mdl=fitglm(predictortable_full_keep(:,[3,2,1,5,12]),'Intercept',false);
scatter(log(mdl.Coefficients.pValue([1:end])),[1:9],'m','filled')
mdl=fitglm(predictortable_full_keep(:,[3,2,1,5,13]),'Intercept',false);
scatter(log(mdl.Coefficients.pValue([1:end])),[1:9],'g','filled')
mdl=fitglm(predictortable_full_keep(:,[3,2,1,5,14]),'Intercept',false);
scatter(log(mdl.Coefficients.pValue([1:end])),[1:9],'k','filled')
mdl=fitglm(predictortable_full_keep(:,[3,2,1,5,15]),'Intercept',false);
scatter(log(mdl.Coefficients.pValue([1:end])),[1:9],'b','filled')
plot([log(.05),log(.05)],[0,10],'k','LineStyle','--')
hold off
title('Predictors p val')
ylim([0 10])
yticks(1:9)
set(gca,'YTickLabel',mdl.CoefficientNames([1:end]))
xlabel('logp')
legend('lung', 'brain/cns', 'bone','liver','lymph node','p=.05','Location','bestoutside')

%%
[c ia ic]=unique(predictortable_full_keep.pdhb,'stable');
predictortable_unique=predictortable_full_keep(ia,:);
predictortable_unique.lung=logical(double(predictortable_unique.lung)-1);
%% predict
mdl15=fitglm(predictortable_unique(:,[5,11]),'Distribution','binomial');
[X,Y,T,AUC] = perfcurve(predictortable_unique.lung,mdl15.Fitted.Probability,1); %.78
figure, plot(X,Y)

%% binarize
figure, histogram(ratio)
figure,scatter([1:99],sort(ratio)) %2.77 or 3 as inflection
%select lung and equal non-lung at random
lungf=find(predictortable_unique.lung);
lungnotf=find(~predictortable_unique.lung);
X=[];
Y=[];
for j=1:10
sel=randi([1 94],11,1);
predictortable_sel=predictortable_unique([lungf;sel],:);
ratio_bin=predictortable_sel.ldhpdh>2.77;
testresult=[];
for i=1:height(predictortable_sel)
mdl16=fitglm(ratio_bin([1:i-1,i+1:end]),predictortable_sel.lung([1:i-1,i+1:end]),'Distribution','binomial');
testresult(i)=predict(mdl16,ratio_bin(i));
end
index = double(table2array(predictortable_sel(:,11)))==round(testresult)';
X=[X;double(table2array(predictortable_sel(:,11)))];
Y=[Y,round(testresult)];
end
predaccuracy=sum(X==Y')/length(Y); %76%
C = confusionmat(X,Y');
figure,confusionchart(C)

