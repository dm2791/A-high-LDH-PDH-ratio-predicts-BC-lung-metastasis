%% load data
metdata_clinical=readtable('data_clinical_patient.txt');
metdata_rna=readtable('data_RNA_Seq_v2_expression_median.txt');
%% select patients
patients=metdata_rna.Properties.VariableNames;
patients=erase(patients,'MBC_');
patients=extractBefore(patients,'_Tumor');

patient_select=ismember(metdata_clinical.x_PatientIdentifier,patients);
patient_select=metdata_clinical(patient_select,:);
%% make table
patient_select=ismember(metdata_clinical.x_PatientIdentifier,patients);
patient_select=metdata_clinical(patient_select,:);
combined_data=[patient_select.x_PatientIdentifier,patient_select.MedRMetastaticSitesAtMetastaticDiagnosis,patient_select.MedRLungMetsAtMetsDx,patient_select.MedRBrain_CNSMetsAtMetsDx];

ldhb=[];
ldha=[];
pdha1=[];
pdha2=[];
for i=1:113
    ind=find(strcmp(patients,combined_data(i)));
    ldhb(i)=table2array(metdata_rna(4124,ind(1)));
    ldha(i)=table2array(metdata_rna(7217,ind(1)));
    pdha1(i)=table2array(metdata_rna(6812,ind(1)));
    pdha2(i)=table2array(metdata_rna(11518,ind(1)));
end
ldhb=ldhb';
ldha=ldha';
pdha1=pdha1';
pdha2=pdha2';
ratio=(ldhb+ldha)./(pdha1+pdha2);
%%
lung=strcmp(combined_data(:,3),'YES');
brain=strcmp(combined_data(:,4),'YES');
none=strcmp(combined_data(:,2),'N/A');
other2=strcmp(combined_data(:,3),'NO') & strcmp(combined_data(:,4),'NO'); %not lung or brain

lungdata=[ldhb(lung),ldha(lung),pdha1(lung),pdha2(lung)];
braindata=[ldhb(brain),ldha(brain),pdha1(brain),pdha2(brain)];
otherdata2=[ldhb(other2),ldha(other2),pdha1(other2),pdha2(other2)];
lung_ratio=ratio(lung);
brain_ratio=ratio(brain);
other_ratio2=ratio(other2);

[h,p1]=ttest2(lung_ratio,other_ratio2);
[h,p2]=ttest2(other_ratio2,brain_ratio);
grp = [ones(1,length(other_ratio2)),ones(1,length(brain_ratio))*2,ones(1,length(lung_ratio))*3];
figure, boxplot([other_ratio2;brain_ratio;lung_ratio],grp)
%% load data rna
metdata_rna=readtable('data_RNA_Seq_v2_expression_median.txt');
geneNames = metdata_rna.Hugo_Symbol;

patients=metdata_rna.Properties.VariableNames;
patients=erase(patients,'MBC_');
patients=extractBefore(patients,'_Tumor');
RNAseq_participants=patients;
RNAseq_participants=RNAseq_participants(3:end)';
%remove duplicates, choose first of each
[c ia]=unique(RNAseq_participants);
a=sort(ia);
RNAseq_participants=RNAseq_participants(a);
metdata_rna=metdata_rna(:,a+2);

mExpression=table2array(metdata_rna);
%% overlay top ldh/pdh ratios on umap

[B,I] = sort(ratio);
[score, umap, clusterIdentifiers] = run_umap(mExpression');

%plot color range
figure,scatter(score(:, 1), score(:, 2), [], ratio, 'filled')
set(gca, 'CLim', [0, 7])
%circle lung
hold on
scatter(score(find(lung),1),score(find(lung),2),100,'ko'); %lung
hold off