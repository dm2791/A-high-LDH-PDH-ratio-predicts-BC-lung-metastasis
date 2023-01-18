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

%% select patients
patients=metdata_rna.Properties.VariableNames;
patients=erase(patients,'MBC_');
patients=extractBefore(patients,'_Tumor');

patient_select=ismember(metdata_clinical.x_PatientIdentifier,patients);
patient_select=metdata_clinical(patient_select,:);

%% make table
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
%% plot
lung=strcmp(combined_data(:,3),'YES');
brain=strcmp(combined_data(:,4),'YES');
none=strcmp(combined_data(:,3),'N/A');
other2=strcmp(combined_data(:,3),'NO') & strcmp(combined_data(:,4),'NO'); %not lung or brain

lungdata=[ldhb(lung),ldha(lung),pdha1(lung),pdha2(lung)];
braindata=[ldhb(brain),ldha(brain),pdha1(brain),pdha2(brain)];
nonedata=[ldhb(none),ldha(none),pdha1(none),pdha2(none)];
otherdata2=[ldhb(other2),ldha(other2),pdha1(other2),pdha2(other2)];
lung_ratio=(lungdata(:,1)+lungdata(:,2))./(lungdata(:,3)+lungdata(:,4));
brain_ratio=(braindata(:,1)+braindata(:,2))./(braindata(:,3)+braindata(:,4));
none_ratio=(nonedata(:,1)+nonedata(:,2))./(nonedata(:,3)+nonedata(:,4));
other_ratio2=(otherdata2(:,1)+otherdata2(:,2))./(otherdata2(:,3)+otherdata2(:,4));

grp = [ones(1,length(other_ratio2)),ones(1,length(brain_ratio))*2,ones(1,length(lung_ratio))*3];
figure, boxplot([other_ratio2;brain_ratio;lung_ratio],grp)
