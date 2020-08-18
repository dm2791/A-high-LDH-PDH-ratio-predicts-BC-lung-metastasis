% analyze metastatic breast cancer project
%% load data clinical
metdata_clinical=readtable('data_clinical_patient.txt');
metdata_clinical(1:4,:)=[];
clinicalMetadata=metdata_clinical;
clinicalMetadata.participant = metdata_clinical.x_PatientIdentifier;
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
%% load data mutations
metdata_mutations=readtable('data_mutations_extended.txt');
mutationTable=metdata_mutations(:,[1,10,17]);
patients=mutationTable.Tumor_Sample_Barcode;
patients=erase(patients,'MBC-');
patients=extractBefore(patients,'-Tumor');
mutationTable.participant = patients;

%% unstack the mutation scores
disp('unstacking full mutation matrix')
mutationTable.ones = ones(height(mutationTable), 1);
tic
mutationMatrix = unstack(mutationTable(:, {'participant' 'Hugo_Symbol' 'ones'}), 'ones', 'Hugo_Symbol');
toc
m = mutationMatrix{:, 2:end};
m(isnan(m)) = 0;
mutationMatrix{:, 2:end} = m;
disp('done.');

%% find patients in all three sets, and keep only those
patient = intersect(clinicalMetadata.participant, mutationMatrix.participant);
patient = intersect(patient, RNAseq_participants);
% remove unique from transcriptomics
mExpression(:, ~ismember(RNAseq_participants, patient)) = [];
RNAseq_participants(~ismember(RNAseq_participants, patient)) = [];
% remove unique from clinical metadata
clinicalMetadata(~ismember(clinicalMetadata.participant, patient), :) = [];
% remove unique from genomics
mutationMatrix(~ismember(mutationMatrix.participant, patient), :) = [];

%% reorder patients in clinical to match genomics
[~, patientReorder] = ismember(mutationMatrix.participant, clinicalMetadata.participant);
clinicalMetadata = clinicalMetadata(patientReorder, :);

%reorder patients in rna to match genomics
[~, patientReorder] = ismember(mutationMatrix.participant, RNAseq_participants);
RNAseq_participants=RNAseq_participants(patientReorder);
metdata_rna = metdata_rna(:,patientReorder);
mExpression=table2array(metdata_rna);
%% Where do cancers metastasize?
clinicalMetadata.metSiteExtended = clinicalMetadata.MedRMetastaticSitesAtMetastaticDiagnosis;
idxOther = contains(clinicalMetadata.metSiteExtended, 'N/A', 'IgnoreCase', true);
tblTopMetsThisCancer = groupcounts(clinicalMetadata, 'metSiteExtended');
tblTopMetsThisCancer = sortrows(tblTopMetsThisCancer, 'GroupCount', 'descend');

metSite_split=clinicalMetadata.metSiteExtended;
metSite_split_all=[];
for i=1:length(metSite_split)
    getsplit=split(metSite_split(i),",");
    getsplit=strtrim(getsplit);
    metSite_split_all=[metSite_split_all;getsplit];
end
metSite_split_all=table(metSite_split_all);
tblTopMetsThisCancer = groupcounts(metSite_split_all, 'metSite_split_all');
tblTopMetsThisCancer = sortrows(tblTopMetsThisCancer, 'GroupCount', 'descend');
%combine brain, brain/cns
tblTopMetsThisCancer.GroupCount(7)=3;
tblTopMetsThisCancer(15,:)=[];
%remove n/a
tblTopMetsThisCancer(1,:)=[];
%remove pleural effusion for simplicity since it's not exactly an organ
tblTopMetsThisCancer(5,:)=[];
tblTopMetsThisCancer = sortrows(tblTopMetsThisCancer, 'GroupCount', 'descend');
%% plot locations
figure(200)
bar(tblTopMetsThisCancer.GroupCount)
xlabel('Metastasis site')
ylabel('Cases')
h = gca;
h.XTick = 1:(height(tblTopMetsThisCancer));
h.XTickLabel = tblTopMetsThisCancer.metSite_split_all;
h.XTickLabelRotation = 90;
title('BRCA')

%% top mutations in this cancer
m = mutationMatrix{:, 2:end};
m(isnan(m)) = 0;
mutationMatrix{:, 2:end} = m;
disp('done.');
[nCases, orderTopMutations] = sort(sum(m>0), 'descend');

% bar plot with top most frequent mutations
ngenes = 50;

figure(201)
bar(nCases(1:ngenes))
xlabel('Gene mutated')
ylabel('Cases')
h = gca;
h.XTick = 1:ngenes;
h.XTickLabel = mutationMatrix.Properties.VariableNames(orderTopMutations(1:ngenes)+1);
h.XTickLabelRotation = 90;
title('BRCA')

%% find top sites
bone=contains(metSite_split,'BONE');
liver=contains(metSite_split,'LIVER');
lung=contains(metSite_split,'LUNG');
brain=contains(metSite_split,'BRAIN');
rln=contains(metSite_split,'REGIONAL LYMPH NODE');
metmatrix=[bone,liver,lung,brain];
metmatrix=[metmatrix,sum(metmatrix(:,1:4),2)];

%% add stacked mutation info to top met sites
pik3ca=mutationMatrix{:, 'PIK3CA'}==1 & mutationMatrix{:, 'TP53'}~=1;
tp53=mutationMatrix{:, 'TP53'}==1 & mutationMatrix{:, 'PIK3CA'}~=1;
pik3ca_tp53=mutationMatrix{:, 'TP53'}==1 & mutationMatrix{:, 'PIK3CA'}==1;
none_pik3ca_tp53=mutationMatrix{:, 'TP53'}~=1 & mutationMatrix{:, 'PIK3CA'}~=1;
metmatrix2=[bone,liver,rln,lung,brain];

pik3cakeep=metmatrix2(find(pik3ca),:);
pik3caprop=[sum(pik3cakeep,1)]; %bone, liver, rln, lung, brain
tp53keep=metmatrix2(find(tp53),:);
tp53prop=[sum(tp53keep,1)];
bothmutkeep=metmatrix2(find(pik3ca_tp53),:);
bothmutprop=[sum(bothmutkeep,1)];
nonemutkeep=metmatrix2(find(none_pik3ca_tp53),:);
nonemutprop=[sum(nonemutkeep,1)];

y = [pik3caprop; tp53prop; bothmutprop; nonemutprop];
figure,bar(y,'stacked')
legend({'bone','liver','regional lymph node','lung','brain/cns'},'Location','northeast')
xlabel('Mutation Status')
ylabel('Cases')
h = gca;
h.XTickLabel = {'PIK3CA','TP53','BOTH','NEITHER'};
h.XTickLabelRotation = 90;
title('BRCA')

%% add stacked hormone receptor info to top met sites
hormone=[];
hormone=strcmp(clinicalMetadata.PRDEverTripleNegative,'YES'); %tripleneg
hormone=[hormone,strcmp(clinicalMetadata.PRDEverHER2Positive,'YES')&strcmp(clinicalMetadata.PRDEverHormoneReceptorPositive,'NO');]; %her2+
hormone=[hormone,strcmp(clinicalMetadata.PRDEverHormoneReceptorPositive,'YES')&strcmp(clinicalMetadata.PRDEverHER2Positive,'NO');]; %er/pr+
hormone=[hormone,strcmp(clinicalMetadata.PRDEverHER2Positive,'YES')&strcmp(clinicalMetadata.PRDEverHormoneReceptorPositive,'YES');]; %er/pr/her2+


tripleneg=hormone(:,1);
her2=hormone(:,2);
erpr=hormone(:,3);
erprher2=hormone(:,4);

triplenegkeep=metmatrix2(find(tripleneg),:);
triplenegprop=[sum(triplenegkeep,1)]; %bone, liver, rln, lung, brain
her2keep=metmatrix2(find(her2),:);
her2prop=[sum(her2keep,1)];
erprkeep=metmatrix2(find(erpr),:);
erprprop=[sum(erprkeep,1)];
erprher2keep=metmatrix2(find(erprher2),:);
erprher2prop=[sum(erprher2keep,1)];

y = [triplenegprop; her2prop; erprprop; erprher2prop];
figure,bar(y,'stacked')
legend({'bone','liver','regional lymph node','lung','brain/cns'},'Location','northeast')
xlabel('Hormone Receptor Status')
ylabel('Cases')
h = gca;
h.XTickLabel = {'TripleNegative','Her2+','ER/PR+','ER/PR/Her2+'};
h.XTickLabelRotation = 90;
title('BRCA')

%% archetypes
[W2,H2] = nnmf(mExpression',6);
[M,arch]=max(W2,[],2);
archetype=zeros(113,6);
for i=1:113
    archetype(i,arch(i))=1;
end

arch1=archetype(:,1);
arch2=archetype(:,2);
arch3=archetype(:,3);
arch4=archetype(:,4);
arch5=archetype(:,5);
arch6=archetype(:,6);

arch1keep=metmatrix2(find(arch1),:);
arch1prop=[sum(arch1keep,1)]; %bone, liver, rln, lung, brain
arch2keep=metmatrix2(find(arch2),:);
arch2prop=[sum(arch2keep,1)];
arch3keep=metmatrix2(find(arch3),:);
arch3prop=[sum(arch3keep,1)];
arch4keep=metmatrix2(find(arch4),:);
arch4prop=[sum(arch4keep,1)];
arch5keep=metmatrix2(find(arch5),:);
arch5prop=[sum(arch5keep,1)];
arch6keep=metmatrix2(find(arch6),:);
arch6prop=[sum(arch6keep,1)];

y = [arch1prop; arch2prop; arch3prop; arch4prop; arch5prop;arch6prop];
figure,bar(y,'stacked')
legend({'bone','liver','regional lymph node','lung','brain/cns'},'Location','northeast')
xlabel('Transcriptomic Archetype')
ylabel('Cases')
h = gca;
h.XTickLabel = {'Arch1','Arch2','Arch3','Arch4','Arch5','Arch6'};
h.XTickLabelRotation = 90;
title('BRCA')

%plot archetypes on umap
expression_plus=[mExpression';H2]; %add perfect fit archetype samples to expression set
[W4,H4] = nnmf(expression_plus,6);
[M4,arch4]=max(W4,[],2);
archetype2=zeros(119,6);
for i=1:119
    archetype2(i,arch4(i))=1;
end
%z-score data
expression_plusZ = zscore(expression_plus,0,2);
[score, umap, clusterIdentifiers] = run_umap(expression_plusZ);

figure
h = gscatter(score(:, 1), score(:, 2), archetype2);
hold on
scatter(score(114:119,1),score(114:119,2),200,'kx');
hold off
legend({'arch1','arch2','arch3','arch4','arch5','arch6','archetype centers'})
xlabel('UMAP1')
ylabel('UMAP2')

%% make LDH and PDH tables
patients=metdata_rna.Properties.VariableNames;
patients=erase(patients,'MBC_');
patients=extractBefore(patients,'_Tumor');
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

%% overlay top ldh/pdh ratios on umap

[B,I] = sort(ratio);

%plot color range
figure,scatter(score(1:113, 1), score(1:113, 2), [], ratio, 'filled')
set(gca, 'CLim', [0, 8])
%circle lung
hold on
scatter(score(find(lung),1),score(find(lung),2),100,'ko'); %lung
hold off