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
% remove from transcriptomics
mExpression(:, ~ismember(RNAseq_participants, patient)) = [];
RNAseq_participants(~ismember(RNAseq_participants, patient)) = [];
% remove from clinical metadata
clinicalMetadata(~ismember(clinicalMetadata.participant, patient), :) = [];
% remove from genomics
mutationMatrix(~ismember(mutationMatrix.participant, patient), :) = [];

%% reorder patients in clinical to match genomics
[~, patientReorder] = ismember(mutationMatrix.participant, clinicalMetadata.participant);
clinicalMetadata = clinicalMetadata(patientReorder, :);

%reorder patients in rna to match genomics
[~, patientReorder] = ismember(mutationMatrix.participant, RNAseq_participants);
RNAseq_participants=RNAseq_participants(patientReorder);
metdata_rna = metdata_rna(:,patientReorder);
mExpression=table2array(metdata_rna);
%% ldh/pdh ratio
ldha=strcmp(geneNames,'LDHA');
ldhb=strcmp(geneNames,'LDHB');
pdha1=strcmp(geneNames,'PDHA1');
pdha2=strcmp(geneNames,'PDHA2');
ratiocheck=(mExpression(ldha,:)+mExpression(ldhb,:))./(mExpression(pdha1,:)+mExpression(pdha2,:));
ratiocheck=ratiocheck';
%% Where do cancers metastasize?
clinicalMetadata.metSiteExtended = clinicalMetadata.MedREverMetastaticSites;
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
tblTopMetsThisCancer.GroupCount(7)=7;
tblTopMetsThisCancer(11,:)=[];
%remove n/a
tblTopMetsThisCancer(2,:)=[];
tblTopMetsThisCancer = sortrows(tblTopMetsThisCancer, 'GroupCount', 'descend');

%% plot locations
figure(200)
bar(tblTopMetsThisCancer.GroupCount)
xlabel('Metastasis site')
ylabel('Cases')
h = gca;
h.XTick = 1:(height(tblTopMetsThisCancer)-1);
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

%% add stacked mutation info to top met sites
bone=contains(metSite_split,'BONE');
liver=contains(metSite_split,'LIVER');
lung=contains(metSite_split,'LUNG');
brain=contains(metSite_split,'BRAIN');
metmatrix=[bone,liver,lung,brain];
metmatrix=[metmatrix,sum(metmatrix(:,1:4),2)];

mut_top=mutationMatrix{:, 'TP53'}==1 & mutationMatrix{:, 'PIK3CA'}~=1;
mut_top=[mut_top,mutationMatrix{:, 'PIK3CA'}==1 & mutationMatrix{:, 'TP53'}~=1];
mut_top=[mut_top,mutationMatrix{:, 'PIK3CA'}.*mutationMatrix{:, 'TP53'}==1]; 
rln=contains(metSite_split,'REGIONAL LYMPH NODE');

bonekeep=mut_top(find(bone),:);
boneprop=[sum(bonekeep,1),length(bonekeep)-sum(sum(bonekeep,1))]; %p53, pik3ca,both,none
brainkeep=mut_top(find(brain),:);
brainprop=[sum(brainkeep,1),length(brainkeep)-sum(sum(brainkeep,1))]; 
lungkeep=mut_top(find(lung),:);
lungprop=[sum(lungkeep,1),length(lungkeep)-sum(sum(lungkeep,1))]; 
liverkeep=mut_top(find(liver),:);
liverprop=[sum(liverkeep,1),length(liverkeep)-sum(sum(liverkeep,1))]; 
rlnkeep=mut_top(find(rln),:);
rlnprop=[sum(rlnkeep,1),length(rlnkeep)-sum(sum(rlnkeep,1))]; 

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

bonekeep=hormone(find(bone),:);
boneprop=[sum(bonekeep,1)]; %tripleneg,her2,er/pr,er/pr/her2
brainkeep=hormone(find(brain),:);
brainprop=[sum(brainkeep,1)]; 
lungkeep=hormone(find(lung),:);
lungprop=[sum(lungkeep,1)]; 
liverkeep=hormone(find(liver),:);
liverprop=[sum(liverkeep,1)]; 
rlnkeep=hormone(find(rln),:);
rlnprop=[sum(rlnkeep,1)]; 

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
%% nnmf
[W2,H2] = nnmf(mExpression',6);
[M,arch]=max(W2,[],2);
archetype=zeros(113,6);
for i=1:113
    archetype(i,arch(i))=1;
end
bonekeep=archetype(find(bone),:);
boneprop=[sum(bonekeep,1)]; %tripleneg,her2,er/pr,er/pr/her2
brainkeep=archetype(find(brain),:);
brainprop=[sum(brainkeep,1)]; 
lungkeep=archetype(find(lung),:);
lungprop=[sum(lungkeep,1)]; 
liverkeep=archetype(find(liver),:);
liverprop=[sum(liverkeep,1)]; 
rlnkeep=archetype(find(rln),:);
rlnprop=[sum(rlnkeep,1)]; 

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

%plot archetypes on PCA
expression_plus=[mExpression';H2]; %add perfect fit archetype samples to expression set
D = pdist(expression_plus);
Z=squareform(D);
[coeff,score,latent,tsquared,explained,mu] = pca(Z);
[W4,H4] = nnmf(expression_plus,6);
[M4,arch4]=max(W4,[],2);
archetype2=zeros(119,6);
for i=1:119
    archetype2(i,arch4(i))=1;
end

expression_plusZ = zscore(expression_plus,0,2);
D2 = pdist(expression_plusZ);
Z2=squareform(D2);
[coeff,score,latent,tsquared,explained,mu] = pca(Z2);
[score, umap, clusterIdentifiers] = run_umap(expression_plusZ);

figure
h = gscatter(score(:, 1), score(:, 2), archetype2);
hold on
scatter(score(114,1),score(114,2),200,'mx');
scatter(score(115,1),score(115,2),200,'cx');
scatter(score(116,1),score(116,2),200,'bx');
scatter(score(117,1),score(117,2),200,'gx');
scatter(score(118,1),score(118,2),200,'yx');
scatter(score(119,1),score(119,2),200,'rx');
hold off
legend({'arch1','arch2','arch3','arch4','arch5','arch6','archetype centers'})
xlabel('UMAP1')
ylabel('UMAP2')

