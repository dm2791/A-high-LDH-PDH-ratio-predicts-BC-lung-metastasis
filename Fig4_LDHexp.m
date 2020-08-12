% analyze RNA expression data 231 cells
%need to download chip file; see Read Me for instructions
%% read files
chipfile=readtable('HG-U133A.na36annot.txt');
annotation=chipfile(26:end,1:19);

%annotate lung IDs to gene IDs
lunggsea=readtable('RNA_lungGSEA2.xlsx');
getIDslung=intersect(lunggsea.x_1_2,annotation.x__ForInformationAboutTheAnnotationFileContent);
rowsl=ismember(annotation.x__ForInformationAboutTheAnnotationFileContent,getIDslung);
geneIDlung=annotation.Var15(rowsl);

braingsea=readtable('RNA_brainGSEA2.xlsx');
getIDsbrain=intersect(braingsea.x_1_2,annotation.x__ForInformationAboutTheAnnotationFileContent);
rowsb=ismember(annotation.x__ForInformationAboutTheAnnotationFileContent,getIDsbrain);
geneIDbrain=annotation.Var15(rowsb);

%% 
glyclist=[{'LDHA'},{'LDHB'}];
glycbrain=[];
glyclung=[];
for i=1:length(glyclist)
    glycbrain=[glycbrain;find(strcmp(geneIDbrain,glyclist(i)))];
    glyclung=[glyclung;find(strcmp(geneIDlung,glyclist(i)))];
end

glycbrain=glycbrain+2;
glyclung=glyclung+2;

glycbraindata=table2array(braingsea(glycbrain,4:end));
glycbraindata=cellfun(@str2num,glycbraindata);
glyclungdata=table2array(lunggsea(glyclung,4:end));
glyclungdata=cellfun(@str2num,glyclungdata);

%FC
brainfc=mean(glycbraindata(1:2,4:5),2)./mean(glycbraindata(1:2,1:3),2);
lungfc=mean(glyclungdata(1:2,3:6),2)./mean(glyclungdata(1:2,1:2),2);
[h,pb]=ttest2(glycbraindata(1:2,4:5)',glycbraindata(1:2,1:3)');
[h,pl]=ttest2(glyclungdata(1:2,3:6)',glyclungdata(1:2,1:2)');

labels=geneIDlung(glyclung(1:2)-2);

figure
bar([mean(glycbraindata(1:2,1:3),2),mean(glycbraindata(1:2,4:5),2)])
figure
bar([mean(glyclungdata(1:2,1:2),2),mean(glyclungdata(1:2,3:6),2)])
