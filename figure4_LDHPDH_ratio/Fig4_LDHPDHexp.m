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

