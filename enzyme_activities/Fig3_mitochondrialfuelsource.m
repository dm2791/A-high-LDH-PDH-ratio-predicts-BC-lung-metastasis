% Analyze Seahorse Mito Flex assay

%% extract data
dependency=readmatrix('dependencymat.csv');
capacity=readmatrix('capacitymat.csv');
%% using only last measurement in interval
baselined=dependency(9:12,:);
baselinec=capacity(9:12,:);

drug1d=dependency(33:36,:);
drug1c=capacity(33:36,:);

drug2d=dependency(57:60,:);
drug2c=capacity(57:60,:);

dependencypercent=((baselined-drug1d)./(baselined-drug2d))*100;
capacitypercent=(1-((baselinec-drug1c)./(baselinec-drug2c)))*100;

depmean=mean(dependencypercent,1);
depstd=std(dependencypercent,1);
capmean=mean(capacitypercent,1);
capstd=std(capacitypercent,1);
flexibilitymean=capmean-depmean;

%% moles of usage 
dependency_total=baselined-drug1d;
glut=dependency_total(:,1:3); 
gluc=dependency_total(:,7:9);

meanglut=mean(glut);
meangluc=mean(gluc);
FCglut=[meanglut(2)/meanglut(1),meanglut(3)/meanglut(1),meanglut(3)/meanglut(2)];
FCgluc=[meangluc(2)/meangluc(1),meangluc(3)/meangluc(1),meangluc(3)/meangluc(2)];
[h,pbq]=ttest2(glut(:,1),glut(:,2));
[h,plq]=ttest2(glut(:,1),glut(:,3));
[h,blq]=ttest2(glut(:,2),glut(:,3));
[h,pbg]=ttest2(gluc(:,1),gluc(:,2));
[h,plg]=ttest2(gluc(:,1),gluc(:,3));
[h,blg]=ttest2(gluc(:,2),gluc(:,3));

figure, bar(meangluc)
figure, bar(meanglut)
