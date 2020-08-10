%% analyze hexokinase data

hexokinase_50min=xlsread('hexokinase50min.xlsx','A38:BA43');
standard=[];
j=1;
for i=3:8
    standard(j,2)=mean([hexokinase_50min(5,i),hexokinase_50min(5,i+7)]);
    j=j+1;
end
%subtract blank from standard
blank=standard(1,2);
standard(:,2)=standard(:,2)-standard(1,2);
standard(:,1)=[0;2.5;5;7.5;10;12.5];
figure,scatter(standard(:,1),standard(:,2))
% y = 0.061*x - 0.011

%subtract blank from positive control wells
pos=[hexokinase_50min(:,9),hexokinase_50min(:,16)];
pos=pos-blank;
pos=mean(pos,2);

%subtract blank from cell wells
celldata=hexokinase_50min(:,17:end);
celldata=celldata-blank;

%subtract cell background
j=1;
for i=25:30
    background(:,j)=mean([celldata(:,i),celldata(:,i+6)],2);
    j=j+1;
end

for i=1:24
    num=rem(i,6);
    if num==0
        num=6;
    end
    if any(background(:,num)<0)
    else
        celldata(:,i)=celldata(:,i)-background(:,num);
    end
end

%calculate activity
parentalhex=[];
brm2hex=[];
lm2hex=[];
parentalhalfhex=[];
brm2halfhex=[];
lm2halfhex=[];

for i=1:6:24
    parentalhex=[parentalhex;(celldata(end,i)-celldata(1,i))];
    brm2hex=[brm2hex;(celldata(end,i+1)-celldata(1,i+1))];
    lm2hex=[lm2hex;(celldata(end,i+2)-celldata(1,i+2))];
    parentalhalfhex=[parentalhalfhex;(celldata(end,i+3)-celldata(1,i+3))];
    brm2halfhex=[brm2halfhex;(celldata(end,i+4)-celldata(1,i+4))];
    lm2halfhex=[lm2halfhex;(celldata(end,i+5)-celldata(1,i+5))];
end

tdiff=hexokinase_50min(end,1)-hexokinase_50min(1,1);

parentalhex=(((parentalhex+.011)/.061)/tdiff)/.0375;
brm2hex=(((brm2hex+.011)/.061)/tdiff)/.0375;
lm2hex=(((lm2hex+.011)/.061)/tdiff)/.0375;
parentalhalfhex=(((parentalhalfhex+.011)/.061)/tdiff)/0.01875; %diluted cell volume
brm2halfhex=(((brm2halfhex+.011)/.061)/tdiff)/0.01875;
lm2halfhex=(((lm2halfhex+.011)/.061)/tdiff)/0.01875;
%above is divided by .0375 or .01875 to get activity in nmol product per million cells (per second)


%plot
x = 1:3;
data = [mean(parentalhalfhex);mean(brm2halfhex);mean(lm2halfhex)];
err = [std(parentalhalfhex)/2, std(brm2halfhex)/2, std(lm2halfhex)/2];

figure,bar(x,data)   
title('hexokinase')
ylabel('nmol/million cells/second')

hold on
er = errorbar(x,data,err,err);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
hold off

[h,ppb] = ttest(parentalhalfhex,brm2halfhex); %.197 .4 for half
[h,ppl] = ttest(parentalhalfhex,lm2halfhex); %.00003 .00007forhalf
[h,pbl] = ttest(lm2halfhex,brm2halfhex); %.00002 .00006for half
FC_BP=mean(brm2halfhex)/mean(parentalhalfhex); %1.05
FC_LP=mean(lm2halfhex)/mean(parentalhalfhex); %7.43
FC_LB=mean(lm2halfhex)/mean(brm2halfhex); %7.04


%regression
tabledata=[parentalhalfhex;brm2halfhex;lm2halfhex];
celltype=[repmat('P',[4,1]);repmat('B',[4,1]);repmat('L',[4,1])];
tablehex=table(celltype,tabledata);
celltype=cellstr(celltype);
celltypecat=categorical(celltype);
growth.celltypecat=celltypecat;
cells = ordinal(celltypecat);
getlevels(cells);
cells= reorderlevels(cells,{'P','B', 'L'});
getlevels(cells);
tablehex.cells=cells;
mdl=fitlme(tablehex,'tabledata~cells');
FCB=(mdl.Coefficients(1,2).Estimate+mdl.Coefficients(2,2).Estimate)/mdl.Coefficients(1,2).Estimate; %1.05
pvalB=mdl.Coefficients(2,6).pValue;
FCL=(mdl.Coefficients(1,2).Estimate+mdl.Coefficients(3,2).Estimate)/mdl.Coefficients(1,2).Estimate; %7.43
pvalL=mdl.Coefficients(3,6).pValue;

%% analyze pfk data

pfk_5min=xlsread('PFK5min.xlsx','A38:BA43');
standard=[];
j=1;
for i=3:8
    standard(j,2)=mean([pfk_5min(5,i),pfk_5min(5,i+7)]);
    j=j+1;
end
%subtract blank from standard
blank=standard(1,2);
standard(:,2)=standard(:,2)-standard(1,2);
standard(:,1)=[0;2.5;5;7.5;10;12.5];
figure,scatter(standard(:,1),standard(:,2))
% y = 0.037*x - 0.022 

%subtract blank from positive control wells
pos=[pfk_5min(:,9),pfk_5min(:,16)];
pos=pos-blank;
pos=mean(pos,2);

%subtract blank from cell wells
celldata=pfk_5min(:,17:end);
celldata=celldata-blank;

%subtract cell background
j=1;
for i=25:30
    background(:,j)=mean([celldata(:,i),celldata(:,i+6)],2);
    j=j+1;
end

for i=1:24
    num=rem(i,6);
    if num==0
        num=6;
    end
    if any(background(:,num)<0)
    else
        celldata(:,i)=celldata(:,i)-background(:,num);
    end
end

%calculate activity
parentalpfk=[];
brm2pfk=[];
lm2pfk=[];
parentalhalfpfk=[];
brm2halfpfk=[];
lm2halfpfk=[];


for i=1:6:24
    parentalpfk=[parentalpfk;(celldata(1,i))];
    brm2pfk=[brm2pfk;(celldata(1,i+1))];
    lm2pfk=[lm2pfk;(celldata(1,i+2))];
    parentalhalfpfk=[parentalhalfpfk;(celldata(1,i+3))];
    brm2halfpfk=[brm2halfpfk;(celldata(1,i+4))];
    lm2halfpfk=[lm2halfpfk;(celldata(1,i+5))];
end
tdiff=5*60;

parentalpfk=(((parentalpfk+.022)/.037)/tdiff)/.0375;
brm2pfk=(((brm2pfk+.022)/.037)/tdiff)/.0375;
lm2pfk=(((lm2pfk+.022)/.037)/tdiff)/.0375;
parentalhalfpfk=(((parentalhalfpfk+.022)/.037)/tdiff)/0.01875;
brm2halfpfk=(((brm2halfpfk+.022)/.037)/tdiff)/0.01875;
lm2halfpfk=(((lm2halfpfk+.022)/.037)/tdiff)/0.01875;
%above is divided by .0375 or .01875 to get activity in nmol product per million cells (per second)


%plot
x = 1:3;
data = [mean(parentalhalfpfk);mean(brm2halfpfk);mean(lm2halfpfk)];
err = [std(parentalhalfpfk)/2, std(brm2halfpfk)/2, std(lm2halfpfk)/2];

figure,bar(x,data)         
title('PFK')
ylabel('nmol/million cells/second')

hold on
er = errorbar(x,data,err,err);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
hold off

[h,ppb] = ttest(parentalhalfpfk,brm2halfpfk); %.003
[h,ppl] = ttest(parentalhalfpfk,lm2halfpfk); %.003
[h,pbl] = ttest(lm2halfpfk,brm2halfpfk); %.004
FC_BP=mean(brm2halfpfk)/mean(parentalhalfpfk); %1.8
FC_LP=mean(lm2halfpfk)/mean(parentalhalfpfk); %3.3
FC_LB=mean(lm2halfpfk)/mean(brm2halfpfk); %1.7

% regression
tabledata=[parentalhalfpfk;brm2halfpfk;lm2halfpfk];
celltype=[repmat('P',[4,1]);repmat('B',[4,1]);repmat('L',[4,1])];
tablepfk=table(celltype,tabledata);
celltype=cellstr(celltype);
celltypecat=categorical(celltype);
growth.celltypecat=celltypecat;
cells = ordinal(celltypecat);
getlevels(cells);
cells= reorderlevels(cells,{'P','B', 'L'});
getlevels(cells);
tablepfk.cells=cells;
mdl=fitlme(tablepfk,'tabledata~cells');
FCB=(mdl.Coefficients(1,2).Estimate+mdl.Coefficients(2,2).Estimate)/mdl.Coefficients(1,2).Estimate; %1.8
pvalB=mdl.Coefficients(2,6).pValue;
FCL=(mdl.Coefficients(1,2).Estimate+mdl.Coefficients(3,2).Estimate)/mdl.Coefficients(1,2).Estimate; %3.3
pvalL=mdl.Coefficients(3,6).pValue;

%% analyze g3pdh data

g3pdh_0min=xlsread('G3PDHREPEAT0MIN.xlsx','A38:BA43');
standard=[];
j=1;
for i=3:8
    standard(j,2)=mean([g3pdh_0min(5,i),g3pdh_0min(5,i+7)]);
    j=j+1;
end
%subtract blank from standard
blank=standard(1,2);
standard(:,2)=standard(:,2)-standard(1,2);
standard(:,1)=[0;2.5;5;7.5;10;12.5];
figure,scatter(standard(:,1),standard(:,2))
% y = 0.052*x - 0.0038 

%subtract blank from positive control wells
pos=[g3pdh_0min(:,9),g3pdh_0min(:,16)];
pos=pos-blank;
pos=mean(pos,2);

%subtract blank from cell wells
celldata=g3pdh_0min(:,17:end);
celldata=celldata-blank;

%subtract cell background
j=1;
for i=25:30
    background(:,j)=mean([celldata(:,i),celldata(:,i+6)],2);
    j=j+1;
end

for i=1:24
    num=rem(i,6);
    if num==0
        num=6;
    end
    if any(background(:,num)<0)
    else
        celldata(:,i)=celldata(:,i)-background(:,num);
    end
end

%calculate activity
parentalg3p=[];
brm2g3p=[];
lm2g3p=[];
parentalhalfg3p=[];
brm2halfg3p=[];
lm2halfg3p=[];

for i=1:6:24
    parentalg3p=[parentalg3p;(celldata(end,i)-celldata(1,i))];
    brm2g3p=[brm2g3p;(celldata(end,i+1)-celldata(1,i+1))];
    lm2g3p=[lm2g3p;(celldata(end,i+2)-celldata(1,i+2))];
    parentalhalfg3p=[parentalhalfg3p;(celldata(end,i+3)-celldata(1,i+3))];
    brm2halfg3p=[brm2halfg3p;(celldata(end,i+4)-celldata(1,i+4))];
    lm2halfg3p=[lm2halfg3p;(celldata(end,i+5)-celldata(1,i+5))];
end
tdiff=g3pdh_0min(end,1)-g3pdh_0min(1,1);

parentalg3p=(((parentalg3p+.0038)/.052)/tdiff)/.0375;
brm2g3p=(((brm2g3p+.0038)/.052)/tdiff)/.0375;
lm2g3p=(((lm2g3p+.0038)/.052)/tdiff)/.0375;
parentalhalfg3p=(((parentalhalfg3p+.0038)/.052)/tdiff)/0.01875;
brm2halfg3p=(((brm2halfg3p+.0038)/.052)/tdiff)/0.01875;
lm2halfg3p=(((lm2halfg3p+.0038)/.052)/tdiff)/0.01875;
%above is divided by .0375 or .01875 to get activity in nmol product per million cells (per second)


%plot
x = 1:3;
data = [mean(parentalhalfg3p);mean(brm2halfg3p);mean(lm2halfg3p)];
err = [std(parentalhalfg3p)/2, std(brm2halfg3p)/2, std(lm2halfg3p)/2];

figure,bar(x,data)  
title('G3PDH')
ylabel('nmol/million cells/second')

hold on
er = errorbar(x,data,err,err);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
hold off

[h,ppb] = ttest(parentalhalfg3p,brm2halfg3p); %.002
[h,ppl] = ttest(parentalhalfg3p,lm2halfg3p); %.00001
[h,pbl] = ttest(lm2halfg3p,brm2halfg3p); %.00003

FC_BP=mean(brm2halfg3p)/mean(parentalhalfg3p); %1.4
FC_LP=mean(lm2halfg3p)/mean(parentalhalfg3p); %2.3
FC_LB=mean(lm2halfg3p)/mean(brm2halfg3p); %1.6

% regression
tabledata=[parentalhalfg3p;brm2halfg3p;lm2halfg3p];
celltype=[repmat('P',[4,1]);repmat('B',[4,1]);repmat('L',[4,1])];
tableg3p=table(celltype,tabledata);
celltype=cellstr(celltype);
celltypecat=categorical(celltype);
growth.celltypecat=celltypecat;
cells = ordinal(celltypecat);
getlevels(cells);
cells= reorderlevels(cells,{'P','B', 'L'});
getlevels(cells);
tableg3p.cells=cells;
mdl=fitlme(tableg3p,'tabledata~cells');
FCB=(mdl.Coefficients(1,2).Estimate+mdl.Coefficients(2,2).Estimate)/mdl.Coefficients(1,2).Estimate; %1.3
pvalB=mdl.Coefficients(2,6).pValue;
FCL=(mdl.Coefficients(1,2).Estimate+mdl.Coefficients(3,2).Estimate)/mdl.Coefficients(1,2).Estimate; %2.3
pvalL=mdl.Coefficients(3,6).pValue;


%% analyze pk data

pk_fl0min=xlsread('PK0min.xlsx','A38:BA43');
standard=[];
j=1;
for i=3:8
    standard(j,2)=mean([pk_fl0min(5,i),pk_fl0min(5,i+7)]);
    j=j+1;
end
%subtract blank from standard
blank=standard(1,2);
standard(:,2)=standard(:,2)-standard(1,2);
standard(:,1)=[0;2.5;5;7.5;10;12.5];
figure,scatter(standard(:,1),standard(:,2))
% y = 0.1*x - 0.039 

%subtract blank from positive control wells
pos=[pk_fl0min(:,9),pk_fl0min(:,16)];
pos=pos-blank;
pos=mean(pos,2);

%subtract blank from cell wells
celldata=pk_fl0min(:,17:end);
celldata=celldata-blank;

%subtract cell background
j=1;
for i=25:30
    background(:,j)=mean([celldata(:,i),celldata(:,i+6)],2);
    j=j+1;
end

for i=1:24
    num=rem(i,6);
    if num==0
        num=6;
    end
    if any(background(:,num)<0)
    else
        celldata(:,i)=celldata(:,i)-background(:,num);
    end
end

%calculate activity
parentalpkfl=[];
brm2pkfl=[];
lm2pkfl=[];
parentalhalfpkfl=[];
brm2halfpkfl=[];
lm2halfpkfl=[];

for i=1:6:24
    parentalpkfl=[parentalpkfl;(celldata(end,i)-celldata(1,i))];
    brm2pkfl=[brm2pkfl;(celldata(end,i+1)-celldata(1,i+1))];
    lm2pkfl=[lm2pkfl;(celldata(end,i+2)-celldata(1,i+2))];
    parentalhalfpkfl=[parentalhalfpkfl;(celldata(end,i+3)-celldata(1,i+3))];
    brm2halfpkfl=[brm2halfpkfl;(celldata(end,i+4)-celldata(1,i+4))];
    lm2halfpkfl=[lm2halfpkfl;(celldata(end,i+5)-celldata(1,i+5))];
end
tdiff=pk_fl0min(end,1)-pk_fl0min(1,1);


parentalpkfl=(((parentalpkfl+.039)/.1)/tdiff)/.0375;
brm2pkfl=(((brm2pkfl+.039)/.1)/tdiff)/.0375;
lm2pkfl=(((lm2pkfl+.039)/.1)/tdiff)/.0375;
parentalhalfpkfl=(((parentalhalfpkfl+.039)/.1)/tdiff)/0.01875;
brm2halfpkfl=(((brm2halfpkfl+.039)/.1)/tdiff)/0.01875;
lm2halfpkfl=(((lm2halfpkfl+.039)/.1)/tdiff)/0.01875;
%above is divided by .0375 or .01875 to get activity in nmol product per million cells (per second)


%plot
x = 1:3;
data = [mean(parentalhalfpkfl);mean(brm2halfpkfl);mean(lm2halfpkfl)];
err = [std(parentalhalfpkfl)/2, std(brm2halfpkfl)/2, std(lm2halfpkfl)/2];

figure,bar(x,data)          
title('PK')
ylabel('nmol/million cells/second')

hold on
er = errorbar(x,data,err,err);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
hold off

[h,ppb] = ttest(parentalhalfpkfl,brm2halfpkfl); %.007 .8for0min
[h,ppl] = ttest(parentalhalfpkfl,lm2halfpkfl); %.03 .04for0min
[h,pbl] = ttest(lm2halfpkfl,brm2halfpkfl); %.01 .03for0min

FC_BP=mean(brm2halfpkfl)/mean(parentalhalfpkfl); %.98
FC_LP=mean(lm2halfpkfl)/mean(parentalhalfpkfl); %1.19
FC_LB=mean(lm2halfpkfl)/mean(brm2halfpkfl); %1.22

% regression
tabledata=[parentalhalfpkfl;brm2halfpkfl;lm2halfpkfl];
celltype=[repmat('P',[4,1]);repmat('B',[4,1]);repmat('L',[4,1])];
tablepkfl=table(celltype,tabledata);
celltype=cellstr(celltype);
celltypecat=categorical(celltype);
growth.celltypecat=celltypecat;
cells = ordinal(celltypecat);
getlevels(cells);
cells= reorderlevels(cells,{'P','B', 'L'});
getlevels(cells);
tablepkfl.cells=cells;
mdl=fitlme(tablepkfl,'tabledata~cells');
FCB=(mdl.Coefficients(1,2).Estimate+mdl.Coefficients(2,2).Estimate)/mdl.Coefficients(1,2).Estimate; %.98
pvalB=mdl.Coefficients(2,6).pValue; %.8
FCL=(mdl.Coefficients(1,2).Estimate+mdl.Coefficients(3,2).Estimate)/mdl.Coefficients(1,2).Estimate; %1.2
pvalL=mdl.Coefficients(3,6).pValue; %.008
