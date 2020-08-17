% analyze growth rates of cell lines, data from 50mM LA normoxia titration
% (plain media condition)

%% input data
	datagrowthinput=[595	311	288	406	177	226	423	211	147
	1051	319	232	661	99	139	496	161	258
	1828	661	486	1299	157	255	991	351	518
	3141	1282	1003	2391	408	594	1817	726	979
	5730	2436	1906	5258	1440	1842	4229	2468	2549];

datagrowth=[];
for i=1:9
    datagrowth=[datagrowth;datagrowthinput(:,i)];
end

t=[(0:3),5]';
time=repmat(t,[9,1]);
lt=length(t);

wellID=[1:9*lt]';

replicate1=ones(lt,1);
replicate2=ones(lt,1)*2;
replicate3=ones(lt,1)*3;
replicates=[replicate1;replicate2;replicate3];
replicate=repmat(replicates,[3,1]);       

parental=repmat('P',[lt*3,1]);
brm2=repmat('B',[lt*3,1]);
lm2=repmat('L',[lt*3,1]);
celltype=[parental;brm2;lm2];

growthtable=table(wellID,time,celltype,replicate,datagrowth);

celltype=cellstr(celltype);
celltypecat=categorical(celltype);
growthtable.celltypecat=celltypecat;

cells = ordinal(celltypecat);
getlevels(cells);
cells= reorderlevels(cells,{'P','B', 'L'});
getlevels(cells);
growthtable.cells=cells;


growthtable.logdata = log(growthtable.datagrowth);
%% regression

mdl = fitlme(growthtable,'logdata~time*cells+(1|wellID)');

figure
plot(t,mean(datagrowthinput(:,1:3),2))
hold on
plot(t,mean(datagrowthinput(:,4:6),2))
plot(t,mean(datagrowthinput(:,7:9),2))
hold off
