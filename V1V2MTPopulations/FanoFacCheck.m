function [] = FanoFacCheck(correct_screen_disparity,subsample)

close all

% Define path to saved distribution data
splPath  = regexp(which('FanoFacCheck'),filesep,'split');
rootDir  = [filesep,fullfile(splPath{1:numel(splPath)-1}),filesep];
analyDir = [rootDir,'analysisFiles',filesep];

% Load in data
[MT] = loadDataMT(correct_screen_disparity);
[V1] = loadDataV1V2(1,correct_screen_disparity);
[V2] = loadDataV1V2(2,correct_screen_disparity);

% Subsample data if requested
if subsample
    load([analyDir,'eccInds'])

    MT = MT(eccInds.MT);
    V1 = V1(eccInds.V1);
    V2 = V2(eccInds.V2);
end

% Toggle save on/off
saveOn = 1;


%% Concatenate data 
% (mean & variance for each disparity presented to each cell - 
% (i.e. one cell's data appears more than once)

MTdat = [];
V1dat = [];
V2dat = [];

MTids = [];
V1ids = [];
V2ids = [];

for ii = 1:numel(MT)
    MTdat = [MTdat; MT{ii}.dat(:,2:4)];
    MTids = [MTids; ii*ones(size(MT{ii}.dat(:,2:4),1),1)];
end

for ii = 1:numel(V1)
    V1dat = [V1dat; V1{ii}.dat(:,2:4)];
    V1ids = [V1ids; ii*ones(size(V1{ii}.dat(:,2:4),1),1)];
end

for ii = 1:numel(V2)
    V2dat = [V2dat; V2{ii}.dat(:,2:4)];
    V2ids = [V2ids; ii*ones(size(V2{ii}.dat(:,2:4),1),1)];
end

% Cull data based on number of repeats
minReps = 5;
maxVar  = 1000;

V1dat = V1dat(V1dat(:,3) >= minReps,:);
V2dat = V2dat(V2dat(:,3) >= minReps,:);
MTdat = MTdat(MTdat(:,3) >= minReps,:);

V1ids = V1ids(V1dat(:,3) >= minReps,:);
V2ids = V2ids(V2dat(:,3) >= minReps,:);
MTids = MTids(MTdat(:,3) >= minReps,:);

V1dat = V1dat(V1dat(:,2) <= maxVar,:);
V2dat = V2dat(V2dat(:,2) <= maxVar,:);
MTdat = MTdat(MTdat(:,2) <= maxVar,:);

V1ids = V1ids(V1dat(:,2) <= maxVar,:);
V2ids = V2ids(V2dat(:,2) <= maxVar,:);
MTids = MTids(MTdat(:,2) <= maxVar,:);

% Cull data with zero responses
V1dat = V1dat(V1dat(:,2) > 0,:);
V2dat = V2dat(V2dat(:,2) > 0,:);
MTdat = MTdat(MTdat(:,2) > 0,:);

V1ids = V1ids(V1dat(:,2) > 0,:);
V2ids = V2ids(V2dat(:,2) > 0,:);
MTids = MTids(MTdat(:,2) > 0,:);

% Find number of cells from each area
MTcells = numel(unique(MTids));
V1cells = numel(unique(V1ids));
V2cells = numel(unique(V2ids));


%% Find best fitting power law

catDat = {V1dat,V2dat,MTdat};

pLaw   = @(a,x,b) a*x.^b;
p0     = [1 1];
opts   = optimoptions('fmincon','Display','off');
A      = [0,-1];
b      = -1;

for ii = 1:3
    data = catDat{ii};

    objFxn = @(p) sum( ( log(pLaw(p(1),data(:,1),p(2))) - log(data(:,2)) ).^2 );
    
    pFit(ii,:) = fmincon(objFxn,p0,A,b,[],[],[],[],[],opts);
end


%% Plot
supp = linspace(1e-1,1e3,100);

xtix = [0.01 0.1 1 10 100 1000];
for ii = 1:numel(xtix)
    xLab{ii} = num2str(xtix(ii));
end

f1 = figure;
f1.Position = [2000 100 1700 450];

% V1 distribution
subplot(1,3,1);
hold on;

scatter(V1dat(:,1),V1dat(:,2),70,[0.6 0.6 0.6]);
plot(supp,pLaw(pFit(1,1),supp,pFit(1,2)),'r','LineWidth',2);
plot([1e-1 1000],[1e-1 1000],'--k','LineWidth',2);
set(gca,'plotboxaspectratio',[1 1 1],'xtick',xtix,'xticklabel',xLab,'xlim',[1e-1 1000],'ylim',[1e-1 1000],'fontsize',15,'yscale','log','xscale','log');
xlabel('mean spike count');
ylabel('spike count variance');
title(['V1 (n=',num2str(sum(eccInds.V1,2)),')']);

% V2 distribution
subplot(1,3,2);
hold on;

scatter(V2dat(:,1),V2dat(:,2),70,[0.6 0.6 0.6]);
plot(supp,pLaw(pFit(2,1),supp,pFit(2,2)),'r','LineWidth',2);
plot([1e-1 1000],[1e-1 1000],'--k','LineWidth',2);
set(gca,'plotboxaspectratio',[1 1 1],'xtick',xtix,'xticklabel',xLab,'xlim',[1e-1 1000],'ylim',[1e-1 1000],'fontsize',15,'yscale','log','xscale','log');
xlabel('mean spike count');
ylabel('spike count variance');
title(['V2 (n=',num2str(sum(eccInds.V2,2)),')']);

% MT distribution
subplot(1,3,3);
hold on;

scatter(MTdat(:,1),MTdat(:,2),70,[0.6 0.6 0.6]);
plot(supp,pLaw(pFit(3,1),supp,pFit(3,2)),'r','LineWidth',2);
plot([1e-1 1000],[1e-1 1000],'--k','LineWidth',2);
set(gca,'plotboxaspectratio',[1 1 1],'xtick',xtix,'xticklabel',xLab,'xlim',[1e-1 1000],'ylim',[1e-1 1000],'fontsize',15,'yscale','log','xscale','log');
xlabel('mean spike count');
ylabel('spike count variance');
title(['MT (n=',num2str(sum(eccInds.MT,2)),')']);

% Plot best-fitting power law
f2 = figure;
f2.Position = [2000 400 1200 450];

subplot(1,2,1);
hold on;
X = categorical({'V1','V2','MT'});
X = reordercats(X,{'V1','V2','MT'});
Y = pFit(:,1);
bar(X,Y,'FaceColor',[0.5 0.5 0.5],'EdgeColor','none');
plot([0 4],[1 1],'--k','linewidth',2);
ylabel('Slope');
set(gca,'plotboxaspectratio',[1 1 1],'xlim',{'V1','MT'},'fontsize',15);

subplot(1,2,2);
hold on;
Y = pFit(:,2);
bar(X,Y,'FaceColor',[0.5 0.5 0.5],'EdgeColor','none');
plot([0 4],[1 1],'--k','linewidth',2);
ylabel('Power law');
set(gca,'plotboxaspectratio',[1 1 1],'xlim',{'V1','MT'},'fontsize',15);

%% Save figures

if saveOn
    saveas(f1,[rootDir,'plots/FanoFactor_scatter.svg']);
    saveas(f2,[rootDir,'plots/FanoFactor_fits.svg']);
end

