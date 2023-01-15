function [dispDat,realDat] = collectBootstrapRuns(numBoots,numHistBins)

%% Collect BORIS stats bootstrapping runs generated on HPC

% Define path to saved distribution data
splPath = regexp(which('collectBootstrapRuns'),filesep,'split');
rootDir = [filesep,fullfile(splPath{1:numel(splPath)-1}),filesep];
imStatsDir = [rootDir,'savedImageStats_BORISdataset',filesep,'borisStats',filesep];
saveDir    = [rootDir,'savedImageStats_BORISdataset',filesep];

% Define two behaviorally-ID'd image sets
imgSet    = {'sando','walking'};

% Define VF subregions (one for filenames, one for struct fieldnames)
% imRegion  = {'Upper VF_','Lower VF_','Central_','Peripheral_',''};
% imRegionF = {'UpperVF','LowerVF','Central','Peripheral','all'};
imRegion  = {''};
imRegionF = {'all'};

% Define suffices we used to sort the saved image stats from the HPC
% suff      = {'_m1','_m2','_m2','_m1',''};
suff      = {''};

% Confidence interval labels for each VF subregion
% imRegionLCI = {'UpperVFLCI','LowerVFLCI','CentralLCI','PeripheralLCI','allLCI'};
% imRegionUCI = {'UpperVFUCI','LowerVFUCI','CentralUCI','PeripheralUCI','allUCI'};
imRegionLCI = {'allLCI'};
imRegionUCI = {'allUCI'};

% Define how we defined the KSDs
KSDid     = {'Circ','V1','V2','MT','V1V2Rect'};


%% Loop over all saved .mat files parallel computing generated

% Define binning for how we calculate CIs
numCIBins = 10;

% Initialize loop variables
missingRuns = [];
missCnt = 1;
dispDat = struct;
realDat = struct;

% Loop over image sets (sando/walking)
for jj = 1:numel(imgSet)

    % Loop over KSDs
    for kk = 1:numel(KSDid)

        % Loop over VF subregions
        for ll = 1:numel(imRegion)

            % Load in the real stats
            realDatStruc = load([imStatsDir,'dispHist',KSDid{kk},'_',imRegion{ll},imgSet{jj},'.mat']);
            fieldName    = fieldnames(realDatStruc);
            thisDat      = getfield(realDatStruc,fieldName{1});

            realDat = setfield(realDat,imgSet{jj},KSDid{kk},imRegionF{ll},thisDat);

            % Loop over bootstraps
            theseRuns = nan(numBoots,numHistBins);

            for ii = 1:numBoots
                try
                    % Collect the data from this run's structure field
                    datStruc  = load([imStatsDir,'dispHist',KSDid{kk},'_',imRegion{ll},imgSet{jj},num2str(ii),'.mat']);
                    fieldName = fieldnames(datStruc);
                    thisRun   = getfield(datStruc,fieldName{1});

                    theseRuns(ii,:) = thisRun;

                catch
                    % ID any bootstrapping runs we might have lost due to a
                    % cluster issue so we can rerun them
                    missingRuns(missCnt).filename = ['dispHist',KSDid{kk},'_',imRegion{ll},imgSet{jj},num2str(ii)];
                    missCnt = missCnt + 1;

                end

            end

            % Pack this data into a new structure field
            dispDat = setfield(dispDat,imgSet{jj},KSDid{kk},imRegionF{ll},theseRuns);

            % Get 95% CIs for each of the histogram bins
            alpha = 0.05;

            for ii = 1:numHistBins

                [thisKSD,theseVals] = ksdensity(theseRuns(:,ii));

                thisCumDens = cumsum(thisKSD/sum(thisKSD));

                [~,LCIind] = min(abs(thisCumDens - alpha/2));
                [~,UCIind] = min(abs(thisCumDens - (1 - alpha/2)));

                theseLCI(ii) = theseVals(LCIind);
                theseUCI(ii) = theseVals(UCIind);

            end

            dispDat = setfield(dispDat,imgSet{jj},KSDid{kk},imRegionLCI{ll},theseLCI);
            dispDat = setfield(dispDat,imgSet{jj},KSDid{kk},imRegionUCI{ll},theseUCI);

        end

    end

end

% Report rumber of missing runs
if ~isempty(missingRuns)
    disp(['Total number of missing runs: ',num2str(missCnt)]);
    disp('- The following runs are missing: ');

    for ii = 1:missCnt
        disp(missingRuns(ii).filename);
    end

    save([saveDir,'missingRuns'],'missingRuns');
end

% Save the collected stats data
save([saveDir,'dispHistStruct'],'dispDat','realDat');

end

