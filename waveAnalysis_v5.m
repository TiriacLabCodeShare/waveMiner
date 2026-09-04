%% Wave analysis pipeline v5
% By Kaylee Odum

% v2 - 250811 - incorporates local synchrony analysis (adapted by Alex Tiriac)
% v3 - 260617 - incorporates STTC analysis and optimizes output variable creation
% v4 - 260630 - optimizes code to be significantly faster
% v5 - 260708 - incorporates IBI, NBI, and II

clearvars, close all

%% User inputs

% Data table name:
dataTableName = "metaDataTable2.xlsx";

% To account for variability in var names:
baseline = "C"; % name of baseline recording in "condition" column
experimental = "post"; % name of +manipulation recording in "condition" column

% For processing raw data:
sampleRate = 20000;  % sample rate in Hz
activityThreshold = 0.5; % activity threshold for processing raw data
framesPerSec = 1; % 1 second intervals (used for all analyses except STTC)

% For wave and STTC analysis:
analyzeFiltData = 1; % toggle to analyze wave data with median filter applied (standard/default); set to 0 if you want to analyze raw wave data
pixelWidth = 87.5;
pixelHeight = 87.5;
pixelArea = pixelWidth*pixelHeight; % area of each pixel
kernelSize = 1; % +/- this value around current pixel for vector flow analysis
dt = 0.5; % in seconds (for STTC)

tWindow = 1; % +/- 1 frame

% For saving the results variable:
saveResults = 1; % 1 if you want to save the results 
tableFormat = 1; % 1 if you want the results var to be a table, 0 if you want a structure (structure = smaller file size)
todaysDate = "260824";



% Directories for finding and saving files:

% windows OS 
% dataTablePath = "C:\Users\yeagerkm\OneDrive - Vanderbilt\Tiriac Lab\Toolbox Project\5_increasedNoise\2_code"; 
% rawDataPath = "C:\Users\yeagerkm\OneDrive - Vanderbilt\Tiriac Lab\Toolbox Project\5_increasedNoise\1_data";
% waveVarPath = "C:\Users\yeagerkm\OneDrive - Vanderbilt\Tiriac Lab\Toolbox Project\5_increasedNoise\3_waveDetectionResults\1_waveData_mat";
% saveResultsLocation = "C:\Users\yeagerkm\OneDrive - Vanderbilt\Tiriac Lab\Toolbox Project\5_increasedNoise\4_waveAnalysisResults"; 

% mac OS 
dataTablePath = "/Users/kayleeodum/Library/CloudStorage/OneDrive-Vanderbilt/Tiriac Lab/Toolbox Project/3_wavesAcrossDevelopment/2_code"; 
rawDataPath = "/Users/kayleeodum/Library/CloudStorage/OneDrive-Vanderbilt/Tiriac Lab/Toolbox Project/3_wavesAcrossDevelopment/1_data";
waveVarPath = "/Users/kayleeodum/Library/CloudStorage/OneDrive-Vanderbilt/Tiriac Lab/Toolbox Project/3_wavesAcrossDevelopment/3_waveDetectionResults/1_waveData_mat";
saveResultsLocation = "/Users/kayleeodum/Library/CloudStorage/OneDrive-Vanderbilt/Tiriac Lab/Toolbox Project/3_wavesAcrossDevelopment/4_waveAnalysisResults";

%% Load the file table

% load meta data table
cd(dataTablePath)    
fileTable = readtable(dataTableName);

% extract file table columns into their own array variables
expNums = fileTable.expNum; % experiment number
fileNames = string(fileTable.h5FileName);
waveDataVarNames = string(fileTable.waveVarFileName);
ages = fileTable.age;
sexes = string(fileTable.sex);
eyes = string(fileTable.eye); % L or R eye
genotypes = string(fileTable.genotype); % + or - genotype depending on experiment
conditions = string(fileTable.condition); % baseline (pre) or +manipulation (post) group
timeVars = fileTable.recordingTime; % recording length in seconds
midlines = fileTable.midline; % horizontal midline of retina on MEA chip

%% Run the pipeline
tic 

numFiles = size(fileTable,1);

% initialize per file results vars
fileNum = nan(numFiles,1);
numWaves = nan(numFiles,1);
retinalArea = nan(numFiles,1);
IWIperROI = cell(numFiles,1);
avgIWI = nan(numFiles,1); % min/wave
frequency = nan(numFiles,1); % waves/min
avgDuration = nan(numFiles,1); % seconds
avgArea = nan(numFiles,1); % mm²
avgLocalSynch = nan(numFiles,1);
avgSpeed = nan(numFiles,1); % µm/s
avgTopSpeed = nan(numFiles,1);
avgDistance = nan(numFiles,1); % µm
flowMatrixX = cell(numFiles,1);
flowMatrixY = cell(numFiles,1);
IBI = nan(numFiles,1);
NBI = nan(numFiles,1);
AIFgradient = cell(numFiles,1);
STTC = cell(numFiles,1);
channelDistances = cell(numFiles,1);

% initialize per wave results vars
wFileNum = [];
wExpNum = [];
wAge = [];
wSex = {};
wEye = {};
wGenotype = {};
wCondition = {};
wMidline = [];
waveMatrix = {};
wStartFrame = [];
wEndFrame = [];
wDuration = []; % seconds
wAreaMatrix = {};
wVoxels = [];
wPixels = [];
wArea = []; % mm²
wRetAreaRatio = [];
wInitCoords = {};
wInitRegion = {};
wxVectors = {};
wyVectors = {};
wxVectorSum = [];
wyVectorSum = [];
wDirTheta = [];
wDirRho = [];
wSpeed = []; % µm/sec
wTopSpeed = [];
wDistTraveled = []; % µm
wLocalSynch = [];
wAIF = [];

% initialize index var for all wave vars
widx = 1;

%%%%% Analyze each file %%%%%
for currentFile = 1:numFiles
    currentFile

% extract meta data variables for current file
    expNum = expNums(currentFile);
    fileName = fileNames(currentFile);
    age = ages(currentFile);
    sex = sexes(currentFile);
    eye = eyes(currentFile);
    genotype = genotypes(currentFile);
    condition = conditions(currentFile);
    timeVar = timeVars(currentFile);
    midline = midlines(currentFile);

% load and process raw data
    % load the data
    cd(rawDataPath)
    spikeData = h5read(fileName,'/data_store/data0000/spikes');
    mappingData = h5read(fileName,'/data_store/data0000/settings/mapping');

    % Replace frameno with spikeTime (it's more intuitive)
    [spikeData.('spikeTime')] = spikeData.('frameno');
    spikeData = rmfield(spikeData,'frameno');

    % Start spikeTime at zero and convert into seconds
    spikeData.spikeTime = double(spikeData.spikeTime - min(spikeData.spikeTime))/sampleRate;
    
    % calc whole spikeRate for plot
    wholeRetSpikeRate = histcounts(spikeData.spikeTime,0:(1/framesPerSec):timeVar);

    % align mapping data to actual data
    channelIDs = unique(spikeData.channel);
    [channelIDs,~,ib] = intersect(channelIDs,mappingData.channel);
    mappingData.channel = mappingData.channel(ib);
    mappingData.electrode = mappingData.electrode(ib);
    mappingData.x = mappingData.x(ib);
    mappingData.y = mappingData.y(ib);
    
    % get rid of non mapping channels from data
    numUniqueChannelsInData = unique(spikeData.channel);
    for n = 1:length(numUniqueChannelsInData)
        if ismember(numUniqueChannelsInData(n), channelIDs) == 0
            chInd = find(spikeData.channel == numUniqueChannelsInData(n));
            spikeData.amplitude(chInd) = [];
            spikeData.channel(chInd) = [];
            spikeData.spikeTime(chInd) = [];
        end
    end
    
    spikeTimes = spikeData.spikeTime;
    channelIDs = unique(spikeData.channel); 
    numChannels = length(channelIDs);

    % calculate spikeRate per channel
    spikeRatePerChannel = zeros(numChannels,(timeVar*framesPerSec));
    for n = 1:numChannels
        currentChannel = channelIDs(n);
        spikeTimesCurrCh = spikeData.spikeTime(spikeData.channel==currentChannel);
        spikeRatePerChannel(n,:) = histcounts(spikeTimesCurrCh, 0:(1/framesPerSec):timeVar);
    end 

    % turn rasterplots into xyt matrix
    burstMatrix = spikeRatePerChannel > activityThreshold;
    
    uniqueXs = unique(mappingData.x);
    uniqueYs = unique(mappingData.y);
    numXs = length(uniqueXs);
    numYs = length(uniqueYs);
    
    spikeRateMovie = zeros(numYs,numXs,timeVar);
    burstMovie = zeros(numYs,numXs,timeVar);
    
    for n = 1:numChannels
        currentX = mappingData.x(n);
        currentY = mappingData.y(n);

        xInd = find(currentX == uniqueXs);
        yInd = find(currentY == uniqueYs);

        burstMovie(yInd,xInd,:) = burstMatrix(n,:);
        spikeRateMovie(yInd,xInd,:) = spikeRatePerChannel(n,:);
    end

    if strcmp(eye,"R") == 1
        burstMovie = fliplr(burstMovie);
        spikeRateMovie = fliplr(spikeRateMovie);
    end

% load wave data variables 
    cd(waveVarPath)
    load(waveDataVarNames(currentFile));
    if analyzeFiltData==1
        waveData = waveData_mfilt; % will use filtered wave data
        clear waveData_mfilt
    else
        clear waveData_mfilt % will use raw wave data (which is already named waveData when loaded)
    end

    % make sure empty pixels are nan, not 0
    waveData(waveData==0) = nan;

% align wave data matrix (temporal will always be on the left, assuming ventral is down)
    if strcmp(fileTable.eye(currentFile),'R')
        waveData = fliplr(waveData);
    end
    
% calculate total wave coverage of retina
    if strcmp(condition,baseline)
        coverageTempMat = sum(waveData,3,"omitnan"); % make it 2D
        coverageTempMat(coverageTempMat>0) = 1; % make each active pixel = 1
        totRetPixels = sum(coverageTempMat(coverageTempMat==1));
        totRetArea = (totRetPixels*pixelArea)/1000000; % convert to microns², then to mm²
    else
        totRetArea = retinalArea(expNums==expNum & strcmp(conditions,baseline));
    end

% create ID and size variables
    waveIDs = unique(waveData(~isnan(waveData)));
    nWaves = length(waveIDs);
    [ny,nx,nt] = size(waveData);
    minutes = timeVar/60;

% store wave summary meta data for current file
    fileNum(currentFile) = currentFile;
    numWaves(currentFile) = nWaves;
    retinalArea(currentFile) = totRetArea;

    %%%%% Analyze each wave %%%%%
    for currentWave = 1:nWaves
        % currentWave
        waveID = waveIDs(currentWave);

        % store meta data about the recording that the wave belongs to
        wFileNum(widx,1) = currentFile;
        wExpNum(widx,1) =  expNum;
        wAge(widx,1) = age;
        wSex{widx,1} = sex;
        wEye{widx,1} = eye;
        wGenotype{widx,1} = genotype;
        wCondition{widx,1} = condition;
        wMidline(widx,1) = midline;
        
     % Wave isolation + calculate onset, end, and duration
        tempWaveData = waveData; 
        tempWaveData(tempWaveData~=waveID) = 0; % isolate current wave
        tempWaveData(tempWaveData>0) = 1; % make all active vox=1 instead of wave ID
        sumFrames = sum(tempWaveData, [1, 2]); % sum along 1st and 2nd dims to make it 1D going in t dim 
        activeFrames = find(sumFrames>0);
        wStart = activeFrames(1); 
        wEnd = activeFrames(end);
        wDur = length(activeFrames); 
        waveMat = tempWaveData(:,:,wStart:wEnd); % isolated wave matrix
        wVox = nnz(waveMat); % number of active voxels
        
        % store data
        wStartFrame(widx,1) = wStart;
        wEndFrame(widx,1) = wEnd;
        wDuration(widx,1) = wDur;
        waveMatrix{widx,1} = waveMat;
        wVoxels(widx,1) = wVox;

    % Area
        wAreaMat = sum(waveMat,3); % make wave matrix 2D
        wAreaMatrix{widx,1} = wAreaMat(:,:); 
        wPix = nnz(wAreaMat); % number of active pixels
        wPixels(widx,1) = wPix; 
        wA = (wPix*pixelArea)/1000000; % convert to microns², then to mm²
        wArea(widx,1) = wA;
        wRetAreaRatio(widx,1) = wA/totRetArea; % proportion of retina that wave area covers

    % Wave initiation site
        firstFrame = waveMat(:,:,1); 
        [yi,xi] = find(firstFrame); % find coordinates of active pixels
        initCoords = [mean(xi), mean(yi)]; % calc centroid
        wInitCoords{widx,1} = initCoords;

        if initCoords(1,1) < midline
            wInitRegion{widx,1} = "Temporal";
        elseif initCoords(1,1) > midline
            wInitRegion{widx,1} = "Nasal";
        else
            wInitRegion{widx,1} = "Midline";
        end
    
    % Vector Fields
        % preallocate centroid matrices and pixel windows
        centersMatX = nan(ny,nx,wDur);
        centersMatY = nan(ny,nx,wDur);

        xWindows = cell(nx,1);
        for x = 1:nx
            xWindows{x} = max(1,x-kernelSize):min(nx,x+kernelSize);
        end
        
        yWindows = cell(ny,1);
        for y = 1:ny
            yWindows{y} = max(1,y-kernelSize):min(ny,y+kernelSize);
        end
        
        % calculate centroids for each pixel of wave over time
        for t = 1:wDur
            for x = 1:nx
                for y = 1:ny

                    xWin = xWindows{x};
                    yWin = yWindows{y};

                    kernel = waveMat(yWin,xWin,t);

                    [X,Y] = meshgrid(xWin,yWin);

                    activeMask = kernel>0;        
                    N = nnz(activeMask);

                    if N==0
                        centersMatX(y,x,t)=nan;
                        centersMatY(y,x,t)=nan;
                    else
                        centersMatX(y,x,t)=sum(X(activeMask))/N;
                        centersMatY(y,x,t)=sum(Y(activeMask))/N;
                    end

                end
            end
        end
        % calculate diff for each pixel of wave over time
        vx3D = diff(centersMatX,[],3);
        vy3D = diff(centersMatY,[],3);

        % make vectors of each pixel 2D
        vx2D = flipud(sum(vx3D,3,"omitnan")); % flip both x and y matrices up and down so coordinates align with image
        vy2D = flipud(sum(vy3D,3,"omitnan"))*-1; % inverse only y matrix so vectors are aligned with image
        
        % make array of vectors and sum them
        vx1D = vx2D(~isnan(vx2D));
        vy1D = vy2D(~isnan(vy2D));

        vxSum = sum(vx1D);
        vySum = sum(vy1D);

        % store vector data into results variables
        tempVecMatX(:,:,currentWave) = vx2D;
        tempVecMatY(:,:,currentWave) = vy2D;
        wxVectors{widx,1} = vx2D(:,:);
        wyVectors{widx,1} = vy2D(:,:);
        wxVectorSum(widx,1) = vxSum;
        wyVectorSum(widx,1) = vySum;

        % convert vectors into polar coordinates and store them
        [wDirTheta(widx,1),wDirRho(widx,1)] = cart2pol(vxSum,vySum);
        
    % Speed
        % diff in magnitude between each frame per pixel (isolates wave fronts, since they are moving)
        vxwf = (diff(vx3D,[],3))*pixelWidth; % x (in µm)
        vywf = (diff(vy3D,[],3))*pixelHeight; % y (in µm)
        wfSpds = hypot(vxwf,vywf); % hypot = wave front magnitude
        wfSpds(wfSpds==0) = nan; 

        % calculate speeds
        avgSpdPerFrame = mean(wfSpds,[1 2],"omitnan"); % avg speed of each frame in microns
        topSpdPerFrame = max(wfSpds,[],[1 2],"omitnan"); % top speed of each frame in microns

        wSpd = (mean(avgSpdPerFrame,"omitnan")); % average speed of wave
        wTopSpd = (mean(topSpdPerFrame,"omitnan")); % top speed of wave

        % store speed results
        wSpeed(widx,1) = wSpd;
        wTopSpeed(widx,1) = wTopSpd;

    % Distance traveled (based on wave speed)
        wDist = wSpd*wDur; % (µm/s)*s = µm
        wDistTraveled(widx,1) = wDist; % store distance traveled data
        
    % Local Synchrony
        timeWindow = 3;
        waveFrontVar = 1;
        tidx = 1;
        for waveTimePt = wStart+timeWindow:wEnd-timeWindow
        
            currentFrame = tempWaveData(:,:,waveTimePt);
            pastFrame = max(tempWaveData(:,:,waveTimePt-timeWindow:waveTimePt-2), [], 3);
            waveFrontFrame = currentFrame-pastFrame;
            waveFrontFrame(waveFrontFrame>0) = 1;
            waveFrontFrame(waveFrontFrame<1) = 0;
            futureFrame = max(tempWaveData(:,:,waveTimePt+2:waveTimePt+timeWindow), [], 3);
            futureFrameOnly = futureFrame-currentFrame;
            futureFrameOnly(futureFrameOnly>0) = 1;
            futureFrameOnly(futureFrameOnly<1) = 0;
        
            if nnz(futureFrameOnly) == 0 
                continue
            end
        
            if waveFrontVar == 1
                currentFrame = waveFrontFrame;
            end
        
            currentFrame(currentFrame == 0) = nan;
            actInWave = spikeRateMovie(:,:,waveTimePt).*currentFrame;
            meanActInWave = mean(actInWave,'all','omitnan');
            
            futureFrameOnly(futureFrameOnly == 0) = nan;
            actOutWave = spikeRateMovie(:,:,waveTimePt).*futureFrameOnly;
            meanActOutWave = mean(actOutWave,'all','omitnan');
        
            % store the results for the current time point 
            wActRatios(tidx,1) = (meanActInWave-meanActOutWave)/(meanActInWave+meanActOutWave);   

            tidx = tidx+1;
        end

        wLocalSynch(widx,1) = mean(wActRatios,'omitnan'); % store local synch data

    % add a row to AIF var
        wAIF(widx,1) = nan;

    % go to next wave index 
        widx = widx+1;

    end
  
% calculate IWI and frequency
    pixIWI = nan(ny,nx);
    emptyPix = nan(ny,nx); % needed to use as a filter later in STTC analysis
    for xx = 1:nx
        for yy = 1:ny
            
            currentPix = waveData(yy,xx,:);
            pixUniqueIDs = unique(currentPix(~isnan(currentPix)));            

            starts = nan(length(pixUniqueIDs),1);
            for u = 1:length(pixUniqueIDs) % get the start frame for each wave that passes through current pixel
                tempMat = currentPix;
                tempMat(tempMat~=pixUniqueIDs(u)) = 0;
                tempMat(tempMat>0) = 1;
                if tempMat(:,:,1) == 1
                    starts(u) = 1;
                else
                    diffs = diff(tempMat);
                    starts(u) = find(diffs==1,1,"first");
                end
            end
            
            starts = sort(starts); % sort so the diff is calculated correctly
    
            iwi = diff(starts); % calc IWI
            iwi = abs(iwi);
            avgPixIWI = mean(iwi,"omitnan");
            pixIWI(yy,xx) = avgPixIWI;

            if isnan(avgPixIWI)
                emptyPix(yy,xx) = 1; % identify pixels that have no waves (to later filter from STTC analysis)
            else
                emptyPix(yy,xx) = 0;
            end

        end
    end
    
    pixIWI = reshape(pixIWI,[ny*nx,1]);
    pixIWI(isnan(pixIWI)) = []; % get rid of pixels that have no waves
    pixIWI = pixIWI/60; % convert to minutes
    recAvgIWI = mean(pixIWI,"omitnan"); % mean IWI 
    recAvgFreq = 1/recAvgIWI; % calc frequency

    % store IWI and freq
    IWIperROI{currentFile} = pixIWI;
    avgIWI(currentFile) = round(recAvgIWI,4,"significant");    
    frequency(currentFile) = round(recAvgFreq,4,"significant");

% calculate averages for current file
    avgDuration(currentFile) = round(mean(wDuration(wFileNum==currentFile),"omitnan"),4,"significant");
    avgArea(currentFile) = round(mean(wArea(wFileNum==currentFile),"omitnan"),4,"significant"); 
    avgLocalSynch(currentFile) = round(mean(wLocalSynch(wFileNum==currentFile),"omitnan"),4,"significant");
    avgSpeed(currentFile) = round(mean(wSpeed(wFileNum==currentFile),"omitnan"),4,"significant");
    avgTopSpeed(currentFile) = round(mean(wTopSpeed(wFileNum==currentFile),"omitnan"),4,"significant");
    avgDistance(currentFile) = round(mean(wDistTraveled(wFileNum==currentFile),"omitnan"),4,"significant");

% calculate vector sums of all waves for current file
    globalVxSum = sum(tempVecMatX,3);
    globalVySum = sum(tempVecMatY,3);
    flowMatrixX{currentFile} = globalVxSum(:,:);
    flowMatrixY{currentFile} = globalVySum(:,:);

% calculate initiation bias index
    [sideCounts,sides] = groupcounts(string(wInitRegion(wFileNum==currentFile)));

    numN = sideCounts(strcmp(sides,"Nasal"));
        if isempty(numN)
            numN = 0;
        end

    numT = sideCounts(strcmp(sides,"Temporal"));
        if isempty(numT)
            numT = 0;
        end

    IBI(currentFile) = numN/(numN+numT); % value closer to 1 = bias closer to nasal

% calculate nasal bias index (prop bias)
    % convert from radians to degrees
    thetaDeg = rad2deg(wDirTheta(wFileNum==currentFile)); 
    negValuesIdx = find(thetaDeg<0);
    thetaDeg(negValuesIdx) = thetaDeg(negValuesIdx)+360;

    % count num waves in each quadrant
    numNQ = length(find(thetaDeg>=0 & thetaDeg<=90 | thetaDeg>=270));
    numTQ = length(find(thetaDeg>=90 & thetaDeg<=270));

    % calc NBI
    NBI(currentFile) = numNQ./(numNQ+numTQ);

% calculate activity impact factor
    voxTable = table;
    voxTable.waveIdx = find(wFileNum==currentFile);
    voxTable.numVoxels = wVoxels(wFileNum==currentFile);

    sVoxTable = sortrows(voxTable,"numVoxels","ascend");
    sVoxTable.AIF = cumsum(sVoxTable.numVoxels/sum(sVoxTable.numVoxels));
    usVoxTable = sortrows(sVoxTable,"waveIdx","ascend");

    wAIF(wFileNum==currentFile) = usVoxTable.AIF;

    % calc ratio of waves responsible for x% of impact
    gradient = linspace(0,1,11);

    ratioGradient = table;
    ratioGradient.AIFpercentile = round(gradient'*100,1,'decimal');
    
    for g = 1:length(gradient)
        impact = gradient(g);
        ratioGradient.percentWaves(g) = (height(usVoxTable.AIF(usVoxTable.AIF>=impact))/nWaves)*100;
    end

    AIFgradient{currentFile} = ratioGradient;

% calculate STTC
    % make sure emptyPix aligns with mapping coords
    if strcmp(fileTable.eye(currentFile),'R') 
        emptyPix = fliplr(emptyPix);
    end

    % filter out channels that have no neuron (no waves)
    emptyCh = nan(numChannels,1);
    for n = 1:numChannels
        currentX = mappingData.x(n);
        currentY = mappingData.y(n);

        xInd = find(currentX == uniqueXs);
        yInd = find(currentY == uniqueYs);

        if emptyPix(yInd,xInd)==1
            emptyCh(n) = channelIDs(n);
        end
    end
    emptyCh(isnan(emptyCh)) = [];
    spikeData.spikeTime(ismember(spikeData.channel,emptyCh)) = [];
    spikeData.channel(ismember(spikeData.channel,emptyCh)) = [];
    channelIDs(ismember(channelIDs,emptyCh)) = [];
    numChannels = length(channelIDs) ;

    % organize spike trains for each channel
    [G,channelList] = findgroups(spikeData.channel);
    chSpikeTimes = splitapply(@(x){x},spikeData.spikeTime,G);

    % get all possible channel combinations
    pairIdx = nchoosek(1:numChannels,2);
    numChComb = size(pairIdx,1);

    % find distances b/w pairs
    [~,mapIdx] = ismember(channelIDs,mappingData.channel);

    dx = mappingData.x(mapIdx(pairIdx(:,2))) - mappingData.x(mapIdx(pairIdx(:,1)));
    dy = mappingData.y(mapIdx(pairIdx(:,2))) - mappingData.y(mapIdx(pairIdx(:,1)));

    distances = hypot(dx,dy);

    % calculate TA/TB
    T = zeros(numChannels,1);

    for n = 1:numChannels
        spikes = chSpikeTimes{n};
        
        % create intervals
        dtStarts = max(0, spikes-dt);
        dtEnds   = min(timeVar, spikes+dt);

        % calc unique time spent within spike
        currentStart = dtStarts(1);
        currentEnd   = dtEnds(1);
        timeWithinSpike = 0;
    
        for s = 2:length(dtStarts)  
            if dtStarts(s) <= currentEnd % if overlapping interval
                currentEnd = max(currentEnd, dtEnds(s));
            else
                timeWithinSpike = timeWithinSpike + (currentEnd-currentStart);

                % go to next interval
                currentStart = dtStarts(s);
                currentEnd = dtEnds(s);
            end
        end

        % add last interval
        timeWithinSpike = timeWithinSpike + (currentEnd-currentStart);
        
        % calc T for current channel
        T(n) = timeWithinSpike/timeVar;

    end

    % calc STTC on all combinations
    sttc = nan(numChComb,1);

    parfor n = 1:numChComb

        chA = pairIdx(n,1);
        chB = pairIdx(n,2);

        spikesA = chSpikeTimes{chA};
        spikesB = chSpikeTimes{chB};

        TA = T(chA);
        TB = T(chB);

        % calculate PA
        pointer = 1;
        matches = 0;
        
        for a = 1:length(spikesA)
  
            % advance the pointer past spikes that could never match
            while pointer < length(spikesB) && spikesB(pointer) < spikesA(a)-dt
                pointer = pointer + 1;
            end
        
            % stop when a B spike is approaching A spike
            if abs(spikesB(pointer)-spikesA(a)) <= dt % check if spikes are within +/- dt of each other
                matches = matches + 1;       
            end
        
        end
        
        PA = matches/length(spikesA);

        % calculate PB
        pointer = 1;
        matches = 0;
        
        for b = 1:length(spikesB)
            
            while pointer < length(spikesA) && spikesA(pointer) < spikesB(b)-dt
                pointer = pointer + 1;
            end
        
            if abs(spikesA(pointer)-spikesB(b)) <= dt
                matches = matches + 1;       
            end
        
        end
        
        PB = matches/length(spikesB);

        % calc STTC for current pair
        sttc(n) = 0.5*((PA-TB)/(1-PA*TB)+(PB-TA)/(1-PB*TA));

    end

    % filter out distances <1000µm
    sttc(distances>1000) = [];
    distances(distances>1000) = [];

    % store data
    STTC{currentFile} = sttc;
    channelDistances{currentFile} = distances;

end

% display how long it took to analyze current file
time_to_analyze = toc

%% Store results into output variables

% convert word variables to string format
wSex = string(wSex);
wEye = string(wEye);
wGenotype = string(wGenotype);
wCondition = string(wCondition);
wInitRegion = string(wInitRegion);

if tableFormat == 1

    % per file/recording table
    perFileResults = table(fileNum,expNums,ages,sexes,eyes,genotypes,conditions,midlines,timeVars,numWaves,retinalArea, ...
        AIFgradient,IWIperROI,avgIWI,frequency,avgDuration,avgArea,avgLocalSynch,avgSpeed,avgTopSpeed,avgDistance, ...
        IBI,NBI,flowMatrixX,flowMatrixY,STTC,channelDistances, ...
        'VariableNames',{'fileNumber','expNumber','age','sex','eye','genotype','condition','midline','recLength','numWaves','retinalArea', ...
        'AIFgradient','IWIperROI','avgIWI','frequency','avgDuration','avgArea','avgLSI','avgSpeed','avgTopSpeed','avgDistanceTraveled', ...
        'IBI','NBI','flowMatrixX','flowMatrixY','STTC','channelDistances'});

    % per wave table
    perWaveResults = table(wFileNum,wExpNum,wAge,wSex,wEye,wGenotype,wCondition,wMidline,waveMatrix,wAIF,wVoxels, ...
        wStartFrame,wEndFrame,wDuration,wInitCoords,wInitRegion,wAreaMatrix,wPixels,wArea,wRetAreaRatio,wLocalSynch, ...
        wSpeed,wTopSpeed,wDistTraveled,wxVectors,wyVectors,wxVectorSum,wyVectorSum,wDirTheta,wDirRho, ...
        'VariableNames',{'fileNumber','expNumber','age','sex','eye','genotype','condition','midline','waveMatrix','AIF','numVoxels', ...
        'startFrame','endFrame','duration','initiationCoordinates','initiationRegion','areaMatrix','numPixels','area','waveToRetinaAreaRatio','LSI', ...
        'avgSpeed','topSpeed','distanceTraveled','waveFlowVectorsX','waveFlowVectorsY','flowVectorSumX','flowVectorSumY','directionTheta','directionRho'});
else
    % store as a structure
    waveAnalysisResults = struct;
    waveAnalysisResults.perFile = table(fileNum,expNums,ages,sexes,eyes,genotypes,conditions,midlines,timeVars,numWaves,retinalArea, ...
        AIFgradient,IWIperROI,avgIWI,frequency,avgDuration,avgArea,avgLocalSynch,avgSpeed,avgTopSpeed,avgDistance, ...
        IBI,NBI,flowMatrixX,flowMatrixY,STTC,channelDistances, ...
        'VariableNames',{'fileNumber','expNumber','age','sex','eye','genotype','condition','midline','recLength','numWaves','retinalArea', ...
        'AIFgradient','IWIperROI','avgIWI','frequency','avgDuration','avgArea','avgLSI','avgSpeed','avgTopSpeed','avgDistanceTraveled', ...
        'IBI','NBI','flowMatrixX','flowMatrixY','STTC','channelDistances'});
    waveAnalysisResults.perWave = struct;
        waveAnalysisResults.perWave.metaData = table(wFileNum,wExpNum,wAge,wSex,wEye,wGenotype,wCondition,wMidline, ...
            'VariableNames',{'fileNumber','expNumber','age','sex','eye','genotype','condition','midline'});      
        waveAnalysisResults.perWave.waveMatrix = waveMatrix;
        waveAnalysisResults.perWave.AIF = wAIF;
        waveAnalysisResults.perWave.numVoxels = wVoxels;
        waveAnalysisResults.perWave.startFrame = wStartFrame;
        waveAnalysisResults.perWave.endFrame = wEndFrame;
        waveAnalysisResults.perWave.duration = wDuration;
        waveAnalysisResults.perWave.initiationCoordinates = wInitCoords;
        waveAnalysisResults.perWave.initiationRegion = wInitRegion;
        waveAnalysisResults.perWave.areaMatrix = wAreaMatrix;
        waveAnalysisResults.perWave.numPixels = wPixels;
        waveAnalysisResults.perWave.area = wArea;
        waveAnalysisResults.perWave.waveToRetinaAreaRatio = wRetAreaRatio;
        waveAnalysisResults.perWave.LSI = wLocalSynch;
        waveAnalysisResults.perWave.avgSpeed = wSpeed;
        waveAnalysisResults.perWave.topSpeed = wTopSpeed;
        waveAnalysisResults.perWave.distanceTraveled = wDistTraveled;
        waveAnalysisResults.perWave.waveFlowVectorsX = wxVectors;
        waveAnalysisResults.perWave.waveFlowVectorsY = wyVectors;
        waveAnalysisResults.perWave.flowVectorSumX = wxVectorSum;
        waveAnalysisResults.perWave.flowVectorSumY = wyVectorSum;
        waveAnalysisResults.perWave.directionTheta = wDirTheta;
        waveAnalysisResults.perWave.directionRho = wDirRho;
end

%% Save results variables
   
if saveResults == 1
    cd(saveResultsLocation)
    save(todaysDate+"_waveAnalysisResults.mat","perWaveResults","perFileResults")
    % save(todaysDate+"_waveAnalysisResults.mat","waveAnalysisResults")
    % uisave({'perWaveResults','perFileResults'},'waveAnalysisResults')
end
