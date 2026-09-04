%% Wave analysis pipeline for calcium imaging data 
% By Kaylee Odum
% Created 260430

% v2 - 260526 - KO - changed output variable from table to structure 

clearvars, close all

%% User inputs 

todaysDate = "260623";
fileTableName = "VuofoPaperData.xlsx";

analyzeFiltData = 1; % toggle to analyze wave data with median filter applied (standard/default); set to 0 if you want to analyze raw wave data
pixelWidth = 33.1; % width of pixels in µm
pixelHeight = 35.61; % height of pixels in µm
pixelArea = pixelWidth*pixelHeight; % area of pixels in µm²
kernelSize = 1;
saveResults = 1; % 1 if you want to save the results
tableFormat = 1; % 1 if you want the results var to be a table, 0 if you want a structure (structure = smaller file size)

% mac OS directories
fileTablePath = "/Users/kayleeodum/Library/CloudStorage/OneDrive-Vanderbilt/Tiriac Lab/Toolbox Project/figure5_calciumImaging/2_code";
rawDataPath = "/Users/kayleeodum/Library/CloudStorage/OneDrive-Vanderbilt/Tiriac Lab/Toolbox Project/figure5_calciumImaging/1_data";
waveVarPath = "/Users/kayleeodum/Library/CloudStorage/OneDrive-Vanderbilt/Tiriac Lab/Toolbox Project/figure5_calciumImaging/3_waveDetectionResults/1_waveData_mat";
saveResultsLocation = "/Users/kayleeodum/Library/CloudStorage/OneDrive-Vanderbilt/Tiriac Lab/Toolbox Project/figure5_calciumImaging/4_waveAnalysisResults"; 

% windows OS directories
% fileTablePath = "C:\Users\yeagerkm\OneDrive - Vanderbilt\Tiriac Lab\Toolbox Project\figure5_calciumImaging\2_code";
% rawDataPath = "C:\Users\yeagerkm\OneDrive - Vanderbilt\Tiriac Lab\Toolbox Project\figure5_calciumImaging\1_data";
% waveVarPath = "C:\Users\yeagerkm\OneDrive - Vanderbilt\Tiriac Lab\Toolbox Project\figure5_calciumImaging\3_waveDetectionResults\1_waveData_mat";
% saveResultsLocation = "C:\Users\yeagerkm\OneDrive - Vanderbilt\Tiriac Lab\Toolbox Project\figure5_calciumImaging\4_waveAnalysisResults"; 

%% Load the file table

% load meta data table
cd(fileTablePath)    
fileTable = readtable(fileTableName);

% extract file table columns into their own array variables
expNums = fileTable.expNum; % experiment number
sampleRates = fileTable.Hz;  % in Hz
fileNames = string(fileTable.fileName);
waveDataVarNames = string(fileTable.waveDataVariable);

%% Run wave analysis pipeline

numFiles = size(fileTable,1);

% initialize per file results vars
fileNum = nan(numFiles,1);
timeVars = nan(numFiles,1); % seconds
numWaves = nan(numFiles,1);
retinalArea = nan(numFiles,1);
IWIperROI = cell(numFiles,1);
avgIWI = nan(numFiles,1); % min/wave
frequency = nan(numFiles,1); % waves/min
avgDuration = nan(numFiles,1); % seconds
avgArea = nan(numFiles,1); % mm²
avgLocalSynch = nan(numFiles,1);
avgSpeed = nan(numFiles,1); % µm/s
avgTopSpeed = nan(numFiles,1); % µm/s
avgDistance = nan(numFiles,1); % µm
flowMatrixX = cell(numFiles,1);
flowMatrixY = cell(numFiles,1);

% initialize per wave results vars
wFileNum = [];
wExpNum = [];
wStartFrame = [];
wEndFrame = [];
wDuration = []; % seconds
wAreaMatrix = {};
wVoxels = [];
wPixels = [];
wArea = []; % mm²
wRetAreaRatio = [];
wInitCoords = {};
wxVectors = {};
wyVectors = {};
wxVectorSum = [];
wyVectorSum = [];
wDirTheta = [];
wDirRho = [];
wAvgSpeed = []; % microns/sec
wTopSpeed = []; 
wDistTraveled = [];
wLocalSynch = [];

% initialize index var for all wave vars
widx = 1;

%%%%% Analyze each file %%%%%
for currentFile = 1:numFiles
    currentFile
    tic
 
% load raw data (for local synch)
    cd(rawDataPath)
    rawMovie = tiffreadVolume(fileNames(currentFile));

% load wave data variables 
    cd(waveVarPath) % set path
    load(waveDataVarNames(currentFile));
    
    if analyzeFiltData==1
        waveData = waveData_mfilt; % will use filtered wave data
        clear waveData_mfilt
    else
        clear waveData_mfilt % will use raw wave data (which is already named waveData when loaded)
    end

    % make sure empty pixels are nan, not 0
    waveData(waveData==0) = nan;

% extract meta data variables for current file
    expNum = expNums(currentFile);
    sampleRate = sampleRates(currentFile);

% create ID and size variables
    waveIDs = unique(waveData(~isnan(waveData)));
    nWaves = length(waveIDs);
    numFrames = size(waveData,3);    
    timeVar = numFrames/sampleRate;
    [ny,nx,nt] = size(waveData);

% calculate area of retina covered by waves
    coverageTempMat = sum(waveData,3,"omitnan"); % make it 2D
    coverageTempMat(coverageTempMat>0) = 1; % make each active pixel = 1
    totRetPixels = sum(coverageTempMat(coverageTempMat==1));
    totRetArea = (totRetPixels*pixelArea)/1000000; % convert to microns², then to mm²

% preallocate vector matrix for later calculating sums 
    tempVecMatX = nan(size(waveData,1),size(waveData,2),nWaves);
    tempVecMatY = nan(size(waveData,1),size(waveData,2),nWaves);

% store wave summary meta data for current file
    fileNum(currentFile) = currentFile;
    timeVars(currentFile) = timeVar;
    numWaves(currentFile) = nWaves;
    retinalArea(currentFile) = totRetArea;

    %%%%% Analyze each wave %%%%%
    for currentWave = 1:nWaves
        waveID = waveIDs(currentWave);

        % store meta data about the recording that the wave belongs to
        wFileNum(widx,1) = currentFile;
        wExpNum(widx,1) =  expNum;
        
    % Wave duration, onset and end
        tempWaveData = waveData; 
        tempWaveData(tempWaveData ~= waveID) = 0; % isolate current wave
        sumFrames = sum(tempWaveData, [1, 2]); % sum along 1st and 2nd dims to make it 1D going in t dim 
        activeFrames = find(sumFrames>0);
        wStart = activeFrames(1); 
        wEnd = activeFrames(end);
        wNumFrames = length(activeFrames);
        wDur = (wNumFrames/sampleRate); 
        waveMat = tempWaveData(:,:,wStart:wEnd); % isolated wave matrix
        wVox = nnz(waveMat); % number of active voxels
        % store data
        wStartFrame(widx,1) = wStart;
        wEndFrame(widx,1) = wEnd;
        wDuration(widx,1) = wDur;
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

    % Vector Fields
        % preallocate centroid matrices and pixel windows
        centersMatX = nan(ny,nx,wNumFrames);
        centersMatY = nan(ny,nx,wNumFrames);

        xWindows = cell(nx,1);
        for x = 1:nx
            xWindows{x} = max(1,x-kernelSize):min(nx,x+kernelSize);
        end
        
        yWindows = cell(ny,1);
        for y = 1:ny
            yWindows{y} = max(1,y-kernelSize):min(ny,y+kernelSize);
        end
        
        % calculate centroids for each pixel of wave over time
        for t = 1:wNumFrames
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
        wAvgSpeed(widx,1) = wSpd;
        wTopSpeed(widx,1) = wTopSpd;

    % Distance traveled (based on wave speed)
        wDist = wSpd*wDur; % (µm/s)*s = µm
        wDistTraveled(widx,1) = wDist; % store distance traveled data
    
    % Local Synchrony
        timeWindow = 3;
        waveFrontVar = 1;
        wActRatios = nan(length(activeFrames)-timeWindow,1);
        
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
        
            if sum(futureFrameOnly,3) == 0
                continue
            end
        
            if waveFrontVar == 1
                currentFrame = waveFrontFrame;
            end
        
            currentFrame(currentFrame == 0) = nan;
            realActInWave = rawMovie(:,:,waveTimePt).*currentFrame;
            meanActInWave = mean(realActInWave,'all','omitnan');
            
            futureFrameOnly(futureFrameOnly == 0) = nan;
            realActOutWave = rawMovie(:,:,waveTimePt).*futureFrameOnly;
            meanActOutWave = mean(realActOutWave,'all','omitnan');
        
            % store the results for the current wave 
            wActRatios(waveTimePt-wStart+1,1) = (meanActInWave-meanActOutWave)/(meanActInWave+meanActOutWave);
        end

        wLocalSynch(widx,1) = mean(wActRatios,'omitnan');

    % go to next wave index
        widx = widx+1;

    end % end of per wave loop
    
% calculate IWI and frequency
    pixIWI = nan(ny,nx);
    for xx = 1:nx
        for yy = 1:ny

            currentCoords = waveData(yy,xx,:);
            currUniqueIDs = unique(currentCoords(~isnan(currentCoords)));            

            starts = nan(length(currUniqueIDs),1);
            for u = 1:length(currUniqueIDs) % get the start frame for each wave that passes through current pixel
                tempMat = currentCoords;
                tempMat(tempMat~=currUniqueIDs(u)) = 0;
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
            avgPixIWI = mean(iwi);
            pixIWI(yy,xx) = avgPixIWI;

        end
    end
    
    pixIWI = reshape(pixIWI,[ny*nx,1]);
    pixIWI(isnan(pixIWI)) = []; % get rid of pixels that have no waves 
    pixIWI = pixIWI/(sampleRate*60); % convert to minutes
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
    avgSpeed(currentFile) = round(mean(wAvgSpeed(wFileNum==currentFile),"omitnan"),4,"significant");
    avgTopSpeed(currentFile) = round(mean(wTopSpeed(wFileNum==currentFile),"omitnan"),4,"significant");
    avgDistance(currentFile) = round(mean(wDistTraveled(wFileNum==currentFile),"omitnan"),4,"significant");

% calculate vector sums of all waves for current file
    globalVxSum = sum(tempVecMatX,3);
    globalVySum = sum(tempVecMatY,3);
    flowMatrixX{currentFile} = globalVxSum(:,:);
    flowMatrixY{currentFile} = globalVySum(:,:);
% display how long it took to analyze current file
    time_to_analyze_file = toc

end % end of per file loop

%% store results into output vars


if tableFormat == 1 % store as tables

    % per file/recording
    perFileResults = table(fileNum,expNums,retinalArea,timeVars, ...
        numWaves,IWIperROI,avgIWI,frequency,avgDuration,avgArea,avgLocalSynch,avgSpeed,avgTopSpeed,avgDistance,flowMatrixX,flowMatrixY, ...
        'VariableNames',{'fileNumber','expNumber','retinalArea','recLength', ...
        'numWaves','IWIperROI','avgIWI','frequency','avgDuration','avgArea','avgLocalSynchrony','avgSpeed','avgTopSpeed','avgDistance','flowMatrixX','flowMatrixY'});

    % per wave 
    perWaveResults = table(wFileNum,wExpNum,wStartFrame,wEndFrame,wDuration, ...
        wInitCoords,wAreaMatrix,wVoxels,wPixels,wArea,wRetAreaRatio,wLocalSynch, ...
        wAvgSpeed,wTopSpeed,wDistTraveled,wxVectors,wyVectors,wxVectorSum,wyVectorSum,wDirTheta,wDirRho, ...
        'VariableNames',{'fileNumber','expNumber','startFrame','endFrame','duration', ...
        'initiationCoordinates','areaMatrix','numVoxels','numPixels','area','waveToRetinaAreaRatio','localSynchrony', ...
        'avgSpeed','topSpeed','distanceTraveled','waveFlowVectorsX','waveFlowVectorsY','flowVectorSumX','flowVectorSumY','directionTheta','directionRho'});

else % store as a structure
    waveAnalysisResults = struct;
    waveAnalysisResults.perFile = table(fileNum,expNums,retinalArea,timeVars, ...
        numWaves,IWIperROI,avgIWI,frequency,avgDuration,avgArea,avgLocalSynch,avgSpeed,avgTopSpeed,avgDistance,flowMatrixX,flowMatrixY, ...
        'VariableNames',{'fileNumber','expNumber','retinalArea','recLength', ...
        'numWaves','IWIperROI','avgIWI','frequency','avgDuration','avgArea','avgLocalSynchrony','avgSpeed','avgTopSpeed','avgDistance','flowMatrixX','flowMatrixY'});
    waveAnalysisResults.perWave = struct;  
        waveAnalysisResults.perWave.fileNumber = wFileNum;
        waveAnalysisResults.perWave.expNumber = wExpNum;
        waveAnalysisResults.perWave.startFrame = wStartFrame;
        waveAnalysisResults.perWave.endFrame = wEndFrame;
        waveAnalysisResults.perWave.duration = wDuration;
        waveAnalysisResults.perWave.initiationCoordinates = wInitCoords;
        waveAnalysisResults.perWave.areaMatrix = wAreaMatrix;
        waveAnalysisResults.perWave.numVoxels = wVoxels;
        waveAnalysisResults.perWave.numPixels = wPixels;
        waveAnalysisResults.perWave.area = wArea;
        waveAnalysisResults.perWave.waveToRetinaAreaRatio = wRetAreaRatio;
        waveAnalysisResults.perWave.localSynchrony = wLocalSynch;
        waveAnalysisResults.perWave.avgSpeed = wAvgSpeed;
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
    % save(todaysDate+"_waveAnalysisResults.mat","perWaveResults","perFileResults")
    save(todaysDate+"_waveAnalysisResults.mat","waveAnalysisResults")
end
