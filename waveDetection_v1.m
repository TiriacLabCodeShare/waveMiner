%% Wave detection code
% By Alex Tiriac and Kaylee Odum
%
% 
% v3 - first functional version, can ID waves
% v4 - added a per neuron burst ID threshold (will work on hist data)
% v5 - Implemented Arda's dynamic burst det code
% v6 - Alex major update - modified Arda's burst det code
% v7 - Dynamic baseline over the course of the recording
% v8&9 - Enriching signal + various debugging
% v10 - tiff export


clearvars, close all


%% User inputs

%%% Fixed parameters %%%

fileName =  '250708_P9_M_L_Baseline.h5'; % .h5 file name

% Directories for uploading and saving files
% filePath = "C:\Users\yeagerkm\OneDrive - Vanderbilt\Tiriac Lab\Toolbox Project\3_wavesAcrossDevelopment\1_data";
% saveWaveDataPath = "C:\Users\yeagerkm\OneDrive - Vanderbilt\Tiriac Lab\Toolbox Project\3_wavesAcrossDevelopment\3_waveDetectionResults\1_waveData_mat";
% saveWaveTiffPath = "C:\Users\yeagerkm\OneDrive - Vanderbilt\Tiriac Lab\Toolbox Project\3_wavesAcrossDevelopment\3_waveDetectionResults\2_waveData_tif";

filePath = "/Users/kayleeodum/Library/CloudStorage/OneDrive-Vanderbilt/Tiriac Lab/Toolbox Project/3_wavesAcrossDevelopment/1_data";
saveWaveDataPath = "/Users/kayleeodum/Library/CloudStorage/OneDrive-Vanderbilt/Tiriac Lab/Toolbox Project/3_wavesAcrossDevelopment/3_waveDetectionResults/1_waveData_mat";
saveWaveTiffPath = "/Users/kayleeodum/Library/CloudStorage/OneDrive-Vanderbilt/Tiriac Lab/Toolbox Project/3_wavesAcrossDevelopment/3_waveDetectionResults/2_waveData_tif";

% Recording parameters
timeVar = 3600; % length of recording in seconds
sampleRate = 20000;  % sample rate in Hz


%%% Adjusted parameters %%%

% Dynamic baseline over time
baselineSampling = 10; % calc baseline every xth fraction of rec length

% Burst variables
isiThreshold = 1; % in seconds
isiThresholdDynamic = 4; % 4x increase in spike rate (or shortening of ISI)

% Burst mean filter for low-activity (enriching signal)
meanFiltStep = 0;  % 0 for off, 1 for on. Turn on for low-activity data

% Wave Detection variables
isPartOfWave = 7; % threshold for how many neighbors must have act to be part of a wave


%%% Other parameters %%%
playMovieToggle = 1; % set to 1 if you want to see a movie of waveID
movieSpeed = 0.01; % 1 = 1 frame/s; smaller # = faster speed

%% Loading and preparing the data

cd(filePath)
data = h5read(fileName, '/data_store/data0000/spikes');
mapping = h5read(fileName,'/data_store/data0000/settings/mapping');

% replace frameno with spikeTime (it's more intuitive)
[data.('spikeTime')] = data.('frameno');
data = rmfield(data,'frameno');

% start spikeTime at zero and convert into seconds
data.spikeTime = double(data.spikeTime - min(data.spikeTime))/sampleRate;

% for long data, only analyze up to timeVar
data.channel(data.spikeTime>timeVar) = [];
data.amplitude(data.spikeTime>timeVar) = [];
data.spikeTime(data.spikeTime>timeVar) = [];

% calc whole spikeRate for plot
wholeRetSpikeRate = histcounts(data.spikeTime,0:1:timeVar);

fig0 = figure('Name','stable baseline check');
subplot(3,1,1)
plot(wholeRetSpikeRate,'k')
xlim([0 timeVar])
subplot(3,1,2:3)
scatter(data.spikeTime,data.channel,0.2,"k.")
xlim([0 timeVar])

%% Make sure data is aligned to mapping

% synch data with mapping
uniqueChannels = unique(data.channel);
[uniqueChannels,ia,ib] = intersect(uniqueChannels,mapping.channel);
mapping.channel = mapping.channel(ib);
mapping.electrode = mapping.electrode(ib);
mapping.x = mapping.x(ib);
mapping.y = mapping.y(ib);

% get rid of non mapping channels from data
numUniqueChannelsInData = unique(data.channel);
for n = 1:length(numUniqueChannelsInData)
    if ismember(numUniqueChannelsInData(n), uniqueChannels) == 0
        chInd = find(data.channel == numUniqueChannelsInData(n));
        data.amplitude(chInd) = [];
        data.channel(chInd) = [];
        data.spikeTime(chInd) = [];
    end
end

uniqueChannels = unique(data.channel);

numChannels = length(uniqueChannels);

%% Identify baseline and ret wave periods

% calculate spikeRate per channel
spikeRateChannel = zeros(numChannels,timeVar);
for n = 1:numChannels
    currentChannel = uniqueChannels(n);
    spikeTimesCurrCh = data.spikeTime(data.channel==currentChannel);
    spikeRateChannel(n,:) = histcounts(spikeTimesCurrCh, 0:1:timeVar);
end 

wholeRetSpikeRate = histcounts(data.spikeTime,0:1:timeVar);
baselineInd = nan(10,baselineSampling);
meanBaselineFR = nan(numChannels,baselineSampling);
meanBaselineISI = nan(numChannels,baselineSampling);

bursts.time = [];
bursts.channel = [];

fig1 = figure('Name','Data and baseline checks');
sgtitle("Data and Baseline Checks")
subplot(3,1,1)
plot(wholeRetSpikeRate,'k')
xlim([0 timeVar])
title("raw spike rate")
subplot(3,1,2:3)
scatter(data.spikeTime,data.channel,0.2,"k.")
ylim([0 numChannels])
xlim([0 timeVar])
title("raw spikes/channel")
subplot(3,1,1)
hold on

for n = 1:baselineSampling
    if n == 1
        startTime = 1;
    else
        startTime = timeVar/baselineSampling*(n-1);
    end
    endTime = startTime+(timeVar/baselineSampling);

    [troughs, baselineInd] = findpeaks(wholeRetSpikeRate(startTime:endTime)*(-1),"NPeaks",5,'SortStr','descend');

    baselineInd = baselineInd+(timeVar/baselineSampling)*(n-1);
    baselineInd(baselineInd < 3) = [];
    baselineInd(baselineInd > timeVar-3) = [];

    xline(baselineInd,'r')

    baselineFR = [];
    for m = 1:length(baselineInd)
        baselineFR = [baselineFR, spikeRateChannel(:,baselineInd(m)-2:baselineInd(m)+2)];
    end
    meanBaselineFR(:,n) = mean(baselineFR,2);
    meanBaselineISI(:,n) = 1./mean(baselineFR,2);

        for m = 1:numChannels
        currentChannel = uniqueChannels(m);

        % get spikes corresponding to the current channel
        channelSpikes = data.spikeTime(data.channel == currentChannel & data.spikeTime>startTime & data.spikeTime<endTime );
        
        % calculate inter-spike-interval for this channel
        channelISI = diff(channelSpikes);
        
        % find any ISI that fall within isiThreshold
        if meanBaselineISI(m,n) == inf
            burstIndices = find(channelISI < isiThreshold);
        else
            burstIndices = find(channelISI < meanBaselineISI(m,n)/isiThresholdDynamic);
        end
        
        bursts.time = [bursts.time;channelSpikes(burstIndices)];
        bursts.channel = [bursts.channel;repmat(currentChannel,[length(burstIndices),1])];
    
        end

end


%% Burst detection

% plot to check results of burst detection
fig2 = figure('Name','Data and baseline checks');

subplot(6,1,1)
plot(wholeRetSpikeRate,'k')
title("raw spike rate")

subplot(6,1,2:3)
scatter(data.spikeTime,data.channel,0.2,"k.")
ylim([0 numChannels])
xlim([0 timeVar])
title("raw spikes/channel")

subplot(6,1,4:5)
scatter(bursts.time,bursts.channel,0.2,"k.")
ylim([0 numChannels])
xlim([0 timeVar])
title("burst spikes/channel")

subplot(6,1,6)
plot(histcounts(bursts.time,0:1:timeVar),'k')
title("burst spike rate")


% calculate burstRate per channel
burstRateChannel = zeros(numChannels,timeVar);
for n = 1:numChannels
    currentChannel = uniqueChannels(n);
    burstTimesCurrCh = bursts.time(bursts.channel==currentChannel);
    burstRateChannel(n,:) = histcounts(burstTimesCurrCh, 0:1:timeVar);
end 

% binarize BurstRateChannel
binaryBursts = burstRateChannel > 0;


choice = questdlg('Proceed?','Proceed','Y','N','Y');
if strcmp(choice,'N')
    return
end

%% turn rasterplots into xyt matrix

burstMatrix = binaryBursts;

recTime = size(burstMatrix,2);

uniqueXs = unique(mapping.x);
uniqueYs = unique(mapping.y);
numXs = length(uniqueXs);
numYs = length(uniqueYs);

spikeRateMov = zeros(numYs,numXs,recTime);
burstMov = zeros(numYs,numXs,recTime);

for n = 1:numChannels
    currentX = mapping.x(n);
    currentY = mapping.y(n);

    xInd = find(currentX == uniqueXs);
    yInd = find(currentY == uniqueYs);

    burstMov(yInd,xInd,:) = burstMatrix(n,:);
    spikeRateMov(yInd,xInd,:) = spikeRateChannel(n,:);
end

%% filter and show movies

burstMov_filt = imboxfilt3(burstMov,[3 3 3]);
burstMov_filt = burstMov_filt > 0.12;

if playMovieToggle == 1

    figure
    
        colormap parula
    
    for n = 1:600 % only show first 10 minutes
        subplot(7,1,1:2)      
        imagesc(spikeRateMov(:,:,n))
        clim([0 5]);
        set(gca,'XColor','none','YColor','none')
        title("raw data")
        subplot(7,1,3:4)
        imagesc(burstMov(:,:,n))
        clim([0 1]);
        set(gca,'XColor','none','YColor','none')
        title("burst data")
        subplot(7,1,5:6)
        imagesc(burstMov_filt(:,:,n))
        clim([0 1]);
        set(gca,'XColor','none','YColor','none')
        title("burst data + 3D box filter")
        subplot(7,1,7)
        plot(sum(burstMatrix))
        hold on
        xline(n)
        hold off
        pause(movieSpeed)
    end

end

if meanFiltStep == 1
    burstMov = double(burstMov_filt);
end


%% Wave detection

% inputs
x = ones(3,3,5);
x(:,1,:) = -1;
x(:,2,:) = 0;
x(:,3,:) = 1;
y = ones(3,3,5);
y(1,:,:) = -1;
y(2,:,:) = 0;
y(3,:,:) = 1;
z = ones(3,3,5);
z(:,:,1) = -2;
z(:,:,2) = -1;
z(:,:,3) = 0;
z(:,:,4) = 1;
z(:,:,5) = 2;

winSizes = size(z);
winElements = winSizes(1)*winSizes(2)*winSizes(3);


tic

waveIDcounter = 1; % ID to be given to pixels that are part of the same wave

% pre allocate waveID
waveData = nan(size(burstMov));

wCounterTot = [];

% wave detection
for t = 1:recTime-2 % go through time
    for xx = 2:numXs-1 % go through x dim
        for yy = 2:numYs-1 % go through y dim
            if burstMov(yy,xx,t) == 1 % if pixel has activity
                if isnan(waveData(yy,xx,t)) % and if it is not yet assigned to a wave
                    t
                    listCoords = [yy,xx,t]; % start a list of coordinates starting with current pixel x,y,t
                    wCounter = 1; % while counter that keeps track of how many times we've gone through while loop
                    numUncheckedNb = 1; % check to see how many neighbors left, ends the while loop if this ever reaches 0
                    while numUncheckedNb > 0
                         if listCoords(wCounter,1) < 2 || listCoords(wCounter,1) > numYs-1 || ...
                                 listCoords(wCounter,2) < 2 || listCoords(wCounter,2) > numXs-1 || ...
                                 listCoords(wCounter,3) < 3 || listCoords(wCounter,3) > recTime-2 % check for edge pixels and skips them (could fix this in future version)

                            wCounter = wCounter+1;
                            numUncheckedNb = size(listCoords,1)-wCounter;
                            continue 
                         end
                         tempMat = burstMov(listCoords(wCounter,1)-1:listCoords(wCounter,1)+1,listCoords(wCounter,2)-1:listCoords(wCounter,2)+1,listCoords(wCounter,3)-2:listCoords(wCounter,3)+2); % get a grid of activity centered on current pixel
                         maskMat = sum(tempMat,3) > 0;
                         if sum(maskMat,'all') > isPartOfWave % if sum of activity breaches threshold
                            
                            % Next few lines are matrix math to quickly ID
                            % coordinates of active pixels
                            tempMat(tempMat == 0) = nan;
                            tempX = tempMat.*x;
                            tempX = reshape(tempX,[winElements,1]);
                            tempX(isnan(tempX)) = [];
                            tempY = tempMat.*y;
                            tempY = reshape(tempY,[winElements,1]);
                            tempY(isnan(tempY)) = [];
                            tempZ = tempMat.*z;
                            tempZ = reshape(tempZ,[winElements,1]);
                            tempZ(isnan(tempZ)) = [];
                            
                            listCoords = unique([listCoords;[listCoords(wCounter,1)+tempY, listCoords(wCounter,2)+tempX, listCoords(wCounter,3)+tempZ]],'rows','stable'); % add above coordinates to the ongoing list of active pixels within putative wave
                         end
                         numUncheckedNb = size(listCoords,1)-wCounter; % do we have neighboring pixels that we have not checked yet?
                         wCounter = wCounter+1; % add 1 to the while counter
                    end
                    wCounterTot = [wCounterTot; wCounter];
                    if size(listCoords,1) > 100 % threshold to be big enough to be considered a wave
                        for mm = 1:size(listCoords,1) % this loop assigns a cluster of neighboring pixels the same wave ID
                            waveData(listCoords(mm,1),listCoords(mm,2),listCoords(mm,3)) = waveIDcounter; 
                        end
                        waveIDcounter = waveIDcounter+1; % add 1 to the wave ID
                    end


                end
                
            end


        end
       
    end
end

timeToAnalyze = toc

%% Create list of waveIDs

[numYs, numXs,recTime] = size(waveData);
[waveOccurence, waveID] = groupcounts(reshape(waveData,[numYs*numXs*recTime 1]));
waveOccurence((isnan(waveID))) = [];
waveID((isnan(waveID))) = [];

%% Identify small events

eventSizes = nan(length(waveID),1);
for n = 1:length(waveID)
    currWave = waveID(n);
    tempMat = waveData;
    tempMat(tempMat~=currWave) = nan;
    tempMat(tempMat==currWave) = 1;

    sumMat = sum(tempMat,3,'omitnan');
    areaMat = sumMat > 0;
    eventSizes(n) = sum(areaMat,"all",'omitnan');
end

figure, histogram(eventSizes,100)
title("event size distribution")

%% Filter out small events

waveAreaThreshold = 15; % determine this number based on the natural breakpoint in the wave areas histogram distribution

for n = 1:length(waveID)
    currWave = waveID(n);
    if eventSizes(n) < waveAreaThreshold 
        waveData(waveData == currWave) = nan;
    end
end

%% apply median filter to the data

waveData(isnan(waveData)) = 0;
waveData_mfilt = medfilt3(waveData,[3 3 3]);

%% play waveID movie

if playMovieToggle == 1

    C = colormap(rand(max(waveData_mfilt,[],"all"),3)*0.5+0.5);
    C = [0,0,0;C];
    figure("Name","Movie of Detected Waves","Units", "normalized", "Position",[0.33 0.2 0.33 0.7])
    
    ax(1) = subplot(7,1,1:2);
    colormap(ax(1),"parula")
    clim(ax(1),[0 3]);
    set(ax(1),'XColor','none','YColor','none','XLim',[0 size(waveData_mfilt,2)],'YLim',[0 size(waveData_mfilt,1)])
    title(ax(1),"raw data")
    hold(ax(1),"on")
    ax(2) = subplot(7,1,3:4);
    colormap(ax(2),C)
    clim(ax(2), [1 max(waveData_mfilt,[],"all")]);
    set(ax(2),'XColor','none','YColor','none','XLim',[0 size(waveData_mfilt,2)],'YLim',[0 size(waveData_mfilt,1)])
    title(ax(2),"wave data")
    hold(ax(2),"on")
    ax(3) = subplot(7,1,5:6);
    colormap(ax(3),C)
    clim(ax(3),[0 max(waveData_mfilt,[],"all")])
    set(ax(3),'XColor','none','YColor','none','XLim',[0 size(waveData_mfilt,2)],'YLim',[0 size(waveData_mfilt,1)])
    title(ax(3),"wave data + 3D median filter")
    hold(ax(3),"on")
    ax(4) = subplot(7,1,7);
    set(ax(4),'XLim',[0 timeVar])

    for n = 1:timeVar

        imagesc(ax(1),spikeRateMov(:,:,n))
        imagesc(ax(2),waveData(:,:,n))
        imagesc(ax(3),waveData_mfilt(:,:,n))
        plot(ax(4),sum(burstMatrix))
        hold on
        xline(n)
        hold off

        pause(0.05)
    end

end

%% Save wave data

waveDataName = extractBefore(fileName,'.h5');

choice2 = questdlg('Save wave data?','Save','Y','N','Y');

if strcmp(choice2,'N')
    return
else 

% save .mat variable
    cd(saveWaveDataPath)
    varFileName = waveDataName+"_waveData.mat";
    save(varFileName,"waveData","waveData_mfilt")
    
% export as a tiff
    cd(saveWaveTiffPath)
    tifFileName = waveDataName+"_waveData.tif";
    
    mov2export = waveData_mfilt;
    
    % dimensions
    [X, Y, T] = size(mov2export);
    
    % create a Tiff object for writing
    t = Tiff(tifFileName, 'w');
    
    % set TIFF tags
    tagstruct.ImageLength = X;
    tagstruct.ImageWidth = Y;
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.BitsPerSample = 16; % change as needed
    tagstruct.SamplesPerPixel = 1;
    tagstruct.RowsPerStrip = X;
    tagstruct.Compression = Tiff.Compression.None; % change as needed
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    
    % write each slice to the TIFF file
    for k = 1:T
        % set the current directory for the slice
        t.setTag(tagstruct);
    
        % write the slice data
        t.write(uint16(mov2export(:, :, k))); % convert to uint16 if necessary
    
        % write the next directory for the next slice
        if k < T
            t.writeDirectory();
        end
    end
    
    % close the Tiff object
    t.close();

end
