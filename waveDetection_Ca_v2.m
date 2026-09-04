%% Wave detection code for Calcium imaging
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
% v10_Ca - adapting waveDet10 to calcium imaging data + tiff export


clearvars, close all


%% User inputs

fileName =  'Mouse 2 retina 1_MMStack_Default.ome_8bit_binned4x_dfof_maskApplied.tif'; % file to analyze

%%% Recording parameters %%%
sampleRate = 100;  % sample rate in Hz

%%% Movie playback parameters %%%
playMovieToggle = 1; % set to 1 if you want to see a movie of waveID
movieSpeed = 0.001; % smaller # = faster speed

%%% Var output parameters %%%
saveWavaDataVar = 1; % toggle to save wave data variable (saves waveData and filtered waveData as a MATLAB variable)
exportAsTiff = 1; % toggle to save waveData as a .tif file (only saves filtered waveData; good for viewing movie in FIJI)

%%% Directories for finding and saving files %%%
% windows OS paths
% dataUploadPath = "C:\Users\yeagerkm\OneDrive - Vanderbilt\Tiriac Lab\Toolbox Project\6_calciumImaging\1_data";
% saveWaveDataPath = "C:\Users\yeagerkm\OneDrive - Vanderbilt\Tiriac Lab\Toolbox Project\6_calciumImaging\3_waveDetectionResults\1_waveData";
% saveWaveMoviePath = "C:\Users\yeagerkm\OneDrive - Vanderbilt\Tiriac Lab\Toolbox Project\6_calciumImaging\3_waveDetectionResults\2_waveMovies";

% mac OS paths
dataUploadPath = "/Users/kayleeodum/Library/CloudStorage/OneDrive-Vanderbilt/Tiriac Lab/Toolbox Project/6_calciumImaging/1_data";
saveWaveDataPath = "/Users/kayleeodum/Library/CloudStorage/OneDrive-Vanderbilt/Tiriac Lab/Toolbox Project/6_calciumImaging/3_waveDetectionResults/1_waveData_mat";
saveWaveMoviePath = "/Users/kayleeodum/Library/CloudStorage/OneDrive-Vanderbilt/Tiriac Lab/Toolbox Project/6_calciumImaging/3_waveDetectionResults/2_waveData_tif";

%% load the movie

cd(dataUploadPath)
V = tiffreadVolume(fileName);

[numRows, numCols, timeVar] = size(V);
spikeRateChannel = reshape(V,numRows*numCols,timeVar);

%% Identify baseline and ret wave periods

wholeRetSpikeRate = sum(spikeRateChannel);
[troughs, baselineInd] = findpeaks(wholeRetSpikeRate*(-1),"MinPeakDistance",10,"NPeaks",20,'SortStr','descend');

map = [1 1 1
    0 0 0];

fig1 = figure('Name','Data and baseline checks');
sgtitle('Data and Baseline Checks')
subplot(3,1,1)
plot(wholeRetSpikeRate,'k')
hold on
xline(baselineInd,'r')
subplot(3,1,2:3)
imagesc(spikeRateChannel)
xlim([0 timeVar])
colormap(map)
clim([0 0.2])
set(gca,'Box','on','XTick',[],'YTick',[],'YDir','reverse')

baselineFR = [];
for n = 1:length(baselineInd)
    baselineFR = [baselineFR, spikeRateChannel(:,baselineInd(n)-1:baselineInd(n)+1)];
end
meanBaselineFR = mean(baselineFR,2);

burstRateChannel = double(spikeRateChannel > meanBaselineFR+abs(meanBaselineFR)*5);
burstMov = reshape(burstRateChannel,numRows,numCols,timeVar);

%% show movies and test mean filtering

burstMov_filt = imboxfilt3(burstMov,[3 3 3]); % averages voxels w/i 3x3x3 cube
burstMov_filt = burstMov_filt > 0.2; % only keeps voxels with a mean > 0.X and sets them to 1 (everything else is 0)

if playMovieToggle == 1

    figure  
    colormap parula
    
    for n = 1:timeVar
        subplot(7,1,1:2)
        imagesc(V(:,:,n))
        clim([0 1]);
        set(gca,'XColor','none','YColor','none','XLim',[0 size(V,2)],'YLim',[0 size(V,1)])
        title("raw data")

        subplot(7,1,3:4)
        imagesc(burstMov(:,:,n))
        clim([0 1]);
        set(gca,'XColor','none','YColor','none','XLim',[0 size(V,2)],'YLim',[0 size(V,1)])
        title("burst data")

        subplot(7,1,5:6)
        imagesc(burstMov_filt(:,:,n))
        clim([0 1]);
        set(gca,'XColor','none','YColor','none','XLim',[0 size(V,2)],'YLim',[0 size(V,1)])
        title("burst data + mean filter")

        subplot(7,1,7)
        plot(sum(spikeRateChannel))
        hold on
        xline(n)
        hold off
        pause(movieSpeed)
    end
end

%% apply mean filtering to burstMov

meanFiltStep = 1;
if meanFiltStep == 1
    burstMov = double(burstMov_filt); % overwrite burstMov with burstMov_filt 
end

%% Wave detection

% Cube size
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

isPartOfWave = 3; % threshold for how many neighbors must have act to be part of a wave

% Pre-allocate wave data
waveData = nan(size(burstMov));
[numYs, numXs, recTime] = size(burstMov);

wCounterTot = [];

waveIDcounter = 1; % ID to be given to pixels that are part of the same wave (starts at 1)

% Wave detection
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
[uniqueOccurence, waveID] = groupcounts(reshape(waveData,[numYs*numXs*recTime 1]));
uniqueOccurence((isnan(waveID))) = [];
waveID((isnan(waveID))) = [];

%% Identify small events

waveArea = nan(length(waveID),1);
for n = 1:length(waveID)
    currWave = waveID(n);
    tempMat = waveData;
    tempMat(tempMat~=currWave) = nan;
    tempMat(tempMat==currWave) = 1;

    sumMat = sum(tempMat,3,'omitnan');
    areaMat = sumMat > 0;
    waveArea(n) = sum(areaMat,"all",'omitnan');
end

figure, histogram(waveArea,100)

%% Filter out small events

waveAreaThreshold = 0; % determine this number based on the natural breakpoint in the wave areas histogram distribution

for n = 1:length(waveID)
    currWave = waveID(n);
    if waveArea(n) < waveAreaThreshold 
        waveData(waveData == currWave) = nan;
    end
end

%% Apply a median filter to the data

waveData(isnan(waveData)) = 0;
waveData_mfilt = medfilt3(waveData,[3 3 3]);

%% play waveID movie

if playMovieToggle == 1

    C = colormap(rand(max(waveData_mfilt,[],"all"),3)*0.5+0.5);
    C = [0,0,0;C];
    figure("Units", "normalized", "Position",[0.33 0.2 0.33 0.7])
    
    ax(1) = subplot(7,1,1:2);
    colormap(ax(1),"parula")
    clim(ax(1),[0 0.05]);
    set(ax(1),'XColor','none','YColor','none','XLim',[0 size(V,2)],'YLim',[0 size(V,1)])
    title(ax(1),"raw data")
    hold(ax(1),"on")

    ax(2) = subplot(7,1,3:4);
    colormap(ax(2),C)
    clim(ax(2), [1 max(waveData_mfilt,[],"all")]);
    set(ax(2),'XColor','none','YColor','none','XLim',[0 size(V,2)],'YLim',[0 size(V,1)])
    title(ax(2),"wave data")
    hold(ax(2),"on")

    ax(3) = subplot(7,1,5:6);
    colormap(ax(3),C)
    clim(ax(3),[0 max(waveData_mfilt,[],"all")])
    set(ax(3),'XColor','none','YColor','none','XLim',[0 size(V,2)],'YLim',[0 size(V,1)])
    title(ax(3),"wave data + median filter")
    hold(ax(3),"on")

    ax(4) = subplot(7,1,7);
    
    for n = 1:timeVar

        imagesc(ax(1),flipud(V(:,:,n)))
        imagesc(ax(2),flipud(waveData(:,:,n)))
        imagesc(ax(3),flipud(waveData_mfilt(:,:,n)))
        plot(ax(4),sum(spikeRateChannel))
        hold on
        xline(ax(4),n)
        hold off

        pause(movieSpeed)
    end

end


%% save wave data variable (if toggle = 1)

if saveWavaDataVar == 1
    cd(saveWaveDataPath)
    varFileName = extractBefore(fileName,".tif")+"_waveData.mat";
    save(varFileName,"waveData","waveData_mfilt")
end

%% export as a tiff

if exportAsTiff == 1
    cd(saveWaveMoviePath)

    mov2export = waveData_mfilt; % only save filtered wave data

    % Dimensions
    [X, Y, T] = size(mov2export);
    
    % Specify the filename for the TIFF stack
    tifFileName = extractBefore(fileName,".tif")+"_waveData.tif";
    
    % Create a Tiff object for writing
    t = Tiff(tifFileName, 'w');
    
    % Set TIFF tags
    tagstruct.ImageLength = X;
    tagstruct.ImageWidth = Y;
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.BitsPerSample = 16; % change as needed
    tagstruct.SamplesPerPixel = 1;
    tagstruct.RowsPerStrip = X;
    tagstruct.Compression = Tiff.Compression.None; % change as needed
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    
    % Write each slice to the TIFF file
    for k = 1:T
        % Set the current directory for the slice
        t.setTag(tagstruct);
    
        % Write the slice data
        t.write(uint16(mov2export(:, :, k))); % convert to uint16 if necessary
    
        % Write the next directory for the next slice
        if k < T
            t.writeDirectory();
        end
    end
    
    % Close the Tiff object
    t.close();
end

