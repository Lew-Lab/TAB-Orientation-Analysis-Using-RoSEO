% 191206 - TD - Added full photon histogram (h5)

% 190618 - TD - Added numOfChunk and change of photonLowThreshold

% 190606 - v2 TD - Added mean and std calculation for ROI cross terms

% 190604 Tianben Ding
% Reset fibril ROI of analyzed data of estMMatrixRoSEO_v11.m or later by
% changing threROI

% This code is based on estMMatrixRoSEO_11.m

%% Analysis configuration
clear; close all;clc;

% Get output figure size (for large figures)
P0 = get(0,'ScreenSize');
P1(3:4) = P0(3:4)-150;

% Default figure setting
set(0, 'DefaultFigureColor','White',...
    'DefaultTextFontName','Arial','DefaultAxesFontName','Arial',...
    'DefaultAxesFontSize', 24,'DefaultTextFontSize',24,...
    'DefaultAxesLabelFontSizeMultiplier',1.2,...
    'DefaultAxesTickLength',[0.01 0.01],'DefaultAxesTickDir','out',...
    'DefaultAxesLineWidth',1,'DefaultLineLineWidth', 2);

% Add 'export_fig' and other useful tools into path of MATLAB
folderPath = pwd;
subFolders = {'functions';'utils';'phasemask'};
for f = 1:length(subFolders)
    addpath(fullfile(folderPath,subFolders{f}))
    
    subSubFolders = subdir(fullfile(folderPath,subFolders{f}));
    if ~isempty(subSubFolders)
        for ff = 1:length(subSubFolders)
            addpath(fullfile(folderPath,subFolders{f},subSubFolders{ff}))
        end
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parameters

% data path --
dataFolderPath = '210708_estMMatrixRoSEO_v16';
fileName = '190729 Data12 reg0.25 lineFitWavelet';

numOfChunk = 2;

dirNamePre = '';
saveNameNew = [fileName ',200 photonThre'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get information for saving
t=datetime('today');

dirNameNew=datestr(t,'yymmdd');
dirNameNew=[dirNamePre dirNameNew '_' mfilename];

if exist(dirNameNew,'dir') ~= 7
    mkdir(dirNameNew);
end

save([dirNameNew filesep saveNameNew ' initial setting' '.mat'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([dataFolderPath filesep fileName ' all analyzed results.mat'])

figure('Position',P1);
hold on
hist = histogram2(locPos(:,1),locPos(:,2),Xedges,Yedges,'DisplayStyle','tile','ShowEmptyBins','on','EdgeColor','none');
N1 = hist.Values;
N1 = reshape(N1,[size(N1,1)*size(N1,2),1]);
N1(N1 == 0) = [];
sNum = 0;
indN = 0;
while sNum < length(N1)*satPop
    indN = indN + 1;
    sNum = sNum + sum(N1 == indN);
end
tickMinBoth = 0;
tickMaxBoth = indN;
colormap(colorMap)
caxis([tickMinBoth,tickMaxBoth]);
axis image ij
ax = gca;
rectangle('Position',[ax.XLim(2)-scaleBar-(ax.XLim(2)-ax.XLim(1))/20,...
    ax.YLim(2)-(ax.YLim(2)-ax.YLim(1))/100-(ax.YLim(2)-ax.YLim(1))/20,...
    scaleBar,(ax.YLim(2)-ax.YLim(1))/100],'FaceColor','w','EdgeColor','w');
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'Box', 'off');
colorbar;
title('SR image of RoSEO','FontSize',24)
hold off
hROI0 = gcf;

reconstructed = hist.Values;
reconstructed = reconstructed.';

% ROI selection
binRec = reconstructed > threROI;
binRecFill = imfill(binRec,'holes');
CC = bwconncomp(binRecFill);
numPixels = cellfun(@numel,CC.PixelIdxList);
fibROITemp = zeros(size(reconstructed));
fibROI = zeros(size(reconstructed));
fibEdge = zeros(numOfChunk,2);
fibROIBound = cell(1,numOfChunk);
xROI = cell(1,numOfChunk);
yROI = cell(1,numOfChunk);
for n = 1:numOfChunk
    [~,idx] = max(numPixels);
    fibROITemp(CC.PixelIdxList{idx}) = 1;
    fibROI = (fibROI+fibROITemp) > 0;
    fibEdgeInd = CC.PixelIdxList{idx}(1)-1;
    fibROITemp = ~logical(fibROITemp);
    if mod(fibEdgeInd,size(fibROITemp,1)) == 0
        fibEdgeY = size(fibROITemp,1);
    else
        fibEdgeY = mod(fibEdgeInd,size(fibROITemp,1));
    end
    fibEdgeX = floor(fibEdgeInd/size(fibROITemp,1))+1;
    fibEdge(n,:) = [fibEdgeX fibEdgeY];
    
    fibROIBound{n} = bwtraceboundary(fibROITemp,[fibEdgeY fibEdgeX],'E',4,Inf,'clockwise');
    xROI{n} = fibROIBound{n}(:,2);
    yROI{n} = fibROIBound{n}(:,1);
    
    numPixels(idx) = 0;
    fibROITemp = zeros(size(reconstructed));
end
fibROIAlpha = fibROI.*0.3;

recSat = reconstructed;
recSat(reconstructed>indN) = indN;
recSat = round(recSat./indN.*length(colorMap));
I = ind2rgb(recSat,colorMap);
green = cat(3, zeros(size(recSat)),ones(size(recSat)), zeros(size(recSat)));
figure('Position',P1);
hold on
imagesc(I)
h = imagesc(green);
set(h,'AlphaData',fibROIAlpha);
for n = 1:numOfChunk
    plot(xROI{n},yROI{n},'w','LineWidth',1)
end
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'Box', 'off');
axis image ij
title({['Selected ROI']; ['transparent green: ROI, white line: boundary']},'FontSize',24)
hold off
hROI1 = gcf;

prompt = ['Original photonLowThreshold, threROI, numOfChunk \n'...
    '= ' num2str(photonLowThreshold) ' ' num2str(threROI) ' ' num2str(numOfChunk) '\n'...
    'Is this ROI good? \n'...
    'Yes->> Enter / No->> type "photonLowThreshold threROI numOfChunk" \n'];
threCon = input(prompt,'s');
while ~isempty(threCon)
    threTemp = str2num(threCon);
    close(hROI0,hROI1)
    
    photonLowThreshold = threTemp(1);
    threROI = threTemp(2);
    numOfChunk = threTemp(3);
    
    % apply new thresholding
    locPosTemp = locPos;
    ind = sig < photonLowThreshold;
    locPosTemp(ind,:) = [];
        
    figure('Position',P1);
    hold on
    hist = histogram2(locPosTemp(:,1),locPosTemp(:,2),Xedges,Yedges,'DisplayStyle','tile','ShowEmptyBins','on','EdgeColor','none');
    N1 = hist.Values;
    N1 = reshape(N1,[size(N1,1)*size(N1,2),1]);
    N1(N1 == 0) = [];
    sNum = 0;
    indN = 0;
    while sNum < length(N1)*satPop
        indN = indN + 1;
        sNum = sNum + sum(N1 == indN);
    end
    tickMinBoth = 0;
    tickMaxBoth = indN;
    colormap(colorMap)
    caxis([tickMinBoth,tickMaxBoth]);
    axis image ij
    ax = gca;
    rectangle('Position',[ax.XLim(2)-scaleBar-(ax.XLim(2)-ax.XLim(1))/20,...
        ax.YLim(2)-(ax.YLim(2)-ax.YLim(1))/100-(ax.YLim(2)-ax.YLim(1))/20,...
        scaleBar,(ax.YLim(2)-ax.YLim(1))/100],'FaceColor','w','EdgeColor','w');
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);
    set(gca, 'Box', 'off');
    colorbar;
    title('SR image of RoSEO','FontSize',24)
    hold off
    hROI0 = gcf;
        
    reconstructed = hist.Values;
    reconstructed = reconstructed.';
        
    % ROI selection
    binRec = reconstructed > threROI;
    binRecFill = imfill(binRec,'holes');
    CC = bwconncomp(binRecFill);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    fibROITemp = zeros(size(reconstructed));
    fibROI = zeros(size(reconstructed));
    fibEdge = zeros(numOfChunk,2);
    fibROIBound = cell(1,numOfChunk);
    xROI = cell(1,numOfChunk);
    yROI = cell(1,numOfChunk);
    for n = 1:numOfChunk
        [~,idx] = max(numPixels);
        fibROITemp(CC.PixelIdxList{idx}) = 1;
        fibROI = (fibROI+fibROITemp) > 0;
        fibEdgeInd = CC.PixelIdxList{idx}(1)-1;
        fibROITemp = ~logical(fibROITemp);
        if mod(fibEdgeInd,size(fibROITemp,1)) == 0
            fibEdgeY = size(fibROITemp,1);
        else
            fibEdgeY = mod(fibEdgeInd,size(fibROITemp,1));
        end
        fibEdgeX = floor(fibEdgeInd/size(fibROITemp,1))+1;
        fibEdge(n,:) = [fibEdgeX fibEdgeY];
        
        fibROIBound{n} = bwtraceboundary(fibROITemp,[fibEdgeY fibEdgeX],'E',4,Inf,'clockwise');
        xROI{n} = fibROIBound{n}(:,2);
        yROI{n} = fibROIBound{n}(:,1);
        
        numPixels(idx) = 0;
        fibROITemp = zeros(size(reconstructed));
    end
    fibROIAlpha = fibROI.*0.3;
    
    recSat = reconstructed;
    recSat(reconstructed>indN) = indN;
    recSat = round(recSat./indN.*length(colorMap));
    I = ind2rgb(recSat,colorMap);
    green = cat(3, zeros(size(recSat)),ones(size(recSat)), zeros(size(recSat)));
    figure('Position',P1);
    hold on
    imagesc(I)
    h = imagesc(green);
    set(h,'AlphaData',fibROIAlpha);
    for n = 1:numOfChunk
        plot(xROI{n},yROI{n},'w','LineWidth',1)
    end
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);
    set(gca, 'Box', 'off');
    axis image ij
    title({['Selected ROI']; ['transparent green: ROI, white line: boundary']},'FontSize',24)
    hold off
    hROI1 = gcf;
    
    prompt = ['New photonLowThreshold threROI numOfChunk: \n'...
        '= ' num2str(photonLowThreshold) ' ' num2str(threROI) ' ' num2str(numOfChunk) '\n'...
        'Are you sure? \n Yes->> Enter / No->> type "photonLowThreshold threROI numOfChunk" \n'];
    threCon = input(prompt,'s');
end

M(ind,:) = [];
sig(ind,:) = [];
frm(ind,:) = [];
locPos(ind,:) = [];
in = zeros(size(locPos,1),1);
for n = 1:numOfChunk
    inTemp = inpolygon(locPos(:,1),locPos(:,2),(xROI{n}-0.5).*binSize - (((sizeROI-1)/2+0.5)*pixelSize),(yROI{n}-0.5).*binSize- (((sizeROI-1)/2+0.5)*pixelSize));
    in = (in + inTemp) > 0;
end

figure('Position',P1);
hold on
for n = 1:numOfChunk
    plot((xROI{n}-0.5).*binSize-(((sizeROI-1)/2+0.5)*pixelSize),(yROI{n}-0.5).*binSize-(((sizeROI-1)/2+0.5)*pixelSize))
end
plot(locPos(in,1),locPos(in,2),'r+')
plot(locPos(~in,1),locPos(~in,2),'bo')
grid on
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'Box', 'on');
axis image ij
axis([Xedges(floor(length(Xedges)/2-length(Xedges)*shrinkFOVX/2)+1) Xedges(floor(length(Xedges)/2+length(Xedges)*shrinkFOVX/2))...
    Yedges(floor(length(Yedges)/2-length(Yedges)*shrinkFOVY/2)+1) Yedges(floor(length(Yedges)/2+length(Yedges)*shrinkFOVY/2))])
title({['Selected ROI vs localizations']; ['red cross: localizations inside, blue circle: localizations outside']},'FontSize',24)
hold off
hROI2 = gcf;

export_fig(hROI0,[dirNameNew filesep saveNameNew ' hROI0, reconstruction' saveFormat],'-transparent')
export_fig(hROI1,[dirNameNew filesep saveNameNew ' hROI1, 2Dhist, check ROI vs reconstruction, image' saveFormat],'-transparent')
export_fig(hROI2,[dirNameNew filesep saveNameNew ' hROI2, scatter, check ROI vs localizations, image' saveFormat],'-transparent')

close(hROI0,hROI1,hROI2)
clear hROI0 hROI1 hROI2

% save the reconstructed image as a tif file
[reconstructedIn,~,~] = histcounts2(locPos(in,1),locPos(in,2),Xedges,Yedges);
reconstructedIn = reconstructedIn.';
reconstructedIn(reconstructedIn>indN) = indN;
reconstructedIn = imgaussfilt(reconstructedIn,1);

reconstructedIn = reconstructedIn./max(max(reconstructedIn))*255;
reconstructedIn = uint8(reconstructedIn);
imwrite(reconstructedIn,[dirNameNew filesep saveNameNew ' reconstructedIn' '.tif'],'Compression','none')

figure('Position',P1);
for mInd = 1:6
    subplot(2,3,mInd)
    mHist = histogram(M(in,mInd));
    if mInd < 4
        mHist.BinEdges = 0:0.05:1;
    else
        meanC = mean(M(in,mInd));
        stdC = std(M(in,mInd));
        mHist.BinEdges = -0.5:0.05:0.5;
        leg = legend(['m=' num2str(round(meanC*1000)/1000) ',s=' num2str(round(stdC*1000)/1000)]);
        leg.FontSize = 10;
    end
    xlabel('Estimated value')
    title(['ROI' mLabel(mInd,:)])
end
h3_1 = gcf;
export_fig(h3_1,[dirNameNew filesep saveNameNew ' h3_1 distribution of estimated M in ROI'])
clear h3_1

figure('Position',P1);
histogram(sig(in));
xlabel('estimated emitted photons (maximum collectable photons if molecules oriented in-plane)','FontSize',18)
title(  [ 'ROI Photon number, median=' num2str( round(median(sig(in)),3,'significant') ) ]  );
h4 = gcf;
export_fig(h4,[dirNameNew filesep saveNameNew ' h4 histogram estimated emitted photons' saveFormat],'-transparent')
clear h4

figure('Position',P1);
histogram(sig);
xlabel('estimated emitted photons (maximum collectable photons if molecules oriented in-plane)','FontSize',18)
title(  [ 'Full photon number, median=' num2str( round(median(sig),3,'significant') ) ' mean=' num2str( round(mean(sig),3,'significant') ) ' num=' num2str( length(sig))]  );
h5 = gcf;
export_fig(h5,[dirNameNew filesep saveNameNew ' h5 histogram estimated emitted photons, full' saveFormat],'-transparent')
clear h5

clear I green
save([dirNameNew filesep saveNameNew ' all analyzed results' '.mat'])