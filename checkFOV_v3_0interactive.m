% 200319 - TD v3 - Read ThunderSTORM histogram instead of csv file for
% defining working FOV quickly.

% 200318 - TD v2 - Use reconstruction instead of scatter plots for defining
% FOV.

% 190415 Tianben Ding
% FOV checking loop for RoSE-O analysis
% This code is based on SMRegisterReconstLargeEnsRegine_v22.m.

function  [centerROI,sizeROI,shrinkFOVX,shrinkFOVY,fig] =...
    checkFOV_v3(centerROI,sizeROI,shrinkFOVX,shrinkFOVY,recHist,pixelSize,regPara,...
    tform,visuWinSize,colorMap,satPop,scaleBar)

boundaryXPosSM = regPara.boundaryXPosSM;
imageSizeHSM = regPara.imageSizeHSM;
imageSizeVSM = regPara.imageSizeVSM;

% Select zoomed FOV regions
% xPosRightTemp = num(num(:,2) > pixelSize*boundaryXPosSM,2);
% xPosRightTemp = xPosRightTemp - pixelSize*boundaryXPosSM;
% yPosRightTemp = num(num(:,2) > pixelSize*boundaryXPosSM,3);
% 
% xPosLeftTemp = num(num(:,2) <= pixelSize*boundaryXPosSM,2);
% xPosLeftTemp = -xPosLeftTemp + pixelSize*boundaryXPosSM;
% yPosLeftTemp = num(num(:,2) <= pixelSize*boundaryXPosSM,3);
recHistRight = recHist(:, ((boundaryXPosSM*size(recHist,2)/imageSizeHSM)+1) :end);
recHistLeft = fliplr(recHist(:,1: (boundaryXPosSM*size(recHist,2)/imageSizeHSM) ));

centerROIFlip = centerROI;
centerROIFlip(1) = -centerROI(1) + boundaryXPosSM;

zoomRegion = [centerROIFlip(1)-(sizeROI-1)/2 centerROIFlip(2)-(sizeROI-1)/2;...
    centerROIFlip(1)+(sizeROI-1)/2 centerROIFlip(2)+(sizeROI-1)/2];
zoomRegionRight = transformPointsInverse(tform,zoomRegion.*pixelSize);
zoomRegionRight = [floor(zoomRegionRight(1,1)/pixelSize)+1 floor(zoomRegionRight(2,1)/pixelSize)+1 floor(zoomRegionRight(1,2)/pixelSize)+1 floor(zoomRegionRight(2,2)/pixelSize)+1];
zoomRegion = [zoomRegion(1,1) zoomRegion(2,1) zoomRegion(1,2) zoomRegion(2,2)];

zoomRegionShrink = [centerROIFlip(1)-round(sizeROI*shrinkFOVX/2) centerROIFlip(2)-round(sizeROI*shrinkFOVY/2);...
    centerROIFlip(1)+round(sizeROI*shrinkFOVX/2) centerROIFlip(2)+round(sizeROI*shrinkFOVY/2)];
zoomRegionRightShrink = transformPointsInverse(tform,zoomRegionShrink.*pixelSize);
zoomRegionRightShrink = [floor(zoomRegionRightShrink(1,1)/pixelSize)+1 floor(zoomRegionRightShrink(2,1)/pixelSize)+1 floor(zoomRegionRightShrink(1,2)/pixelSize)+1 floor(zoomRegionRightShrink(2,2)/pixelSize)+1];
zoomRegionShrink = [zoomRegionShrink(1,1) zoomRegionShrink(2,1) zoomRegionShrink(1,2) zoomRegionShrink(2,2)];

Xedges = [ 1  boundaryXPosSM ];
Yedges = [ 1  imageSizeVSM ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scatter plot in both channel
figure('Position',visuWinSize);
subplot(2,2,[1 2])
hold on
imagesc(Xedges,Yedges,recHistLeft)
N1 = recHistLeft;
N1 = reshape(N1,[size(N1,1)*size(N1,2),1]);
N1(N1 == 0) = [];
sNum = 0;
indN = 0;
while sNum < length(N1)*satPop
    indN = indN + 1;
    sNum = sNum + sum(N1 == indN);
end
tickMinBothL = 0;
tickMaxBothL = indN;
colormap(colorMap)
caxis([tickMinBothL,tickMaxBothL]);
colorbar;
axis image ij
axis([(Xedges(end)-Xedges(1))/2-1.5*(Yedges(end)-Yedges(1)) (Xedges(end)-Xedges(1))/2+1.5*(Yedges(end)-Yedges(1)) ...
    Yedges(1) Yedges(end)])
ax = gca;
rectangle('Position',[ax.XLim(2)-4*scaleBar/pixelSize-(ax.XLim(2)-ax.XLim(1))/20,...
    ax.YLim(2)-(ax.YLim(2)-ax.YLim(1))/100-(ax.YLim(2)-ax.YLim(1))/20,...
    4*scaleBar/pixelSize,(ax.YLim(2)-ax.YLim(1))/100],'FaceColor','w','EdgeColor','w');
set(gca, 'XTick', [ceil(ax.XLim(1)/50)*50:50:floor(ax.XLim(2)/50)*50]);
set(gca, 'YTick', [ceil(ax.YLim(1)/50)*50:50:floor(ax.YLim(2)/50)*50]);
xtickpos = get(gca, 'xtick');
ytickpos = get(gca, 'ytick');
for row = 1:length(ytickpos)
    plot([ax.XLim(1) ax.XLim(2)], [ytickpos(row) ytickpos(row)], 'w:','LineWidth',1)
end
for col = 1:length(xtickpos)
    plot([xtickpos(col) xtickpos(col)], [ax.YLim(1) ax.YLim(2)], 'w:','LineWidth',1)
end
% set(gca, 'Box', 'off');
colorbar;
xlabel('x position [camera pix]','FontSize',10)
ylabel('y position [camera pix]','FontSize',10)
title('SR image of left channel','FontSize',10)
set(gca,'FontSize',10)
hold off

subplot(2,2,[3 4])
hold on
imagesc(Xedges,Yedges,recHistRight)
N1 = recHistRight;
N1 = reshape(N1,[size(N1,1)*size(N1,2),1]);
N1(N1 == 0) = [];
sNum = 0;
indN = 0;
while sNum < length(N1)*satPop
    indN = indN + 1;
    sNum = sNum + sum(N1 == indN);
end
tickMinBothR = 0;
tickMaxBothR = indN;
colormap(colorMap)
caxis([tickMinBothR,tickMaxBothR]);
colorbar;
axis image ij
axis([(Xedges(end)-Xedges(1))/2-1.5*(Yedges(end)-Yedges(1)) (Xedges(end)-Xedges(1))/2+1.5*(Yedges(end)-Yedges(1)) ...
    Yedges(1) Yedges(end)])
ax = gca;
rectangle('Position',[ax.XLim(2)-4*scaleBar/pixelSize-(ax.XLim(2)-ax.XLim(1))/20,...
    ax.YLim(2)-(ax.YLim(2)-ax.YLim(1))/100-(ax.YLim(2)-ax.YLim(1))/20,...
    4*scaleBar/pixelSize,(ax.YLim(2)-ax.YLim(1))/100],'FaceColor','w','EdgeColor','w');
set(gca, 'XTick', [ceil(ax.XLim(1)/50)*50:50:floor(ax.XLim(2)/50)*50]);
set(gca, 'YTick', [ceil(ax.YLim(1)/50)*50:50:floor(ax.YLim(2)/50)*50]);
xtickpos = get(gca, 'xtick');
ytickpos = get(gca, 'ytick');
for row = 1:length(ytickpos)
    plot([ax.XLim(1) ax.XLim(2)], [ytickpos(row) ytickpos(row)], 'w:','LineWidth',1)
end
for col = 1:length(xtickpos)
    plot([xtickpos(col) xtickpos(col)], [ax.YLim(1) ax.YLim(2)], 'w:','LineWidth',1)
end
% set(gca, 'Box', 'off');
colorbar;
xlabel('x position [camera pix]','FontSize',10)
ylabel('y position [camera pix]','FontSize',10)
title('SR image of right channel','FontSize',10)
set(gca,'FontSize',10)
hold off

fig.h0 = gcf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scatter plot in both channel (RoSE-O working region)
figure('Position',visuWinSize);
subplot(2,2,1)
hold on
imagesc(Xedges,Yedges,recHistLeft)
colormap(colorMap)
caxis([tickMinBothL,tickMaxBothL]);
colorbar;
axis image ij
axis(zoomRegion)
ax = gca;
rectangle('Position',[ax.XLim(2)-scaleBar/pixelSize-(ax.XLim(2)-ax.XLim(1))/20,...
    ax.YLim(2)-(ax.YLim(2)-ax.YLim(1))/100-(ax.YLim(2)-ax.YLim(1))/20,...
    scaleBar/pixelSize,(ax.YLim(2)-ax.YLim(1))/100],'FaceColor','w','EdgeColor','w');
% set(gca, 'Box', 'off');
colorbar;
xlabel('x position [camera pix]','FontSize',10)
ylabel('y position [camera pix]','FontSize',10)
title({['RoSE-O working region confirmation']; ['Localized left channel emitters']},'FontSize',10)
set(gca,'FontSize',10)
hold off

subplot(2,2,2)
hold on
imagesc(Xedges,Yedges,recHistRight)
colormap(colorMap)
caxis([tickMinBothR,tickMaxBothR]);
colorbar;
axis image ij
axis(zoomRegionRight)
ax = gca;
rectangle('Position',[ax.XLim(2)-scaleBar/pixelSize-(ax.XLim(2)-ax.XLim(1))/20,...
    ax.YLim(2)-(ax.YLim(2)-ax.YLim(1))/100-(ax.YLim(2)-ax.YLim(1))/20,...
    scaleBar/pixelSize,(ax.YLim(2)-ax.YLim(1))/100],'FaceColor','w','EdgeColor','w');
% set(gca, 'Box', 'off');
colorbar;
xlabel('x position [camera pix]','FontSize',10)
ylabel('y position [camera pix]','FontSize',10)
title({['RoSE-O working region confirmation']; ['Localized right channel emitters']},'FontSize',10)
set(gca,'FontSize',10)
hold off

subplot(2,2,3)
hold on
imagesc(Xedges,Yedges,recHistLeft)
colormap(colorMap)
caxis([tickMinBothL,tickMaxBothL]);
colorbar;
axis image ij
axis(zoomRegionShrink)
ax = gca;
rectangle('Position',[ax.XLim(2)-scaleBar/pixelSize-(ax.XLim(2)-ax.XLim(1))/20,...
    ax.YLim(2)-(ax.YLim(2)-ax.YLim(1))/100-(ax.YLim(2)-ax.YLim(1))/20,...
    scaleBar/pixelSize,(ax.YLim(2)-ax.YLim(1))/100],'FaceColor','w','EdgeColor','w');
% set(gca, 'Box', 'off');
colorbar;
xlabel('x position [camera pix]','FontSize',10)
ylabel('y position [camera pix]','FontSize',10)
title({['Zoom region confirmation']; ['Localized left channel emitters']},'FontSize',10)
set(gca,'FontSize',10)
hold off

subplot(2,2,4)
hold on
imagesc(Xedges,Yedges,recHistRight)
colormap(colorMap)
caxis([tickMinBothR,tickMaxBothR]);
colorbar;
axis image ij
axis(zoomRegionRightShrink)
ax = gca;
rectangle('Position',[ax.XLim(2)-scaleBar/pixelSize-(ax.XLim(2)-ax.XLim(1))/20,...
    ax.YLim(2)-(ax.YLim(2)-ax.YLim(1))/100-(ax.YLim(2)-ax.YLim(1))/20,...
    scaleBar/pixelSize,(ax.YLim(2)-ax.YLim(1))/100],'FaceColor','w','EdgeColor','w');
% set(gca, 'Box', 'off');
colorbar;
xlabel('x position [camera pix]','FontSize',10)
ylabel('y position [camera pix]','FontSize',10)
title({['Zoom region confirmation']; ['Localized right channel emitters']},'FontSize',10)
set(gca,'FontSize',10)
hold off

fig.h0_1 = gcf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0
prompt = ['Original centerROIx centerROIy sizeROI [pixel] shrinkFOVX shrinkFOVY [0-1]: \n'...
    '-> ' num2str(centerROI) ' ' num2str(sizeROI) ' ' num2str(shrinkFOVX) ' ' num2str(shrinkFOVY) '\n'...
    'Correlated zoomRegion (xStart xEnd yStart yEnd) [pixel]: \n'...
    num2str(zoomRegion) '\n'...
    'Is the zoom region reasonable? \n'...
    'Yes->> Enter / No->> zoomRegion (xStart xEnd yStart) [pixel] shrinkFOVX shrinkFOVY [0-1]: \n'];
zoomCon = input(prompt,'s');
while ~isempty(zoomCon)
    zoomTemp = str2num(zoomCon);
    close(fig.h0_1)
    
    centerROI = [-round((zoomTemp(2)+zoomTemp(1))/2)+boundaryXPosSM zoomTemp(3)+round((zoomTemp(2)-zoomTemp(1))/2)];
    sizeROI = zoomTemp(2)-zoomTemp(1)+1;
    shrinkFOVX = zoomTemp(4);
    shrinkFOVY = zoomTemp(5);
    
    centerROIFlip = centerROI;
    centerROIFlip(1) = -centerROI(1) + boundaryXPosSM;
    
    zoomRegion = [centerROIFlip(1)-(sizeROI-1)/2 centerROIFlip(2)-(sizeROI-1)/2;...
        centerROIFlip(1)+(sizeROI-1)/2 centerROIFlip(2)+(sizeROI-1)/2];
    zoomRegionRight = transformPointsInverse(tform,zoomRegion.*pixelSize);
    zoomRegionRight = [floor(zoomRegionRight(1,1)/pixelSize)+1 floor(zoomRegionRight(2,1)/pixelSize)+1 floor(zoomRegionRight(1,2)/pixelSize)+1 floor(zoomRegionRight(2,2)/pixelSize)+1];
    zoomRegion = [zoomRegion(1,1) zoomRegion(2,1) zoomRegion(1,2) zoomRegion(2,2)];
    
    zoomRegionShrink = [centerROIFlip(1)-round(sizeROI*shrinkFOVX/2) centerROIFlip(2)-round(sizeROI*shrinkFOVY/2);...
        centerROIFlip(1)+round(sizeROI*shrinkFOVX/2) centerROIFlip(2)+round(sizeROI*shrinkFOVY/2)];
    zoomRegionRightShrink = transformPointsInverse(tform,zoomRegionShrink.*pixelSize);
    zoomRegionRightShrink = [floor(zoomRegionRightShrink(1,1)/pixelSize)+1 floor(zoomRegionRightShrink(2,1)/pixelSize)+1 floor(zoomRegionRightShrink(1,2)/pixelSize)+1 floor(zoomRegionRightShrink(2,2)/pixelSize)+1];
    zoomRegionShrink = [zoomRegionShrink(1,1) zoomRegionShrink(2,1) zoomRegionShrink(1,2) zoomRegionShrink(2,2)];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Scatter plot in both channel (RoSE-O working region, double check)
    figure('Position',visuWinSize);
    subplot(2,2,1)
    hold on
    imagesc(Xedges,Yedges,recHistLeft)
    colormap(colorMap)
    caxis([tickMinBothL,tickMaxBothL]);
    colorbar;
    axis image ij
    axis(zoomRegion)
    ax = gca;
    rectangle('Position',[ax.XLim(2)-scaleBar/pixelSize-(ax.XLim(2)-ax.XLim(1))/20,...
        ax.YLim(2)-(ax.YLim(2)-ax.YLim(1))/100-(ax.YLim(2)-ax.YLim(1))/20,...
        scaleBar/pixelSize,(ax.YLim(2)-ax.YLim(1))/100],'FaceColor','w','EdgeColor','w');
    % set(gca, 'Box', 'off');
    colorbar;
    xlabel('x position [camera pix]','FontSize',10)
    ylabel('y position [camera pix]','FontSize',10)
    title({['RoSE-O working region confirmation']; ['Localized left channel emitters']},'FontSize',10)
    set(gca,'FontSize',10)
    hold off
    
    subplot(2,2,2)
    hold on
    imagesc(Xedges,Yedges,recHistRight)
    colormap(colorMap)
    caxis([tickMinBothR,tickMaxBothR]);
    colorbar;
    axis image ij
    axis(zoomRegionRight)
    ax = gca;
    rectangle('Position',[ax.XLim(2)-scaleBar/pixelSize-(ax.XLim(2)-ax.XLim(1))/20,...
        ax.YLim(2)-(ax.YLim(2)-ax.YLim(1))/100-(ax.YLim(2)-ax.YLim(1))/20,...
        scaleBar/pixelSize,(ax.YLim(2)-ax.YLim(1))/100],'FaceColor','w','EdgeColor','w');
    % set(gca, 'Box', 'off');
    colorbar;
    xlabel('x position [camera pix]','FontSize',10)
    ylabel('y position [camera pix]','FontSize',10)
    title({['RoSE-O working region confirmation']; ['Localized right channel emitters']},'FontSize',10)
    set(gca,'FontSize',10)
    hold off
    
    subplot(2,2,3)
    hold on
    imagesc(Xedges,Yedges,recHistLeft)
    colormap(colorMap)
    caxis([tickMinBothL,tickMaxBothL]);
    colorbar;
    axis image ij
    axis(zoomRegionShrink)
    ax = gca;
    rectangle('Position',[ax.XLim(2)-scaleBar/pixelSize-(ax.XLim(2)-ax.XLim(1))/20,...
        ax.YLim(2)-(ax.YLim(2)-ax.YLim(1))/100-(ax.YLim(2)-ax.YLim(1))/20,...
        scaleBar/pixelSize,(ax.YLim(2)-ax.YLim(1))/100],'FaceColor','w','EdgeColor','w');
    % set(gca, 'Box', 'off');
    colorbar;
    xlabel('x position [camera pix]','FontSize',10)
    ylabel('y position [camera pix]','FontSize',10)
    title({['Zoom region confirmation']; ['Localized left channel emitters']},'FontSize',10)
    set(gca,'FontSize',10)
    hold off
    
    subplot(2,2,4)
    hold on
    imagesc(Xedges,Yedges,recHistRight)
    colormap(colorMap)
    caxis([tickMinBothR,tickMaxBothR]);
    colorbar;
    axis image ij
    axis(zoomRegionRightShrink)
    ax = gca;
    rectangle('Position',[ax.XLim(2)-scaleBar/pixelSize-(ax.XLim(2)-ax.XLim(1))/20,...
        ax.YLim(2)-(ax.YLim(2)-ax.YLim(1))/100-(ax.YLim(2)-ax.YLim(1))/20,...
        scaleBar/pixelSize,(ax.YLim(2)-ax.YLim(1))/100],'FaceColor','w','EdgeColor','w');
    % set(gca, 'Box', 'off');
    colorbar;
    xlabel('x position [camera pix]','FontSize',10)
    ylabel('y position [camera pix]','FontSize',10)
    title({['Zoom region confirmation']; ['Localized right channel emitters']},'FontSize',10)
    set(gca,'FontSize',10)
    hold off
    fig.h0_1 = gcf;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    prompt = ['Correlated centerROI, sizeROI [pixel] and shrinkFOVX, shrinkFOVY [0-1]: \n'...
        '-> ' num2str(centerROI) ' ' num2str(sizeROI) ' ' num2str(shrinkFOVX) ' ' num2str(shrinkFOVY) '\n'...
        'Are you sure? \n Yes->> Enter / No->> zoomRegion (xStart xEnd yStart) [pixel] shrinkFOVX shrinkFOVY [0-1]: \n'];
    zoomCon = input(prompt,'s');
end
end
end