% 200305 - v4 TD - Remove localizations outside of reasonable region after
% drift correction

% 200123 - v3 TD - Reset output of estMMatrixRoSEO_v11 or later.

% 190829 - v2 TD - Use super-resolved images instead of diffraction-limited
% ones. Also included correction of FOV rotation.

% 190826 Tianben Ding 
% Drift correlation by cross-correlation of wavelet filtered diffraction-
% limited images
% Ref of Wavelet filtering: M. Ovesny et al. ThunderSTORM: a comprehensive
% ImageJ plugin for PALM and STORM data analysis and super-resolution
% imaging - Methodology and Algorithms, ver 1.2

% This code is based on colocalizationAnalysis_WO3DefectWire_v5_1.m

clear; close all;clc; % clear everything

% Get output figure size (for large figures)
P0 = get(0,'ScreenSize');
P1(3:4) = P0(3:4)-150;

% Default figure setting
set(0, 'DefaultFigureColor','White',...
    'DefaultTextFontName','Arial','DefaultAxesFontName','Arial',...
    'DefaultAxesFontSize', 18,'DefaultTextFontSize',18,...
    'DefaultAxesLabelFontSizeMultiplier',1,...
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
fFormat = '.tif'; % file format
binSize = 20; % pixel size on objective plane (nm)

% Address and file name of reconstruction
fileAddressSRT1 = '210708_estMMatrixRoSEO_v16';% data folder containing the first time point
fileAddressSRT2 = '210708_estMMatrixRoSEO_v16';% data folder containing the later time point

fileNameSRT1 = '190729 Data12 reg0.25 lineFitWavelet';% file name of your first time point
fileNameSRT2 = '190729 Data12 reg0.25 lineFitWavelet';% file name of your later time point

fileAddressSRT1Val = '210709_estMMatrixRoSEO_resetROI_v2';% data folder containing your resetROI_v2 resutls on the first time point
fileNameSRT1Val = '190729 Data12 reg0.25 lineFitWavelet,200 photonThre';
% only needed for validation of the drift correction using detected ROI
% boundaries. Need to be an output of "extMMatricRoSEO_resetROI_v2" on the
% same data set as "fileNameSRT1".

pixelNum42DGaussFit = 3; %Pixel radius for 2D Gaussian fit, from center to edge, x and y directions

% Parameters for wavelet filtering with B-spline
scaleB = 2;
orderB = 3;

% save name
saveName = '190729 Data12vs12 test of code';
saveFormat = '.png';

% Input part end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get information for save
t=datetime('today');

dirName=datestr(t,'yymmdd');
dirName=[dirName '_' mfilename];

if exist(dirName,'dir') ~= 7
    mkdir(dirName);
end

save([dirName '\' saveName ' configuration' '.mat'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start analysis
% Read analyzed information
t = Tiff([fileAddressSRT1 filesep fileNameSRT1 ' reconstructed.tif'],'r');
SRImgT1 = t.read();
SRImgT1 = double(SRImgT1);
t.close();

t = Tiff([fileAddressSRT2 filesep fileNameSRT2 ' reconstructed.tif'],'r');
SRImgT2 = t.read();
SRImgT2 = double(SRImgT2);
t.close();

load([fileAddressSRT1 filesep fileNameSRT1 ' all analyzed results.mat'],'shrinkFOVX','shrinkFOVY')
% SRImgT1Temp = SRImgT1;
% SRImgT2Temp = SRImgT2;
sizeT1 = size(SRImgT1);
SRImgT1Temp = SRImgT1(round(sizeT1(1)/2-sizeT1(1)*shrinkFOVY/2):round(sizeT1(1)/2+sizeT1(1)*shrinkFOVY/2),...
    round(sizeT1(2)/2-sizeT1(2)*shrinkFOVX/2):round(sizeT1(2)/2+sizeT1(2)*shrinkFOVX/2));
sizeT2 = size(SRImgT2);
SRImgT2Temp = SRImgT2(round(sizeT2(1)/2-sizeT2(1)*shrinkFOVY/2):round(sizeT2(1)/2+sizeT2(1)*shrinkFOVY/2),...
    round(sizeT2(2)/2-sizeT2(2)*shrinkFOVX/2):round(sizeT2(2)/2+sizeT2(2)*shrinkFOVX/2));

figure('Position',P1,'Visible','on');
subplot(1,2,1)
imagesc(SRImgT1Temp)
axis image ij
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'Box', 'off');
title('T1')
subplot(1,2,2)
imagesc(SRImgT2Temp)
axis image ij
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'Box', 'off');
title('T2')
hFOV1 = gcf;

prompt = ['Original shrinkFOVX shrinkFOVY: \n'...
    '= ' num2str(shrinkFOVX) ' ' num2str(shrinkFOVY) '\n'...
    'Is this FOV good? \n'...
    'Yes->> Enter / No->> type "shrinkFOVX, shrinkFOVY [0-1]" \n'];
shrinkCon = input(prompt,'s');
while ~isempty(shrinkCon)
    shrinkTemp = str2num(shrinkCon);
    close(hFOV1)
    
    shrinkFOVX = shrinkTemp(1);
    shrinkFOVY = shrinkTemp(2);
    
    SRImgT1Temp = SRImgT1(round(sizeT1(1)/2-sizeT1(1)*shrinkFOVY/2):round(sizeT1(1)/2+sizeT1(1)*shrinkFOVY/2),...
        round(sizeT1(2)/2-sizeT1(2)*shrinkFOVX/2):round(sizeT1(2)/2+sizeT1(2)*shrinkFOVX/2));
    SRImgT2Temp = SRImgT2(round(sizeT2(1)/2-sizeT2(1)*shrinkFOVY/2):round(sizeT2(1)/2+sizeT2(1)*shrinkFOVY/2),...
        round(sizeT2(2)/2-sizeT2(2)*shrinkFOVX/2):round(sizeT2(2)/2+sizeT2(2)*shrinkFOVX/2));
    
    figure('Position',P1,'Visible','on');
    subplot(1,2,1)
    imagesc(SRImgT1Temp)
    axis image ij
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);
    set(gca, 'Box', 'off');
    title('T1')
    subplot(1,2,2)
    imagesc(SRImgT2Temp)
    axis image ij
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);
    set(gca, 'Box', 'off');
    title('T2')
    hFOV1 = gcf;
    
    prompt = ['New shrinkFOVX shrinkFOVY: \n'...
    '= ' num2str(shrinkFOVX) ' ' num2str(shrinkFOVY) '\n'...
    'Are you sure? \n Yes->> Enter / No->> type "shrinkFOVX, shrinkFOVY [0-1]" \n'];
    shrinkCon = input(prompt,'s');
end
export_fig(hFOV1,[dirName '\' saveName ' hFOV1, image, check of shrinked FOV' saveFormat],'-transparent')
close(hFOV1)
clear hFOV1

SRImgT1 = SRImgT1Temp;
SRImgT2 = SRImgT2Temp;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Check read images
% figure('Position',P1);
% subplot(1,2,1)
% imagesc(SRImgT1)
% axis image ij
% title('SR image of T1')
% subplot(1,2,2)
% imagesc(SRImgT2)
% axis image ij
% title('SR image of T2')
% set(gca,'Color','black')
% hold off;
% h0 = gcf;
% export_fig(h0,[dirName '\' saveName ' h0, image, SR images of T1 and T2' saveFormat],'-transparent')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Drift Correction and ROI selection
% Wavelet filtering with B-spline
% calculate kernels
kernelSize = 2*ceil(orderB*scaleB/2) - 1;
kernelInd = 1:kernelSize;
kernelInput = kernelInd-(kernelSize+1)/2;

kernel1 = B_spline_basis(kernelInput/scaleB + orderB/2,orderB);
normKernel = 1/sum(kernel1);
kernel1 = normKernel*kernel1;
kernel1 = transpose(kernel1);

kernel2 = reshape([kernel1.';zeros(size(kernel1.'))],[],1);
kernel2 = kernel2(1:end-1);

% wavelet level
sizeKernel1 = size(kernel1,1);
sizeKernel2 = size(kernel2,1);

convImg0 = SRImgT1;
convImg1 = conv2(kernel1,kernel1.',SRImgT1);
convImg1 = convImg1((sizeKernel1-1)/2+1:end-(sizeKernel1-1)/2,(sizeKernel1-1)/2+1:end-(sizeKernel1-1)/2);
convImg2 = conv2(kernel2,kernel2.',convImg1);
convImg2 = convImg2((sizeKernel2-1)/2+1:end-(sizeKernel2-1)/2,(sizeKernel2-1)/2+1:end-(sizeKernel2-1)/2);
waveletL1T1 = convImg0 - convImg1;
waveletL2T1 = convImg1 - convImg2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Position',P1,'Visible','off');
imagesc(SRImgT1);
axis image ij
colorbar
title('original SR T1 image')
hDriftC1 = gcf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Position',P1,'Visible','off');
imagesc(waveletL1T1);
axis image ij
colorbar
title('level 1 of wavelet transform: T1')
hDriftC2 = gcf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Position',P1,'Visible','off');
imagesc(waveletL2T1);
axis image ij
colorbar
title('level 2 of wavelet transform: T1')
hDriftC3 = gcf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

convImg0 = SRImgT2;
convImg1 = conv2(kernel1,kernel1.',SRImgT2);
convImg1 = convImg1((sizeKernel1-1)/2+1:end-(sizeKernel1-1)/2,(sizeKernel1-1)/2+1:end-(sizeKernel1-1)/2);
convImg2 = conv2(kernel2,kernel2.',convImg1);
convImg2 = convImg2((sizeKernel2-1)/2+1:end-(sizeKernel2-1)/2,(sizeKernel2-1)/2+1:end-(sizeKernel2-1)/2);
waveletL1T2 = convImg0 - convImg1;
waveletL2T2 = convImg1 - convImg2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Position',P1,'Visible','off');
imagesc(SRImgT2);
axis image ij
colorbar
title('original SR T2 image')
hDriftC4 = gcf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Position',P1,'Visible','off');
imagesc(waveletL1T2);
axis image ij
colorbar
title('level 1 of wavelet transform: T2')
hDriftC5 = gcf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Position',P1,'Visible','off');
imagesc(waveletL2T2);
axis image ij
colorbar
title('level 2 of wavelet transform: T2')
hDriftC6 = gcf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

export_fig(hDriftC1,[dirName '\' saveName ' hDriftC1, image, orifinal DL T1 image' saveFormat],'-transparent')
export_fig(hDriftC2,[dirName '\' saveName ' hDriftC2, image, L1 wavelet trans T1' saveFormat],'-transparent')
export_fig(hDriftC3,[dirName '\' saveName ' hDriftC3, image, L2 wavelet trans T1' saveFormat],'-transparent')
export_fig(hDriftC4,[dirName '\' saveName ' hDriftC4, image, orifinal DL T2 image' saveFormat],'-transparent')
export_fig(hDriftC5,[dirName '\' saveName ' hDriftC5, image, L1 wavelet trans T2' saveFormat],'-transparent')
export_fig(hDriftC6,[dirName '\' saveName ' hDriftC6, image, L2 wavelet trans T2' saveFormat],'-transparent')
close(hDriftC1,hDriftC2,hDriftC3,hDriftC4,hDriftC5,hDriftC6)
clear hDriftC1 hDriftC2 hDriftC3 hDriftC4 hDriftC5 hDriftC6

% calculate auto- and cross-correlation of the level 2 wavelet transforms
aCorrT1 = xcorr2(waveletL2T1);
% aCorrFA = xcorr2(dLImgFA);
cCorr = xcorr2(waveletL2T1,waveletL2T2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure output, cross-correlation
figure('Position',P1,'Visible','off');
imagesc(aCorrT1);
axis image ij
colorbar
title('Auto-correlation of T1')
hDriftC7 = gcf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure output, cross-correlation
figure('Position',P1,'Visible','off');
imagesc(cCorr);
axis image ij
colorbar
title('Cross-correlation of two SR images')
hDriftC8 = gcf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract the center part of each correlation image
[maxEachColumn,maxIndEachColumn] = max(aCorrT1,[],1);
[~,maxIndColumn] = max(maxEachColumn);
cenIndYT1 = maxIndEachColumn(maxIndColumn);
cenIndXT1 = maxIndColumn;
aCorrT1Cen = ...
    aCorrT1(cenIndYT1-pixelNum42DGaussFit:cenIndYT1+pixelNum42DGaussFit,...
    cenIndXT1-pixelNum42DGaussFit:cenIndXT1+pixelNum42DGaussFit);

[maxEachColumn,maxIndEachColumn] = max(cCorr,[],1);
[~,maxIndColumn] = max(maxEachColumn);
cenIndYCCorr = maxIndEachColumn(maxIndColumn);
cenIndXCCorr = maxIndColumn;
cCorrCen = ...
    cCorr(cenIndYCCorr-pixelNum42DGaussFit:cenIndYCCorr+pixelNum42DGaussFit,...
    cenIndXCCorr-pixelNum42DGaussFit:cenIndXCCorr+pixelNum42DGaussFit);

cenImgAxisX = ((-pixelNum42DGaussFit):pixelNum42DGaussFit);
cenImgAxisY = ((-pixelNum42DGaussFit):pixelNum42DGaussFit);
cenImgAxisX = cenImgAxisX.*binSize;
cenImgAxisY = cenImgAxisY.*binSize;

[x,y] = meshgrid(cenImgAxisX,cenImgAxisY);

% 2D Gaussian fit on auto-correlated image
[maxEachColumn,maxIndEachColumn] = max(aCorrT1Cen,[],1);
[maxValue,maxIndColumn] = max(maxEachColumn);

%Inital guess parameters [amplitude,mean in x dir,std in x dir,mean in y dir,std in y dir,clockwise rot]
x0 = [maxValue,x(maxIndEachColumn(maxIndColumn),maxIndColumn),(cenImgAxisX(end)-cenImgAxisX(1))/2,...
    y(maxIndEachColumn(maxIndColumn),maxIndColumn),(cenImgAxisY(end)-cenImgAxisY(1))/2,0];
% Lower bound
lb = [0,cenImgAxisX(1),0,cenImgAxisY(1),0,-pi/2];
% Upper bound
ub = [realmax('double'),cenImgAxisX(end),4*(cenImgAxisX(end)-cenImgAxisX(1)),cenImgAxisY(end),4*(cenImgAxisY(end)-cenImgAxisY(1)),pi/2];
% Least square curve fitting
coord = nan(size(x,1),size(y,2),2);
coord(:,:,1) = x;
coord(:,:,2) = y;
xOutACorr = lsqcurvefit(@twoDGaussFun_withRot,x0,coord,aCorrT1Cen,lb,ub);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure output, validation of the 2D Gaussian fit, auto-correlation
figure('Position',P1,'Visible','off');
surf(x,y,aCorrT1Cen,'FaceColor','none','EdgeColor','interp');
hold on
surf(x,y,twoDGaussFun_withRot(xOutACorr,coord),'EdgeColor','none','FaceAlpha',0.5)
hold off
xlabel('x [nm]')
ylabel('y [nm]')
zlabel('digital count')
title(['Auto-correlation of T1 SR Img and its 2D Gaussian fitting result, cen x=' num2str(xOutACorr(2)) ', cen y=' num2str(xOutACorr(4))])
hDriftC9 = gcf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2D Gaussian fit on cross-correlated image
[maxEachColumn,maxIndEachColumn] = max(cCorrCen,[],1);
[maxValue,maxIndColumn] = max(maxEachColumn);

%Inital guess parameters [amplitude,mean in x dir,std in x dir,mean in y dir,std in y dir,clockwise rot]
x0 = [maxValue,x(maxIndEachColumn(maxIndColumn),maxIndColumn),(cenImgAxisX(end)-cenImgAxisX(1))/2,...
    y(maxIndEachColumn(maxIndColumn),maxIndColumn),(cenImgAxisY(end)-cenImgAxisY(1))/2,0];
% Lower bound
lb = [0,cenImgAxisX(1),0,cenImgAxisY(1),0,-pi/2];
% Upper bound
ub = [realmax('double'),cenImgAxisX(end),4*(cenImgAxisX(end)-cenImgAxisX(1)),cenImgAxisY(end),4*(cenImgAxisY(end)-cenImgAxisY(1)),pi/2];
% Least square curve fitting
coord = nan(size(x,1),size(y,2),2);
coord(:,:,1) = x;
coord(:,:,2) = y;
xOutCCorr = lsqcurvefit(@twoDGaussFun_withRot,x0,coord,cCorrCen,lb,ub);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure output, validation of the 2D Gaussian fit, cross-correlation
figure('Position',P1,'Visible','off');
surf(x,y,cCorrCen,'FaceColor','none','EdgeColor','interp');
hold on
surf(x,y,twoDGaussFun_withRot(xOutCCorr,coord),'EdgeColor','none','FaceAlpha',0.5)
hold off
xlabel('x [nm]')
ylabel('y [nm]')
zlabel('digital count')
title(['Cross-correlation of SR Imgs and its 2D Gaussian fitting result, cen x=' num2str(xOutCCorr(2)) ', cen y=' num2str(xOutCCorr(4))])
hDriftC10 = gcf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate drift between two images
driftX = xOutCCorr(2) - xOutACorr(2) + (cenIndXCCorr-cenIndXT1)*binSize;
driftY = xOutCCorr(4) - xOutACorr(4) + (cenIndYCCorr-cenIndYT1)*binSize ;

% pixel level drift movement
driftXPix = round(driftX/binSize);
driftYPix = round(driftY/binSize);

export_fig(hDriftC7,[dirName '\' saveName ' hDriftC7, image, auto-correlation of T1 DL image' saveFormat],'-transparent')
export_fig(hDriftC8,[dirName '\' saveName ' hDriftC8, image, cross-correlation of 2 DL images' saveFormat],'-transparent')
export_fig(hDriftC9,[dirName '\' saveName ' hDriftC9, surf, 2D Gaussian fit on aCorr of T1' saveFormat],'-transparent')
export_fig(hDriftC10,[dirName '\' saveName ' hDriftC10, surf, 2D Gaussian fit on cCorr' saveFormat],'-transparent')

clear hDriftC7 hDriftC8 hDriftC9 hDriftC10

save([dirName '\' saveName ' output' '.mat'])

%% validation of correction
load([fileAddressSRT1Val filesep fileNameSRT1Val ' all analyzed results.mat'],'xROI','yROI','numOfChunk','colorMap')
load([fileAddressSRT2 filesep fileNameSRT2 ' all analyzed results.mat'],'locPos','Xedges','Yedges','indN')

hist = histogram2(locPos(:,1),locPos(:,2),Xedges,Yedges,'DisplayStyle','tile','ShowEmptyBins','on','EdgeColor','none');

reconstructed = hist.Values;
reconstructed = reconstructed.';

recSat = reconstructed;
recSat(reconstructed>indN) = indN;

locPosC(:,1) = locPos(:,1) + driftX;
locPosC(:,2) = locPos(:,2) + driftY;

histC = histogram2(locPosC(:,1),locPosC(:,2),Xedges,Yedges,'DisplayStyle','tile','ShowEmptyBins','on','EdgeColor','none');

reconstructedC = histC.Values;
reconstructedC = reconstructedC.';

recSatC = reconstructedC;
recSatC(reconstructedC>indN) = indN;

figure('Position',P1,'Visible','on');
subplot(1,2,1)
hold on
imagesc(recSat)
colormap(colorMap)
colorbar
for n = 1:numOfChunk
    plot(xROI{n},yROI{n},'w','LineWidth',1)
end
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'Box', 'off');
axis image ij
title({['SR image and ROI before correction']; ['white line: boundary']},'FontSize',24)
hold off

subplot(1,2,2)
hold on
imagesc(recSatC)
colormap(colorMap)
colorbar
for n = 1:numOfChunk
    plot(xROI{n},yROI{n},'w','LineWidth',1)
end
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'Box', 'off');
axis image ij
title({['Drift corrected ROI']; ['white line: boundary']},'FontSize',24)
hold off

hDriftC11 = gcf;

export_fig(hDriftC11,[dirName '\' saveName ' hDriftC11,with,without correction' saveFormat],'-transparent')

%% correct of analysis output
load([fileAddressSRT1 filesep fileNameSRT1 ' all analyzed results.mat'],'locPos')
locPos1 = locPos;
xLE1 = min(locPos1(:,1));
xRE1 = max(locPos1(:,1));
yBE1 = min(locPos1(:,2));
yTE1 = max(locPos1(:,2));

load([fileAddressSRT2 filesep fileNameSRT2 ' all analyzed results.mat'],'locPos')
locPos2 = locPos;
xLE2 = min(locPos2(:,1))+driftX;
xRE2 = max(locPos2(:,1))+driftX;
yBE2 = min(locPos2(:,2))+driftY;
yTE2 = max(locPos2(:,2))+driftY;

xLE = max(xLE1,xLE2);
xRE = min(xRE1,xRE2);
yBE = max(yBE1,yBE2);
yTE = min(yTE1,yTE2);

% re-save data 2
clearvars -except driftX driftY fileAddressSRT1 fileNameSRT1 fileAddressSRT2 fileNameSRT2 dirName xLE xRE yBE yTE
dirNameDrift = dirName;

% validation of correction
% load([fileAddressSRT1Reset filesep fileNameSRT1 ' all analyzed results.mat'],'xROI','yROI','numOfChunk','colorMap')
load([fileAddressSRT2 filesep fileNameSRT2 ' all analyzed results.mat'])

locPosOri = locPos;
locPos(:,1) = locPos(:,1) + driftX;
locPos(:,2) = locPos(:,2) + driftY;

xLEInd = locPos(:,1) < xLE;
xREInd = locPos(:,1) > xRE;
yBEInd = locPos(:,2) < yBE;
yTEInd = locPos(:,2) > yTE;
eInd = (xLEInd+xREInd+yBEInd+yTEInd) > 0;
locPos(eInd,:) = [];
frm(eInd) = [];
M(eInd,:) = [];
sig(eInd) = [];

figure('Position',P1,'Visible','off');
hold on
hist = histogram2(locPosOri(:,1),locPosOri(:,2),Xedges,Yedges,'DisplayStyle','tile','ShowEmptyBins','on','EdgeColor','none');
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
title('Data2 image of RoSEO original','FontSize',24)
hold off
h1 = gcf;
export_fig(h1,[dirNameDrift filesep saveName ' h1 data 2 reconstruction original' saveFormat],'-transparent')

figure('Position',P1,'Visible','off');
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
title('Data2 image of RoSEO driftCorr after loc removal','FontSize',24)
hold off
h2 = gcf;
export_fig(h2,[dirNameDrift filesep saveName ' h2 data 2 reconstruction driftCorr' saveFormat],'-transparent')

save([dirNameDrift filesep fileNameSRT2 ' all analyzed results.mat'])

% re-save data 1
clearvars -except fileAddressSRT1 fileNameSRT1 dirNameDrift xLE xRE yBE yTE

load([fileAddressSRT1 filesep fileNameSRT1 ' all analyzed results.mat'])

locPosOri = locPos;

xLEInd = locPos(:,1) < xLE;
xREInd = locPos(:,1) > xRE;
yBEInd = locPos(:,2) < yBE;
yTEInd = locPos(:,2) > yTE;
eInd = (xLEInd+xREInd+yBEInd+yTEInd) > 0;
locPos(eInd,:) = [];
frm(eInd) = [];
M(eInd,:) = [];
sig(eInd) = [];

figure('Position',P1,'Visible','off');
hold on
hist = histogram2(locPosOri(:,1),locPosOri(:,2),Xedges,Yedges,'DisplayStyle','tile','ShowEmptyBins','on','EdgeColor','none');
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
title('Data1 image of RoSEO original','FontSize',24)
hold off
h3 = gcf;
export_fig(h3,[dirNameDrift filesep saveName ' h3 data 1 reconstruction original' saveFormat],'-transparent')

figure('Position',P1,'Visible','off');
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
title('Data1 image of RoSEO driftCorr after loc removal','FontSize',24)
hold off
h4 = gcf;
export_fig(h4,[dirNameDrift filesep saveName ' h4 data 1 reconstruction driftCorr' saveFormat],'-transparent')

save([dirNameDrift filesep fileNameSRT1 ' all analyzed results.mat'])