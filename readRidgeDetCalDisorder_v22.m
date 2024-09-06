% 210709 - TD - Use classical secondM2SymmConeWeighted_v7_1 for the second
% moment projection. Use assigned backbone orientations as the reference.

% 201213 - TD - Characterization of eigenvalues of second moment matrices.
% Use secondM2SymmConeWeighted_v11

% 201211 - v22 TD - secondM2SymmConeWeighted_v7 for projection. And
% adjusted for analyzing outputs of estMMatrixRoSEO_v14_1.m or later

% 200625 - v21 TD - Save analized position and orientation data into a
% separate .mat file

% 200521 - v20 TD - Changed quiver plot to line plot for reducing
% computational cost and removing plot shifts in molecule dipole plots
% (edge -> center)

% 200410 - v19 TD - Turn back to refractiveIndxSam for classical Nanoscope.
% Introduce parameters representing numApt.

% 191230 - v18 TD - Added groupPosLineX and groupPosLineY analysis for
% tagging grouped bursts with backbone indices

% 191204 - v17 TD - Commented out estimated theta correction

% 190914 - v16 TD - Save non-theta-filtered localizations

% 190914 - TD - Changed legend format of h0

% 190912 - v15 TD - also analyze off fibril localizations

% 190911 - v14_test TD - test secondM2SymmConeWeighted_v7_1 again

% 190911 - TD - Added psi2RefDegInplaneS analysis and storage

% 190906 - TD - Store photon number in BB analysis

% 190905 - TD - Matched phi and theta projection to
% validationRoSEOEstUsingSim_v9.m

% 190819 - v13 TD - Sliding window analysis along fibril backbone

% 190818 - TD - Added sorted scatter plots. Changed pixel bias of outputs
% of ridge detection (0.5->1 pixels).

% 190814 - v12 TD - Implemented a two step projector based on
% LSProjectorAnalysis_v2

% 190812 - v11 TD - Updated LS projector using
% secondM2SymmConeWeighted_v7_1 and applied theta thresholding. Moved
% nearNei calculation to later

% 190703 - TD - Flip signs of estimated mux and x axis direction of
% backbone orientation after all estimation (This is for visualization).

% 190701 - TD - Flip signs of x axis direction of backbone orientation
% instead of estimated mux signs.

% 190618 - v10 TD - Start estimation of mux, muy, muz and rotMobil from the
% nearest backbone orientations.

% 190618 - TD - Normalized estimated mux, muy and muz within ROI.

% 190617 - v9 TD - Use secondM2SymmConeWeighted_v3.m for angular space
% projection

% 190613 - v8 TD - Added blinking grouping algorithm for calculation
% tau_on. changed indices of pixelated visualization (rgbPlot~) from round
% to floor series. Added visualization of tau_on/photon per
% molecule/brightness[Hz]

% 190606 - v7 TD - Changed \psi _{ref} representation. Old (average
% orientation of neighbors) against (nearest backbone orientation) -> New
% (individual SM orientation) against (nearest backbone orientation). Use
% degree instead of cos(\psi).

% 190604 - v6 TD - Accept output of estMMatrixRoSEO_v11,
% estMMatrixRoSEO_restROI_v1, or their later version. Omitted
% ordered/disordered region comparison. Use the longest detected edge as
% the first reference angle in backbone orientation flipping.

% 190531 - v5 TD - use psi to represent dot product. Hemisphere rotation
% for solving in-plane orientation cancelation in average.

% 190523 - v4 TD - do not consider localizations whose neighbors are few in
% the heterogeneity maps

% 190515 - TD - moved back to the gamma estimation with its upper bound.

% 190514 - v3 TD - use mux, muy, muz to calculate dot product, remove bound
% in gamma estimation, solid angle

% 190508 - v2 TD - various orientation map including cos(phy-mean(phi)),
% cos(mean(phi)-ref), gamma, std(gamma), |gamma-mean(gamma)|, alpha,
% std(alpha), |alpha-mean(alpha)| using dot product of mux and muy (no
% muz) and mean calculation based on localizations nearby.

% 190417 Tianben Ding
% Read ridge detection results from ImageJ macro and calculate fibril
% disorder using estimated angle

%% Analysis configuration
clear; close all;clc;

% Get output figure size (for large figures)
P0 = get(0,'ScreenSize');
P1(3:4) = P0(3:4)-150;

% Default figure setting
set(0, 'DefaultFigureColor','White',...
    'DefaultTextFontName','Arial','DefaultAxesFontName','Arial',...
    'DefaultAxesFontSize', 20,'DefaultTextFontSize',20,...
    'DefaultTextColor','k',...
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
%% Input parameters - fixed parameters
colorMap = parula_hdr;
colorMapAng = parula_bcgyo;
colormapC = cold2hot(200);
colormapArrow = flipud(colormapC);

scaleBarMin = 50; % nm, scalebar size in zoomed-in SR reconstructions/scatter plots

arrowScalePix = 2; % arrow size in arrow visualization
arrowScale = 10;

% refractiveIndxSam=1.334; % refractive index of the sample, water: 1.334, matched: 1.518, air: 1
% sampleRefractiveIndx=1.334;
% numApt = 1.5;

theoSigma = 112;% nm, from backFocalPlaneSimEns_r5_1 for this wavelentgh and isotropic emitter
groupCoe = 3; % grouping coefficient, groupCoe*(localization precision) would be the grouping radii in the localization grouping part
expTime = 0.020; % s, exposure time on a detector

saveFormat = '.png';
%% variable parameters depending on data sets
lineFitNum = 4; % fitting number for obtaining ramp of fibrils

nearNeiDisNew = 30; % nm, radius of circles for calculating mean orientation of neighbors at each localization
lowNei = 0; % threshold out localizations with neighbors less than this value, only remove from the final visualization, all localization are preserved in analysis
lowSeg = 10;% threshold out segmented portions of fibrils with neighbors less than this value

plotLocInd = 1000; % localization index for demonstration region of circles defined by nearNeiDisNew

thetaThreDeg = 70; % deg, threshold out estimated theta smaller than this value, also larger than 180-thetaThreDeg [deg] (mirrored)

interPDist = 20; % nm, interpolation distance

visMargin = 100;

arrowScaleNoHead = 10;

visLineInd = [1 1]; % index of backbones for trajectory visualization at the end of the code, 2 integer at most, should be smaller than the number of backbone lines
% visRot = []; % input of camroll function
visZoom = 3; % input of camzoom function

maxDegPlot = 45; % degree, maximum degree in quiver color plots


% data path --
dataFolderPath = '210709_estMMatrixRoSEO_resetROI_v2';
fileName = '190729 Data12 reg0.25 lineFitWavelet,200 photonThre';
backgFile = ['G:\RawData\190729_2158BBHPTara laser irradiation tests\offsetSubtracted\Data12\BGEst_lineFitWavelet\4RoSEO' '\croppedData4RoSEO.h5'];

dirNamePre = '';
saveNameNewNew = [fileName ', nearNeiDis=' num2str(nearNeiDisNew) ', thThre=' num2str(thetaThreDeg)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get information for saving
t=datetime('today');

dirNameNewNew=datestr(t,'yymmdd');
dirNameNewNew=[dirNamePre dirNameNewNew '_' mfilename];

if exist(dirNameNewNew,'dir') ~= 7
    mkdir(dirNameNewNew);
end

save([dirNameNewNew filesep saveNameNewNew ' initial setting' '.mat'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% read detected fibril backbones and validation
num = readtable([dataFolderPath filesep fileName ' reconstructedIn ridgeDet.csv']);
num = table2array(num(:,[2 4:6]+1));
load([dataFolderPath filesep fileName ' all analyzed results.mat'])

numCRef = unique(num(:,1)); % detect unique ridge index
lineCoor = cell(1,length(numCRef)); % ridge coordinates in nm
lineCoorBin = cell(1,length(numCRef)); % ridge coordinates in bin
lineLength = nan(1,length(numCRef)); % obtain length information of ridge detection
for numInd = 1:length(numCRef)
    lineCoor{numInd} = (  (num(num(:,1)==numCRef(numInd),2:3)+1)  .*binSize - (((sizeROI-1)/2+0.5)*pixelSize) - 0.5*binSize).';
    lineCoorBin{numInd} = num(num(:,1)==numCRef(numInd),2:3).' + 1;
    %     lineCoor{numInd} = (  (num(num(:,1)==numCRef(numInd),2:3)+0.5)  .*binSize - (((sizeROI-1)/2+0.5)*pixelSize)).';
    %     lineCoorPix{numInd} = num(num(:,1)==numCRef(numInd),2:3).' + 0.5;
    %     lineCoorNum(numInd) = size(lineCoor{numInd},2);
    ind = num(:,1) == numCRef(numInd);
    ind = find(ind);
    lineLength(numInd) = num(ind(1),4);
end
% sort the ridgefit in terms of length
% [~,lineCoorInd] = sort(lineCoorNum,'descend');
[~,lineCoorInd] = sort(lineLength,'descend');
lineCoorTemp = lineCoor;
lineCoorPixTemp = lineCoorBin;
lineLengthTemp = lineLength;
for numInd = 1:length(numCRef)
    lineCoor{numInd} = lineCoorTemp{lineCoorInd(numInd)};
    lineCoorBin{numInd} = lineCoorPixTemp{lineCoorInd(numInd)};
    lineLength(numInd) = lineLengthTemp(lineCoorInd(numInd));
end

figure('Position',P1,'Visible','on');
hold on
imagesc(reconstructed)
tickMinBoth = 0;
tickMaxBoth = indN;
colormap(colorMap)
caxis([tickMinBoth,tickMaxBoth]);
axis image ij
ax = gca;
rectangle('Position',[ax.XLim(2)-scaleBar/binSize-(ax.XLim(2)-ax.XLim(1))/20,...
    ax.YLim(2)-(ax.YLim(2)-ax.YLim(1))/100-(ax.YLim(2)-ax.YLim(1))/20,...
    scaleBar/binSize,(ax.YLim(2)-ax.YLim(1))/100],'FaceColor','w','EdgeColor','w');
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'Box', 'off');
% colorbar;
title('SR image of RoSEO with estimated backbone','FontSize',24)
legendName = cell(1,length(numCRef));
for numInd = 1:length(numCRef)
    plot(lineCoorBin{numInd}(1,:).',lineCoorBin{numInd}(2,:).','o-');
    txt = ['\leftarrow ' num2str(numInd)];
    text(lineCoorBin{numInd}(1,1).',lineCoorBin{numInd}(2,1),txt,'Color','w','FontSize',20)
    %     legendName{numInd} = num2str(numInd);
end
% legend(legendName)
hold off
h0 = gcf;
export_fig(h0,[dirNameNewNew filesep saveNameNewNew ' h0, reconst, backbone of fibril' saveFormat],'-transparent')

%% calculate and validate fibril backbone slope
slopeLine = lineCoor; % slope of each lines (ridges) in Cartesian coordinate system
slopeInd = lineCoor; % number of lines contributed to calculate the final slopeCoor
lineCoorAll = []; % non-cell version of lineCoor
lineCoorPixAll = []; % non-cell version of lineCoor
slopeLineAll = []; % non-cell version of slopeLine
lineIndAll = []; % index of lines in all non-cell version arrays

interPointX = cell(1,length(numCRef)); % interpolated x coordinates
interPointY = cell(1,length(numCRef)); % interpolated y coordinates
interPointSlope = cell(1,length(numCRef)); % slopes of interpolated points

refRad = 0;

for lineInd = 1:length(lineCoor)
    if size(lineCoor{lineInd},2) < lineFitNum
        slopeLine{lineInd} = nan(1,size(lineCoor{lineInd},2));
        slopeInd{lineInd} = nan(1,size(lineCoor{lineInd},2));
        lineCoor{lineInd} = nan(2,size(lineCoor{lineInd},2));
        lineCoorBin{lineInd} = nan(2,size(lineCoor{lineInd},2));
    else
        slopeLine{lineInd} = zeros(1,size(lineCoor{lineInd},2));
        slopeInd{lineInd} = zeros(1,size(lineCoor{lineInd},2));
        for plotInd = 1:(length(lineCoor{lineInd})-lineFitNum+1)
            slopeTemp = polyfit(lineCoor{lineInd}( 1,plotInd:(plotInd+lineFitNum-1) ),lineCoor{lineInd}( 2,plotInd:(plotInd+lineFitNum-1) ),1);
            slopeTemp = atan(slopeTemp(1));
            if (slopeTemp - refRad) > pi/2
                slopeTemp = slopeTemp - pi;
            elseif (slopeTemp - refRad) < -pi/2
                slopeTemp = slopeTemp + pi;
            end
            refRad = slopeTemp;
            slopeTemp = repmat(slopeTemp,[1,lineFitNum]);
            slopeLine{lineInd}(plotInd:plotInd+lineFitNum-1) = slopeLine{lineInd}(plotInd:plotInd+lineFitNum-1) + slopeTemp;
            slopeInd{lineInd}(plotInd:plotInd+lineFitNum-1) = slopeInd{lineInd}(plotInd:plotInd+lineFitNum-1) + ones(1,lineFitNum);
        end
        slopeLine{lineInd}(:) = slopeLine{lineInd}(:)./slopeInd{lineInd}(:);
        
        % obtain interpolated positions based on interPDist
        distFromStart = cumsum([0, sqrt( (lineCoor{lineInd}(1,2:end)-lineCoor{lineInd}(1,1:end-1)).^2 +...
            (lineCoor{lineInd}(2,2:end)-lineCoor{lineInd}(2,1:end-1)).^2 )]);
        interPLocs = interPDist:interPDist:distFromStart(end);
        interPInd = interp1(distFromStart, 1:length(distFromStart), interPLocs);
        interPBaseLocs = floor(interPInd);
        interPWeight = interPInd - interPBaseLocs;
        
        interPMask = interPBaseLocs < length(distFromStart);
        finalX = lineCoor{lineInd}(1,end);
        finalY = lineCoor{lineInd}(2,end);
        finalSlope = slopeLine{lineInd}(end);
        finalInd = any(~interPMask) & (lineCoor{lineInd}(1,1)~=finalX|lineCoor{lineInd}(2,1)~=finalY);
        
        interPointX{lineInd} = [lineCoor{lineInd}(1,1),...
            lineCoor{lineInd}(1,interPBaseLocs(interPMask)).*(1-interPWeight(interPMask))+lineCoor{lineInd}(1,interPBaseLocs(interPMask)+1).*interPWeight(interPMask),...
            finalX(finalInd)];
        interPointY{lineInd} = [lineCoor{lineInd}(2,1),...
            lineCoor{lineInd}(2,interPBaseLocs(interPMask)).*(1-interPWeight(interPMask))+lineCoor{lineInd}(2,interPBaseLocs(interPMask)+1).*interPWeight(interPMask),...
            finalY(finalInd)];
        interPointSlope{lineInd} = [slopeLine{lineInd}(1),...
            slopeLine{lineInd}(interPBaseLocs(interPMask)).*(1-interPWeight(interPMask))+slopeLine{lineInd}(interPBaseLocs(interPMask)+1).*interPWeight(interPMask),...
            finalSlope(finalInd)];
    end
    lineCoorAll = [lineCoorAll lineCoor{lineInd}];
    lineCoorPixAll = [lineCoorPixAll lineCoorBin{lineInd}];
    slopeLineAll = [slopeLineAll slopeLine{lineInd}]; % radian
    lineIndAll = [lineIndAll repmat(lineInd,[1,length(slopeLine{lineInd})])];
end

lineCoorAll(:,isnan(lineCoorAll(1,:))) = [];
lineCoorPixAll(:,isnan(lineCoorPixAll(1,:))) = [];
lineIndAll(isnan(slopeLineAll)) = [];
slopeLineAll(isnan(slopeLineAll)) = [];


figure('Position',P1,'Visible','off');
hold on
imagesc(reconstructed)
tickMinBoth = 0;
tickMaxBoth = indN;
colormap(colorMap)
caxis([tickMinBoth,tickMaxBoth]);
axis image ij
ax = gca;
rec = rectangle('Position',[ax.XLim(2)-scaleBar/binSize-(ax.XLim(2)-ax.XLim(1))/20,...
    ax.YLim(2)-(ax.YLim(2)-ax.YLim(1))/100-(ax.YLim(2)-ax.YLim(1))/20,...
    scaleBar/binSize,(ax.YLim(2)-ax.YLim(1))/100],'FaceColor','w','EdgeColor','w');
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'Box', 'off');
colorbar;
title('SR image of RoSEO','FontSize',24)
quiver(lineCoorPixAll(1,:),lineCoorPixAll(2,:),arrowScalePix*cos(slopeLineAll),arrowScalePix*sin(slopeLineAll),...
    'AutoScale','off','LineWidth',3,'Color','r','ShowArrowHead','on');
hold off
h0_1 = gcf;
export_fig(h0_1,[dirNameNewNew filesep saveNameNewNew ' h0_1, reconst, backbone of fibril, with angle' saveFormat],'-transparent')

figure('Position',P1,'Visible','off');
hold on
for lineInd = 1:length(lineCoor)
    plot(interPointX{lineInd},interPointY{lineInd},'o','MarkerSize',5)
    quiver(interPointX{lineInd},interPointY{lineInd},arrowScale*cos(interPointSlope{lineInd}),arrowScale*sin(interPointSlope{lineInd}),...
        'AutoScale','off','LineWidth',3,'Color','r','ShowArrowHead','on');
    plot(lineCoor{lineInd}(1,:),lineCoor{lineInd}(2,:),'+','MarkerSize',10)
end
axis image ij
title('interpolated points vs detected ridge coordinates')
hold off
h0_1_1 = gcf;
export_fig(h0_1_1,[dirNameNewNew filesep saveNameNewNew ' h0_1_1, iterpolation check' saveFormat],'-transparent')

% calculate the nearest backbone direction of each localization
locPosROI = locPos(in,:);

locPosMeshX = repmat(locPosROI(:,1).',size(lineCoorAll,2),1);
locPosMeshY = repmat(locPosROI(:,2).',size(lineCoorAll,2),1);
lineCoorMeshX = repmat(lineCoorAll(1,:).',1,length(locPosMeshX));
lineCoorMeshY = repmat(lineCoorAll(2,:).',1,length(locPosMeshY));

allVectorDis = hypot(locPosMeshX-lineCoorMeshX,locPosMeshY-lineCoorMeshY);

clear locPosMeshX locPosMeshY lineCoorMeshX lineCoorMeshY

[~,minInd] = min(allVectorDis,[],1);

locPosROIRefPhi = slopeLineAll(minInd);
locPosROIRefMu = [cos(locPosROIRefPhi);sin(locPosROIRefPhi);zeros(size(locPosROIRefPhi))];
locPosROIRefMu = locPosROIRefMu.';
locPosROIRefMuR = sqrt(sum(locPosROIRefMu.^2,2));
locPosROIRefMuR = repmat(locPosROIRefMuR,[1,size(locPosROIRefMu,2)]);
locPosROIRefMu = locPosROIRefMu./locPosROIRefMuR;

locPosROIRefMu(:,1) = -locPosROIRefMu(:,1);

% categorize localization based on the nearest backbone orientaion
locPosROIRefInd = lineIndAll(minInd).';

% validation of localization index assignment with respect to backbones
figure('Position',P1,'Visible','off');
subplot(1,2,1)
hold on
scatter(locPosROI(locPosROIRefInd==visLineInd(1),1),locPosROI(locPosROIRefInd==visLineInd(1),2))
plot(lineCoor{visLineInd(1)}(1,:),lineCoor{visLineInd(1)}(2,:),'+-','MarkerSize',10)
legend('loc','backbone')
axis image ij
title(['validate fibril BB=' num2str(visLineInd(1)) ' and localization associated'])
hold off
subplot(1,2,2)
hold on
scatter(locPosROI(locPosROIRefInd==visLineInd(2),1),locPosROI(locPosROIRefInd==visLineInd(2),2))
plot(lineCoor{visLineInd(2)}(1,:),lineCoor{visLineInd(2)}(2,:),'+-','MarkerSize',10)
legend('loc','backbone')
axis image ij
title(['validate fibril BB=' num2str(visLineInd(2)) ' and localizations associated'])
hold off
h0_1_2 = gcf;
export_fig(h0_1_2,[dirNameNewNew filesep saveNameNewNew ' h0_1_2, validation of FB and localization association' saveFormat],'-transparent')

%% group localizations in adjacent frames
% read backg for estimation of physical angles
backg = h5read(backgFile,'/backg');

xPosAll = [];
yPosAll = [];
sigmaAll = [];
idAll = [];
photonAll = [];
backgroundAll = [];
precAll = [];

posOld = [];
photonOld = [];
precOld = [];
onTimeOld = [];

posNonCorrAll = [];
photonNonCorrAll = [];
precNonCorrAll = [];
idNonCorrAll= [];
onTimeNonCorrAll = [];

corrPos = [];

frmROI = frm(in);
sigROI = sig(in);

for h = 1:max(frmROI)
    id = transpose(frmROI(frmROI==h)); % frame index
    
    if ~isempty(id) % if the frame has localizations, run the following process
        xPos = transpose(locPosROI(frmROI==h,1));
        yPos = transpose(locPosROI(frmROI==h,2));
        sigma = transpose(repmat(theoSigma,[length(xPos),1]));
        photon = transpose(sigROI(frmROI==h));
        
        backgXPix = sum(sum(backg(:,1:sizeROI,h)))/sizeROI/sizeROI;
        backgYPix = sum(sum(backg(:,(sizeROI+1):end,h)))/sizeROI/sizeROI;
        backgPix = backgXPix + backgYPix;
        background = transpose(repmat(backgPix,[length(xPos),1]));
        
        % Calculate localization precision based on captured photon number
        % (in least square sense)
        tau = 2*pi*background.*(sigma.^2+pixelSize^2/12)./photon./(pixelSize^2);
        prec = (sigma.^2+pixelSize^2/12)./photon.*(16/9+4*tau);
        prec = sqrt(prec);
        
        % Store localization information
        xPosAll = [xPosAll xPos];
        yPosAll = [yPosAll yPos];
        sigmaAll = [sigmaAll sigma];
        idAll = [idAll id];
        photonAll = [photonAll photon];
        backgroundAll = [backgroundAll background];
        precAll = [precAll prec];
    end
    % Group localizations in adjacent frames
    posNew = [xPosAll(idAll==h);yPosAll(idAll==h)];
    photonNew = photonAll(idAll==h);
    precNew = precAll(idAll==h);
    [posOld,photonOld,precOld,onTimeOld,corrPos,posNonCorrAll,photonNonCorrAll,precNonCorrAll,idNonCorrAll,onTimeNonCorrAll] = ...
        groupLocAdjFrame_v1(posNew,photonNew,precNew,posOld,photonOld,precOld,onTimeOld,corrPos,...
        posNonCorrAll,photonNonCorrAll,precNonCorrAll,idNonCorrAll,onTimeNonCorrAll,groupCoe,h);
end
% store information of localizations in the last frame after the grouping
posCorrAll = [posNonCorrAll posOld];
posCorrAll = posCorrAll.';
photonCorrAll = [photonNonCorrAll photonOld];
photonCorrAll = photonCorrAll.';
precCorrAll = [precNonCorrAll precOld];
precCorrAll = precCorrAll.';
idCorrAll = [idNonCorrAll h*ones(1,size(posOld,2))];
idCorrAll = idCorrAll.';
onTimeCorrAll = [onTimeNonCorrAll onTimeOld];
onTimeCorrAll = onTimeCorrAll.*expTime;
onTimeCorrAll = onTimeCorrAll.';

brightnessHzAll = photonCorrAll./onTimeCorrAll./1000; %kHz

% store grouped localizations into different valuables
locPosGroup = posCorrAll;

figure('Position',P1,'Visible','off');
hold on
histG = histogram2(locPosGroup(:,1),locPosGroup(:,2),Xedges,Yedges,'DisplayStyle','tile','ShowEmptyBins','on','EdgeColor','none');
histValuesG = histG.Values;
N1G = histG.Values;
N1G = reshape(N1G,[size(N1G,1)*size(N1G,2),1]);
N1G(N1G == 0) = [];
sNumG = 0;
indNG = 0;
while sNumG < length(N1G)*satPop
    indNG = indNG + 1;
    sNumG = sNumG + sum(N1G == indNG);
end
tickMinBoth = 0;
tickMaxBoth = indNG;
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
title('SR image of RoSEO (grouped)','FontSize',24)
hold off
h0_4 = gcf;
export_fig(h0_4,[dirNameNewNew filesep saveNameNewNew ' h0_4 reconstruction with grouped localizations' saveFormat],'-transparent')
clear h0_4

% calculate the nearest backbone direction of each grouped localization
locPosMeshX = repmat(locPosGroup(:,1).',size(lineCoorAll,2),1);
locPosMeshY = repmat(locPosGroup(:,2).',size(lineCoorAll,2),1);
lineCoorMeshX = repmat(lineCoorAll(1,:).',1,length(locPosMeshX));
lineCoorMeshY = repmat(lineCoorAll(2,:).',1,length(locPosMeshY));

allVectorDis = hypot(locPosMeshX-lineCoorMeshX,locPosMeshY-lineCoorMeshY);

clear locPosMeshX locPosMeshY lineCoorMeshX lineCoorMeshY

[~,minIndGroup] = min(allVectorDis,[],1);

% categorize grouped localization based on the nearest backbone orientaion
locPosGroupRefInd = lineIndAll(minIndGroup).';

%% estimation of physical angle from estimated M matrix
% projection to physical angles
muxROI = nan(sum(in),1);
muyROI = nan(sum(in),1);
rotMobilROI = nan(sum(in),1);
muzROI = nan(sum(in),1);
% thetaROI = nan(sum(in),1);
% alphaROI = nan(sum(in),1);
% omegaROI = nan(sum(in),1);
% phiROI = nan(sum(in),1);
MROI = M(in,:);
eigVROI = nan(3,sum(in));

% redefine bases and normalization factors for the projection
n1=Nanoscope('pixelSize',pixelSize,'imageSize',sizeROI,'ADcount',1,'emissWavelength',emitWaveL,...
    'sampleRefractiveIndx',sampleRefractiveIndx,'phasemaskpara',phasemaskpara,'numApt',numApt,...
    'ZernikeCoefficient',zernikecoefficient);

Emitter_t.position_para.x = 0.0; % dummy
Emitter_t.position_para.y = 0.0;
Emitter_t.position_para.z = 0.0;

[bx.XX,bx.YY,bx.ZZ, bx.XY,bx.XZ,bx.YZ,...
    by.XX,by.YY,by.ZZ,by.XY,by.XZ,by.YZ]=...
    Nanoscope.computeBases(n1,Emitter_t);

%handle for croping region of interest
up_sample=n1.pixelUpsample;
img_size=n1.imageSize;
N_pupil=size(n1.phaseMask,1);

roi=@(img)img(-up_sample*(img_size-1)/2+N_pupil/2+2:1:up_sample*(img_size-1)/2+N_pupil/2+2,....
    -up_sample*(img_size-1)/2+N_pupil/2+2:1:up_sample*(img_size-1)/2+N_pupil/2+2,:);


bx.XX = roi(bx.XX);
bx.YY = roi(bx.YY);
bx.ZZ = roi(bx.ZZ);
bx.XY = roi(bx.XY);
bx.XZ = roi(bx.XZ);
bx.YZ = roi(bx.YZ);

by.XX = roi(by.XX);
by.YY = roi(by.YY);
by.ZZ = roi(by.ZZ);
by.XY = roi(by.XY);
by.XZ = roi(by.XZ);
by.YZ = roi(by.YZ);

% normalization factor
Emitter_t.polar_para.phiD=pi/2;
Emitter_t.polar_para.thetaD=pi/2;
Emitter_t.position_para.x=0;
Emitter_t.position_para.y=0;
Emitter_t.position_para.z=0;

[brightness_scalingX,brightness_scalingY]=n1.simDipole_novotny(n1,Emitter_t);
brightness_scaling = brightness_scalingX + brightness_scalingY;
sumNorm = sum(sum(roi(brightness_scaling)));

% calculate Hadamard products of basis for FIM calculation in LS projection
XXx = n1.XXxBasis;
YYx = n1.YYxBasis;
ZZx = n1.ZZxBasis;
XYx = n1.XYxBasis;
XZx = n1.XZxBasis;
YZx = n1.YZxBasis;

XXy = n1.XXyBasis;
YYy = n1.YYyBasis;
ZZy = n1.ZZyBasis;
XYy = n1.XYyBasis;
XZy = n1.XZyBasis;
YZy = n1.YZyBasis;

% Hadamard products, for x channel
Bx.aa = (XXx).*(XXx);
Bx.ab = (XXx).*(YYx);
Bx.ac = (XXx).*(ZZx);
Bx.ad = (XXx).*(XYx);
Bx.ae = (XXx).*(XZx);
Bx.af = (XXx).*(YZx);

Bx.bb = (YYx).*(YYx);
Bx.bc = (YYx).*(ZZx);
Bx.bd = (YYx).*(XYx);
Bx.be = (YYx).*(XZx);
Bx.bf = (YYx).*(YZx);

Bx.cc = (ZZx).*(ZZx);
Bx.cd = (ZZx).*(XYx);
Bx.ce = (ZZx).*(XZx);
Bx.cf = (ZZx).*(YZx);

Bx.dd = (XYx).*(XYx);
Bx.de = (XYx).*(XZx);
Bx.df = (XYx).*(YZx);

Bx.ee = (XZx).*(XZx);
Bx.ef = (XZx).*(YZx);

Bx.ff = (YZx).*(YZx);

% for y channel
By.aa = (XXy).*(XXy);
By.ab = (XXy).*(YYy);
By.ac = (XXy).*(ZZy);
By.ad = (XXy).*(XYy);
By.ae = (XXy).*(XZy);
By.af = (XXy).*(YZy);

By.bb = (YYy).*(YYy);
By.bc = (YYy).*(ZZy);
By.bd = (YYy).*(XYy);
By.be = (YYy).*(XZy);
By.bf = (YYy).*(YZy);

By.cc = (ZZy).*(ZZy);
By.cd = (ZZy).*(XYy);
By.ce = (ZZy).*(XZy);
By.cf = (ZZy).*(YZy);

By.dd = (XYy).*(XYy);
By.de = (XYy).*(XZy);
By.df = (XYy).*(YZy);

By.ee = (XZy).*(XZy);
By.ef = (XZy).*(YZy);

By.ff = (YZy).*(YZy);

for l = 1:subLoopNum
    backgTemp = backg(:,:, ((l-1)*anaLoop+1) : (l*anaLoop) );
    firstAnaLoopMInd = zeros(size(frmROI));
    lastAnaLoopMInd = zeros(size(frmROI));
    firstShift = 0;
    lastShift = 0;
    while sum(firstAnaLoopMInd) == 0
        firstAnaLoopMInd = frmROI == ( (l-1)*anaLoop+1 + firstShift );
        firstAnaLoopMInd = find(firstAnaLoopMInd);
        firstShift = firstShift + 1;
    end
    firstAnaLoopMInd = firstAnaLoopMInd(1);
    while sum(lastAnaLoopMInd) == 0
        lastAnaLoopMInd = frmROI == ( l*anaLoop - lastShift );
        lastAnaLoopMInd = find(lastAnaLoopMInd);
        lastShift = lastShift + 1;
    end
    lastAnaLoopMInd = lastAnaLoopMInd(end);
    
    parfor ii = firstAnaLoopMInd:lastAnaLoopMInd
        %         [muxROI(ii),muyROI(ii),muzROI(ii),~]= ...
        %             secondM2SymmConeWeightedMuOnly_v1_1(bx,by,Bx,By,sumNorm,MROI(ii,:),sigROI(ii),backgTemp(:,:,frmROI(ii)-(l-1)*anaLoop),locPosROIRefMu(ii,:));
        %
        %         [rotMobilROI(ii),~]= ...
        %             secondM2SymmConeWeightedGammaOnly_v1_1(bx,by,Bx,By,sumNorm,MROI(ii,:),sigROI(ii),backgTemp(:,:,frmROI(ii)-(l-1)*anaLoop),locPosROIRefMu(ii,:),...
        %             muxROI(ii),muyROI(ii),muzROI(ii));
        %         [muxROI(ii),muyROI(ii),muzROI(ii),~]= ...
        %             secondM2SymmConeWeightedMuOnly_v1_2(bx,by,Bx,By,sumNorm,MROI(ii,:),sigROI(ii),backgTemp(:,:,frmROI(ii)-(l-1)*anaLoop),locPosROIRefMu(ii,:));
        %
        %         [rotMobilROI(ii),~]= ...
        %             secondM2SymmConeWeightedGammaOnly_v1_1(bx,by,Bx,By,sumNorm,MROI(ii,:),sigROI(ii),backgTemp(:,:,frmROI(ii)-(l-1)*anaLoop),locPosROIRefMu(ii,:),...
        %             muxROI(ii),muyROI(ii),muzROI(ii));
        %                 [muxROI(ii),muyROI(ii),muzROI(ii),rotMobilROI(ii),~]= ...
        %                     secondM2SymmConeWeighted_v4(bx,by,Bx,By,MROI(ii,:),sigROI(ii),backgTemp(:,:,frmROI(ii)-(l-1)*anaLoop),locPosROIRefMu(ii,:));
        [muxROI(ii),muyROI(ii),muzROI(ii),rotMobilROI(ii),~]= ...
            secondM2SymmConeWeighted_v7_1(bx,by,Bx,By,sumNorm,MROI(ii,:),sigROI(ii),backgTemp(:,:,frmROI(ii)-(l-1)*anaLoop),locPosROIRefMu(ii,:));
    end
end

clear backgTemp

muxROI = -muxROI;
muR = sqrt(muxROI.^2+muyROI.^2+muzROI.^2);
muxROI = muxROI./muR;
muyROI = muyROI./muR;
muzROI = muzROI./muR;

% muxROIAng = muxROI;
% muyROIAng = muyROI;
% muzROIAng = muzROI;

% flip the estimated each molecule orientation based on its deviation from
% the nearest backbone direction
locPosROIRefMu(:,1) = - locPosROIRefMu(:,1);
dotMu2Ref = muxROI.*locPosROIRefMu(:,1) + muyROI.*locPosROIRefMu(:,2) + muzROI.*locPosROIRefMu(:,3);
muxROI(dotMu2Ref < 0) = -muxROI(dotMu2Ref < 0);
muyROI(dotMu2Ref < 0) = -muyROI(dotMu2Ref < 0);
muzROI(dotMu2Ref < 0) = -muzROI(dotMu2Ref < 0);

thetaROI = acos(muzROI);
% thetaROI(abs(pi/2 - (pi-thetaROI)) < abs(pi/2 - thetaROI)) = pi - thetaROI(abs(pi/2 - (pi-thetaROI)) < abs(pi/2 - thetaROI));
% thetaROI((thetaROI - pi/2) > pi/2) = thetaROI((thetaROI - pi/2) > pi/2) - pi;
% muxROIAng((thetaROI - pi/2) > pi/2) = - muxROIAng((thetaROI - pi/2) > pi/2);
% muyROIAng((thetaROI - pi/2) > pi/2) = - muyROIAng((thetaROI - pi/2) > pi/2);

% thetaROI((thetaROI - pi/2) < -pi/2) = thetaROI((thetaROI - pi/2) < -pi/2) + pi; % do not have this situation due to the range of emitter.theta
% muxROIAng((thetaROI - pi/2) < -pi/2) = - muxROIAng((thetaROI - pi/2) < -pi/2);
% muyROIAng((thetaROI - pi/2) < -pi/2) = - muyROIAng((thetaROI - pi/2) < -pi/2);

% locPosROIRefPhi = atan2(locPosROIRefMu(:,2),locPosROIRefMu(:,1));
locPosROIRefPhi = locPosROIRefPhi.';
phiROI = atan2(muyROI,muxROI);%atan2(muyROIAng,muxROIAng);
phiROI((locPosROIRefPhi - phiROI) > pi) = phiROI((locPosROIRefPhi - phiROI) > pi) + 2*pi;
phiROI((locPosROIRefPhi - phiROI) < -pi) = phiROI((locPosROIRefPhi - phiROI) < -pi) - 2*pi;

alphaROI = acos( (  -1 + sqrt(1+8*rotMobilROI)  )/2);
omegaROI = 2.*4.*pi.*sin(alphaROI/2).^2;

validThetaID1 = rad2deg(thetaROI) >= thetaThreDeg;
validThetaID2 = rad2deg(thetaROI) <= (180-thetaThreDeg);
validThetaID = (validThetaID1 + validThetaID2) >1;

figure('Position',P1,'Visible','off');
subplot(2,2,1)
histogram(rad2deg(alphaROI));
xlabel('[deg]','FontSize',18)
title('ROI \alpha')
subplot(2,2,2)
histogram(rotMobilROI);
title('ROI \gamma')
subplot(2,2,3)
histogram(rad2deg(thetaROI));
xlabel('[deg]','FontSize',18)
title('ROI \theta')
subplot(2,2,4)
histogram(rad2deg(phiROI));
xlabel('[deg]','FontSize',18)
title('ROI \phi')
h1 = gcf;

figure('Position',P1,'Visible','off');
subplot(2,2,1)
histogram(rad2deg(alphaROI(validThetaID)));
xlabel('[deg]','FontSize',18)
title('ROI \alpha')
subplot(2,2,2)
histogram(rotMobilROI(validThetaID));
title('ROI \gamma')
subplot(2,2,3)
histogram(rad2deg(thetaROI(validThetaID)));
xlabel('[deg]','FontSize',18)
title('ROI \theta')
subplot(2,2,4)
histogram(rad2deg(phiROI(validThetaID)));
xlabel('[deg]','FontSize',18)
title('ROI \phi')
h1_1 = gcf;

export_fig(h1,[dirNameNewNew filesep saveNameNewNew ' h1, distribution of estimated angles' saveFormat],'-transparent')
export_fig(h1_1,[dirNameNewNew filesep saveNameNewNew ' h1_1, distribution of estimated angles after theta thre' saveFormat],'-transparent')
clear h1 h1_1

locPosROINonThreROI = locPosROI;
locPosROIRefIndNonThreROI = locPosROIRefInd;
thetaROINonThreROI = thetaROI;

MROI = MROI(validThetaID,:);
thetaROI = thetaROI(validThetaID);
phiROI = phiROI(validThetaID);
rotMobilROI = rotMobilROI(validThetaID);
alphaROI = alphaROI(validThetaID);
omegaROI = omegaROI(validThetaID);
sigROI = sigROI(validThetaID);
muxROI = muxROI(validThetaID);
muyROI = muyROI(validThetaID);
muzROI = muzROI(validThetaID);
locPosROIRefPhi = locPosROIRefPhi(validThetaID);
locPosROIRefMu = locPosROIRefMu(validThetaID,:);
locPosROI = locPosROI(validThetaID,:);
locPosROIRefInd = locPosROIRefInd(validThetaID,:);
frmROI = frmROI(validThetaID);
minInd = minInd(validThetaID).';

% %% check second moment distribution of localizations with omega ~ 4*pi
% figure('Position',P1,'Visible','off');
% for mInd = 1:6
%     subplot(2,3,mInd)
%     mHist = histogram(MROI(:,mInd));
%     if mInd < 4
%         mHist.BinEdges = 0:0.05:1;
%     else
%         meanC = mean(MROI(:,mInd));
%         stdC = std(MROI(:,mInd));
%
%         leg = legend(['m=' num2str(round(meanC*1000)/1000) ',s=' num2str(round(stdC*1000)/1000)]);
%         leg.FontSize = 10;
%     end
%     xlabel('Estimated value')
%     title(['ROI' mLabel(mInd,:)])
% end
% h1_2 = gcf;
% export_fig(h1_2,[dirNameNewNew filesep saveNameNewNew ' h1_2, distribution of estimated M in ROI' saveFormat],'-transparent')
%
% isoInd = omegaROI > 12.55; % ~4*pi
%
% figure('Position',P1,'Visible','off');
% for mInd = 1:6
%     subplot(2,3,mInd)
%     mHist = histogram(MROI(isoInd,mInd));
%     if mInd < 4
%         mHist.BinEdges = 0:0.05:1;
%     else
%         meanC = mean(MROI(isoInd,mInd));
%         stdC = std(MROI(isoInd,mInd));
%
%         leg = legend(['m=' num2str(round(meanC*1000)/1000) ',s=' num2str(round(stdC*1000)/1000)]);
%         leg.FontSize = 10;
%     end
%     xlabel('Estimated value')
%     title(['ROI' mLabel(mInd,:)])
% end
% h1_3 = gcf;
% export_fig(h1_3,[dirNameNewNew filesep saveNameNewNew ' h1_3, distribution of estimated M with omega eq 4pi' saveFormat],'-transparent')
%
% clear h1_2 h1_3

%% estimation of physical angle from estimated M matrix (off ROI)
% projection to physical angles
muxOFF = nan(sum(~in),1);
muyOFF = nan(sum(~in),1);
rotMobilOFF = nan(sum(~in),1);
muzOFF = nan(sum(~in),1);
eigVOFF = nan(3,sum(~in));

locPosOFF = locPos(~in,:);
MOFF = M(~in,:);
frmOFF = frm(~in);
sigOFF = sig(~in);

for l = 1:subLoopNum
    backgTemp = backg(:,:, ((l-1)*anaLoop+1) : (l*anaLoop) );
    firstAnaLoopMInd = zeros(size(frmOFF));
    lastAnaLoopMInd = zeros(size(frmOFF));
    firstShift = 0;
    lastShift = 0;
    while sum(firstAnaLoopMInd) == 0
        firstAnaLoopMInd = frmOFF == ( (l-1)*anaLoop+1 + firstShift );
        firstAnaLoopMInd = find(firstAnaLoopMInd);
        firstShift = firstShift + 1;
    end
    firstAnaLoopMInd = firstAnaLoopMInd(1);
    while sum(lastAnaLoopMInd) == 0
        lastAnaLoopMInd = frmOFF == ( l*anaLoop - lastShift );
        lastAnaLoopMInd = find(lastAnaLoopMInd);
        lastShift = lastShift + 1;
    end
    lastAnaLoopMInd = lastAnaLoopMInd(end);
    
    parfor ii = firstAnaLoopMInd:lastAnaLoopMInd
        [muxOFF(ii),muyOFF(ii),muzOFF(ii),rotMobilOFF(ii),~]= ...
            secondM2SymmConeWeighted_v7(bx,by,Bx,By,sumNorm,MOFF(ii,:),sigOFF(ii),backgTemp(:,:,frmOFF(ii)-(l-1)*anaLoop));
%         [eigVOFF(:,ii),ind] = sort(eigVOFF(:,ii),'descend');
    end
end

clear backgTemp

muxOFF = -muxOFF;
muR = sqrt(muxOFF.^2+muyOFF.^2+muzOFF.^2);
muxOFF = muxOFF./muR;
muyOFF = muyOFF./muR;
muzOFF = muzOFF./muR;

phiOFF = atan2(muyOFF,muxOFF);%atan2(muyROIAng,muxROIAng);
phiOFF(phiOFF < 0) = phiOFF(phiOFF < 0) + pi;
muxOFF(phiOFF < 0) = -muxOFF(phiOFF < 0);
muyOFF(phiOFF < 0) = -muyOFF(phiOFF < 0);
muzOFF(phiOFF < 0) = -muzOFF(phiOFF < 0);

thetaOFF = acos(muzOFF);
% thetaOFF(abs(pi/2 - (pi-thetaOFF)) < abs(pi/2 - thetaOFF)) = pi - thetaOFF(abs(pi/2 - (pi-thetaOFF)) < abs(pi/2 - thetaOFF));

alphaOFF = acos( (  -1 + sqrt(1+8*rotMobilOFF)  )/2);
omegaOFF = 2.*4.*pi.*sin(alphaOFF/2).^2;

validThetaID1 = rad2deg(thetaOFF) >= thetaThreDeg;
validThetaID2 = rad2deg(thetaOFF) <= (180-thetaThreDeg);
validThetaIDOFF = (validThetaID1 + validThetaID2) >1;

figure('Position',P1,'Visible','off');
subplot(2,2,1)
histogram(rad2deg(alphaOFF));
xlabel('[deg]','FontSize',18)
title('OFF \alpha')
subplot(2,2,2)
histogram(rotMobilOFF);
title('OFF \gamma')
subplot(2,2,3)
histogram(rad2deg(thetaOFF));
xlabel('[deg]','FontSize',18)
title('OFF \theta')
subplot(2,2,4)
histogram(rad2deg(phiOFF));
xlabel('[deg]','FontSize',18)
title('OFF \phi')
h1_2 = gcf;

figure('Position',P1,'Visible','off');
subplot(2,2,1)
histogram(rad2deg(alphaOFF(validThetaIDOFF)));
xlabel('[deg]','FontSize',18)
title('OFF \alpha')
subplot(2,2,2)
histogram(rotMobilOFF(validThetaIDOFF));
title('OFF \gamma')
subplot(2,2,3)
histogram(rad2deg(thetaOFF(validThetaIDOFF)));
xlabel('[deg]','FontSize',18)
title('OFF \theta')
subplot(2,2,4)
histogram(rad2deg(phiOFF(validThetaIDOFF)));
xlabel('[deg]','FontSize',18)
title('OFF \phi')
h1_3 = gcf;

export_fig(h1_2,[dirNameNewNew filesep saveNameNewNew ' h1_2, distribution of estimated angles, offROI' saveFormat],'-transparent')
export_fig(h1_3,[dirNameNewNew filesep saveNameNewNew ' h1_3, distribution of estimated angles after theta thre, offROI' saveFormat],'-transparent')
clear h1_2 h1_3

MOFF = MOFF(validThetaIDOFF,:);
thetaOFF = thetaOFF(validThetaIDOFF);
phiOFF = phiOFF(validThetaIDOFF);
rotMobilOFF = rotMobilOFF(validThetaIDOFF);
alphaOFF = alphaOFF(validThetaIDOFF);
omegaOFF = omegaOFF(validThetaIDOFF);
sigOFF = sigOFF(validThetaIDOFF);
muxOFF = muxOFF(validThetaIDOFF);
muyOFF = muyOFF(validThetaIDOFF);
muzOFF = muzOFF(validThetaIDOFF);
locPosOFF = locPosOFF(validThetaIDOFF,:);
frmOFF = frmOFF(validThetaIDOFF);

clear backg
% save([dirNameNewNew filesep saveNameNewNew ' all locPos, M results' '.mat'])

%% calculate all neighbor localizations within nearNeiDis
nearNeiInd = cell(size(locPosROI,1),1);
invalidLocROI = zeros(size(nearNeiInd,1),1);
for i = 1:size(locPosROI,1)
    refLocPosROI = repmat(locPosROI(i,:),[size(locPosROI,1),1]);
    
    locPosROIDis = hypot(refLocPosROI(:,1)-locPosROI(:,1),refLocPosROI(:,2)-locPosROI(:,2));
    
    nearNeiIndTemp = (locPosROIDis < nearNeiDisNew).';
    nearNeiIndTemp(i) = 0;
    if sum(nearNeiIndTemp) < lowNei
        invalidLocROI(i) = 1;
    end
    nearNeiInd{i} = find(nearNeiIndTemp);
    
end
invalidLocROI = logical(invalidLocROI);
% nearNeiIndSize = cellfun(@(x) length(cell2mat(x)),nearNeiInd);
% nearNeiIndSize = max(nearNeiIndSize);

figure('Position',P1,'Visible','off');
hold on
hist = histogram2(locPosROI(~invalidLocROI,1),locPosROI(~invalidLocROI,2),Xedges,Yedges,'DisplayStyle','tile','ShowEmptyBins','on','EdgeColor','none');
histValues = hist.Values;
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
h0_2 = gcf;
export_fig(h0_2,[dirNameNewNew filesep saveNameNewNew ' h0_2 reconstruction with valid localizations' saveFormat],'-transparent')
clear h0_2

figure('Position',P1,'Visible','off');
hold on
scatter(locPosROI(~invalidLocROI,1),locPosROI(~invalidLocROI,2),10,'b','filled','MarkerFaceAlpha',0.3)
quiver(lineCoorAll(1,:),lineCoorAll(2,:),arrowScale*cos(slopeLineAll),arrowScale*sin(slopeLineAll),...
    'AutoScale','off','LineWidth',3,'Color','r','ShowArrowHead','on');
grid on
axis image ij
axis off
% caxis([round(minParaMeanCosPhi*10)/10 round(maxParaMeanCosPhi*10)/10])
% colormap(colormapC);
% colorbar('south')
axis([-sizeROI/2*pixelSize*shrinkFOVX sizeROI/2*pixelSize*shrinkFOVX -sizeROI/2*pixelSize*shrinkFOVY sizeROI/2*pixelSize*shrinkFOVY])
title('\phi arrow map','FontSize',20)
hold off
h0_3 = gcf;
export_fig(h0_3,[dirNameNewNew filesep saveNameNewNew ' h0_3, backbone orientation vs valid localizations' saveFormat],'-transparent')
clear h0_3

%% calculation of deviation (local heterogeneity)
muxROIMean = nan(size(locPosROI,1),1);
muyROIMean = nan(size(locPosROI,1),1);
muzROIMean = nan(size(locPosROI,1),1);
omegaROIMean = nan(size(locPosROI,1),1);

for i = 1:size(locPosROI,1)
    muxROITemp = muxROI(nearNeiInd{i});
    muyROITemp = muyROI(nearNeiInd{i});
    muzROITemp = muzROI(nearNeiInd{i});
    
    % flip the estimated each molecule orientation based on its deviation from
    % the center SM orientation
    dotMu2Center = muxROITemp.*muxROI(i) + muyROITemp.*muyROI(i) + muzROITemp.*muzROI(i);
    muxROITemp(dotMu2Center < 0) = -muxROITemp(dotMu2Center < 0);
    muyROITemp(dotMu2Center < 0) = -muyROITemp(dotMu2Center < 0);
    muzROITemp(dotMu2Center < 0) = -muzROITemp(dotMu2Center < 0);
    
    muxROIMean(i) = mean(muxROITemp);
    muyROIMean(i) = mean(muyROITemp);
    muzROIMean(i) = mean(muzROITemp);
    
    r = sqrt(muxROIMean(i)^2+muyROIMean(i)^2+muzROIMean(i)^2);
    muxROIMean(i) = muxROIMean(i)/r;
    muyROIMean(i) = muyROIMean(i)/r;
    muzROIMean(i) = muzROIMean(i)/r;
    
    omegaROIMean(i) = mean(omegaROI(nearNeiInd{i}));
end

% dot product of individual orientation vs its mean-neighbor orientation
dotMu = muxROI.*muxROIMean + muyROI.*muyROIMean + muzROI.*muzROIMean;

psiDeg = rad2deg(acos(dotMu));

% dot product of mean-neighbor orientation of an individual localization vs
% its nearest fibril backbone orientation
dotMu2Ref = locPosROIRefMu(:,1).*muxROI + locPosROIRefMu(:,2).*muyROI + locPosROIRefMu(:,3).*muzROI;
muxROITemp = muxROI;
muyROITemp = muyROI;
muzROITemp = muzROI;
muxROITemp(dotMu2Ref < 0) = -muxROITemp(dotMu2Ref < 0);
muyROITemp(dotMu2Ref < 0) = -muyROITemp(dotMu2Ref < 0);
muzROITemp(dotMu2Ref < 0) = -muzROITemp(dotMu2Ref < 0);
dotMu2Ref = locPosROIRefMu(:,1).*muxROITemp + locPosROIRefMu(:,2).*muyROITemp + locPosROIRefMu(:,3).*muzROITemp;

psi2RefDeg = rad2deg(acos(dotMu2Ref));

% in-plane version of above two
dotMuInplane = muxROI.*muxROIMean + muyROI.*muyROIMean;
muxROIMeanTemp = muxROIMean;
muyROIMeanTemp = muyROIMean;
muxROIMeanTemp(dotMuInplane < 0) = -muxROIMeanTemp(dotMuInplane < 0);
muyROIMeanTemp(dotMuInplane < 0) = -muyROIMeanTemp(dotMuInplane < 0);
dotMuInplane = muxROI.*muxROIMeanTemp + muyROI.*muyROIMeanTemp;
dotMuInplane = dotMuInplane./hypot(muxROI,muyROI)./hypot(muxROIMeanTemp,muyROIMeanTemp);

psiDegInplane = rad2deg(acos(dotMuInplane));

dotMu2RefInplane = locPosROIRefMu(:,1).*muxROI + locPosROIRefMu(:,2).*muyROI;
muxROITemp = muxROI;
muyROITemp = muyROI;
muxROITemp(dotMu2RefInplane < 0) = -muxROITemp(dotMu2RefInplane < 0);
muyROITemp(dotMu2RefInplane < 0) = -muyROITemp(dotMu2RefInplane < 0);
dotMu2RefInplane = locPosROIRefMu(:,1).*muxROITemp + locPosROIRefMu(:,2).*muyROITemp;
dotMu2RefInplane = dotMu2RefInplane./hypot(locPosROIRefMu(:,1),locPosROIRefMu(:,2))./hypot(muxROITemp,muyROITemp);

psi2RefDegInplane = rad2deg(acos(dotMu2RefInplane));

negSign = phiROI < locPosROIRefPhi;
psi2RefDegInplaneS = psi2RefDegInplane;
psi2RefDegInplaneS(negSign) = -psi2RefDegInplaneS(negSign);% with sign with respect to backbone orientation

devOmega = omegaROI - omegaROIMean;

%% Segmentation analysis
locPosLineX = cell(1,length(lineCoor));
locPosLineY = cell(1,length(lineCoor));
locPosNonThreLineX = cell(1,length(lineCoor));
locPosNonThreLineY = cell(1,length(lineCoor));
phiLine = cell(1,length(lineCoor));
thetaLine = cell(1,length(lineCoor));
psiDegInplaneLine = cell(1,length(lineCoor));
psi2RefDegInplaneLine = cell(1,length(lineCoor));
psi2RefDegInplaneSLine = cell(1,length(lineCoor));
locPosROIRefMuLine = cell(1,length(lineCoor));
omegaLine = cell(1,length(lineCoor));
sigLine = cell(1,length(lineCoor));
groupPosLineX = cell(1,length(lineCoor));
groupPosLineY = cell(1,length(lineCoor));
onTimeLine = cell(1,length(lineCoor));
photonGLine = cell(1,length(lineCoor));

psiDegInplaneLineMed = nan(1,length(lineCoor));
psi2RefDegInplaneLineMed = nan(1,length(lineCoor));
psi2RefDegInplaneSLineMed = nan(1,length(lineCoor));
omegaLineMed = nan(1,length(lineCoor));
sigLineMed = nan(1,length(lineCoor));

psiDegStdSeg = cell(1,length(lineCoor));
psi2RefDegStdSeg = cell(1,length(lineCoor));
psiDegInplaneStdSeg = cell(1,length(lineCoor));
psi2RefDegInplaneStdSeg = cell(1,length(lineCoor));
psi2RefDegInplaneSStdSeg = cell(1,length(lineCoor));
omegaROIStdSeg = cell(1,length(lineCoor));
omegaROIMedSeg = cell(1,length(lineCoor));
% absDevOmegaStdSeg = cell(1,length(lineCoor));
sigROIStdSeg = cell(1,length(lineCoor));
sigROIMedSeg = cell(1,length(lineCoor));

maxPsi2RefDegInplaneStdSeg = zeros(1,length(lineCoor));
minPsi2RefDegInplaneStdSeg = Inf*ones(1,length(lineCoor));
maxOmegaROIMedSeg = zeros(1,length(lineCoor));
minOmegaROIMedSeg = Inf*ones(1,length(lineCoor));

maxPsi2RefInplane_interPoint = cell(1,length(lineCoor));
minPsi2RefInplane_interPoint = cell(1,length(lineCoor));
maxPsi2RefInplane_locPos = cell(1,length(lineCoor));
minPsi2RefInplane_locPos = cell(1,length(lineCoor));
maxPsi2RefInplane_psi2RefInplane = cell(1,length(lineCoor));
minPsi2RefInplane_psi2RefInplane = cell(1,length(lineCoor));
maxPsi2RefInplane_psi2RefInplaneS = cell(1,length(lineCoor));
minPsi2RefInplane_psi2RefInplaneS = cell(1,length(lineCoor));
maxPsi2RefInplane_phi = cell(1,length(lineCoor));
minPsi2RefInplane_phi = cell(1,length(lineCoor));
maxPsi2RefInplane_theta = cell(1,length(lineCoor));
minPsi2RefInplane_theta = cell(1,length(lineCoor));
maxPsi2RefInplane_omega = cell(1,length(lineCoor));
minPsi2RefInplane_omega = cell(1,length(lineCoor));
maxPsi2RefInplane_sig = cell(1,length(lineCoor));
minPsi2RefInplane_sig = cell(1,length(lineCoor));

maxOmega_interPoint = cell(1,length(lineCoor));
minOmega_interPoint = cell(1,length(lineCoor));
maxOmega_locPos = cell(1,length(lineCoor));
minOmega_locPos = cell(1,length(lineCoor));
maxOmega_psi2RefInplane = cell(1,length(lineCoor));
minOmega_psi2RefInplane = cell(1,length(lineCoor));
maxOmega_psi2RefInplaneS = cell(1,length(lineCoor));
minOmega_psi2RefInplaneS = cell(1,length(lineCoor));
maxOmega_phi = cell(1,length(lineCoor));
minOmega_phi = cell(1,length(lineCoor));
maxOmega_theta = cell(1,length(lineCoor));
minOmega_theta = cell(1,length(lineCoor));
maxOmega_omega = cell(1,length(lineCoor));
minOmega_omega = cell(1,length(lineCoor));
maxOmega_sig = cell(1,length(lineCoor));
minOmega_sig = cell(1,length(lineCoor));

for lineInd = 1:length(lineCoor)
    locPosWorkInd = ~invalidLocROI & (locPosROIRefInd == lineInd);
    locPosWorkX = locPosROI(locPosWorkInd,1);
    locPosWorkY = locPosROI(locPosWorkInd,2);
    
    locPosROINonThreX = locPosROINonThreROI((locPosROIRefIndNonThreROI == lineInd),1);
    locPosROINonThreY = locPosROINonThreROI((locPosROIRefIndNonThreROI == lineInd),2);
    
    psiDegWork = psiDeg(locPosWorkInd);
    psi2RefDegWork = psi2RefDeg(locPosWorkInd);
    psiDegInplaneWork = psiDegInplane(locPosWorkInd);
    psi2RefDegInplaneWork = psi2RefDegInplane(locPosWorkInd);
    psi2RefDegInplaneSWork = psi2RefDegInplaneS(locPosWorkInd);
    locPosROIRefMuWork = locPosROIRefMu(locPosWorkInd,:);
    phiROIWork = phiROI(locPosWorkInd);
    thetaROIWork = thetaROI(locPosWorkInd);
    omegaROIWork = omegaROI(locPosWorkInd);
    sigROIWork = sigROI(locPosWorkInd);
    %     absDevOmegaWork = abs(devOmega(locPosWorkInd));
    
    locPosGroupInd = locPosGroupRefInd == lineInd;
    locPosGroupX = locPosGroup(locPosGroupInd,1);
    locPosGroupY = locPosGroup(locPosGroupInd,2);
    onTimeGroup = onTimeCorrAll(locPosGroupInd);
    photonGroup = photonCorrAll(locPosGroupInd);
    
    
    locPosLineX{lineInd} = locPosWorkX;
    locPosLineY{lineInd} = locPosWorkY;
    
    locPosNonThreLineX{lineInd} = locPosROINonThreX;
    locPosNonThreLineY{lineInd} = locPosROINonThreY;
    
    phiLine{lineInd} = phiROIWork;
    thetaLine{lineInd} = thetaROIWork;
    psiDegInplaneLine{lineInd} = psiDegInplaneWork;
    psi2RefDegInplaneLine{lineInd} = psi2RefDegInplaneWork;
    psi2RefDegInplaneSLine{lineInd} = psi2RefDegInplaneSWork;
    locPosROIRefMuLine{lineInd} = locPosROIRefMuWork;
    omegaLine{lineInd} = omegaROIWork;
    sigLine{lineInd} = sigROIWork;
    psiDegInplaneLineMed(lineInd) = median(psiDegInplaneWork);
    psi2RefDegInplaneLineMed(lineInd) = median(psi2RefDegInplaneWork);
    psi2RefDegInplaneSLineMed(lineInd) = median(psi2RefDegInplaneSWork);
    omegaLineMed(lineInd) = median(omegaROIWork);
    sigLineMed(lineInd) = median(sigROIWork);
    
    groupPosLineX{lineInd} = locPosGroupX;
    groupPosLineY{lineInd} = locPosGroupY;
    onTimeLine{lineInd} = onTimeGroup;
    photonGLine{lineInd} = photonGroup;
    
    psiDegStdSeg{lineInd} = nan(1,length(interPointX{lineInd}));
    psi2RefDegStdSeg{lineInd} = nan(1,length(interPointX{lineInd}));
    psiDegInplaneStdSeg{lineInd} = nan(1,length(interPointX{lineInd}));
    psi2RefDegInplaneStdSeg{lineInd} = nan(1,length(interPointX{lineInd}));
    psi2RefDegInplaneSStdSeg{lineInd} = nan(1,length(interPointX{lineInd}));
    omegaROIStdSeg{lineInd} = nan(1,length(interPointX{lineInd}));
    omegaROIMedSeg{lineInd} = nan(1,length(interPointX{lineInd}));
    %     absDevOmegaStdSeg{lineInd} = nan(1,length(interPointX{lineInd}));
    sigROIStdSeg{lineInd} = nan(1,length(interPointX{lineInd}));
    sigROIMedSeg{lineInd} = nan(1,length(interPointX{lineInd}));
    
    for workInd = 1:length(interPointX{lineInd})
        if workInd~=1 && workInd~=length(interPointX{lineInd})
            if interPointY{lineInd}(workInd+1) > interPointY{lineInd}(workInd)
                locPosWorkIndU =...
                    locPosWorkY <= ( -cos(interPointSlope{lineInd}(workInd))/sin(interPointSlope{lineInd}(workInd)).*locPosWorkX...
                    +(interPointY{lineInd}(workInd+1)+cos(interPointSlope{lineInd}(workInd))/sin(interPointSlope{lineInd}(workInd))*interPointX{lineInd}(workInd+1)) );
                locPosWorkIndL =...
                    locPosWorkY >= ( -cos(interPointSlope{lineInd}(workInd))/sin(interPointSlope{lineInd}(workInd)).*locPosWorkX...
                    +(interPointY{lineInd}(workInd-1)+cos(interPointSlope{lineInd}(workInd))/sin(interPointSlope{lineInd}(workInd))*interPointX{lineInd}(workInd-1)) );
            else
                locPosWorkIndL =...
                    locPosWorkY >= ( -cos(interPointSlope{lineInd}(workInd))/sin(interPointSlope{lineInd}(workInd)).*locPosWorkX...
                    +(interPointY{lineInd}(workInd+1)+cos(interPointSlope{lineInd}(workInd))/sin(interPointSlope{lineInd}(workInd))*interPointX{lineInd}(workInd+1)) );
                locPosWorkIndU =...
                    locPosWorkY <= ( -cos(interPointSlope{lineInd}(workInd))/sin(interPointSlope{lineInd}(workInd)).*locPosWorkX...
                    +(interPointY{lineInd}(workInd-1)+cos(interPointSlope{lineInd}(workInd))/sin(interPointSlope{lineInd}(workInd))*interPointX{lineInd}(workInd-1)) );
            end
            locPosWorkIndUL = locPosWorkIndU & locPosWorkIndL;
            
            if sum(locPosWorkIndUL) >= lowSeg
                psiDegStdSeg{lineInd}(workInd) = std(psiDegWork(locPosWorkIndUL));
                psi2RefDegStdSeg{lineInd}(workInd) = std(psi2RefDegWork(locPosWorkIndUL));
                psiDegInplaneStdSeg{lineInd}(workInd) = std(psiDegInplaneWork(locPosWorkIndUL));
                psi2RefDegInplaneStdSeg{lineInd}(workInd) = std(psi2RefDegInplaneWork(locPosWorkIndUL));
                psi2RefDegInplaneSStdSeg{lineInd}(workInd) = std(psi2RefDegInplaneSWork(locPosWorkIndUL));
                omegaROIStdSeg{lineInd}(workInd) = std(omegaROIWork(locPosWorkIndUL));
                omegaROIMedSeg{lineInd}(workInd) = median(omegaROIWork(locPosWorkIndUL));
                %                 absDevOmegaStdSeg{lineInd}(workInd) = std(absDevOmegaWork(locPosWorkIndUL));
                sigROIStdSeg{lineInd}(workInd) = std(sigROIWork(locPosWorkIndUL));
                sigROIMedSeg{lineInd}(workInd) = median(sigROIWork(locPosWorkIndUL));
                
                if psi2RefDegInplaneStdSeg{lineInd}(workInd) > maxPsi2RefDegInplaneStdSeg(lineInd)
                    maxPsi2RefInplane_interPoint{lineInd} = [interPointX{lineInd}(workInd) interPointY{lineInd}(workInd)];
                    maxPsi2RefInplane_locPos{lineInd} = [locPosWorkX(locPosWorkIndUL) locPosWorkY(locPosWorkIndUL)];
                    maxPsi2RefInplane_psi2RefInplane{lineInd} = psi2RefDegInplaneWork(locPosWorkIndUL);
                    maxPsi2RefInplane_psi2RefInplaneS{lineInd} = psi2RefDegInplaneSWork(locPosWorkIndUL);
                    maxPsi2RefInplane_phi{lineInd} = phiROIWork(locPosWorkIndUL);
                    maxPsi2RefInplane_theta{lineInd} = thetaROIWork(locPosWorkIndUL);
                    maxPsi2RefInplane_omega{lineInd} = omegaROIWork(locPosWorkIndUL);
                    maxPsi2RefInplane_sig{lineInd} = sigROIWork(locPosWorkIndUL);
                    maxPsi2RefDegInplaneStdSeg(lineInd) = psi2RefDegInplaneStdSeg{lineInd}(workInd);
                end
                if psi2RefDegInplaneStdSeg{lineInd}(workInd) < minPsi2RefDegInplaneStdSeg(lineInd)
                    minPsi2RefInplane_interPoint{lineInd} = [interPointX{lineInd}(workInd) interPointY{lineInd}(workInd)];
                    minPsi2RefInplane_locPos{lineInd} = [locPosWorkX(locPosWorkIndUL) locPosWorkY(locPosWorkIndUL)];
                    minPsi2RefInplane_psi2RefInplane{lineInd} = psi2RefDegInplaneWork(locPosWorkIndUL);
                    minPsi2RefInplane_psi2RefInplaneS{lineInd} = psi2RefDegInplaneSWork(locPosWorkIndUL);
                    minPsi2RefInplane_phi{lineInd} = phiROIWork(locPosWorkIndUL);
                    minPsi2RefInplane_theta{lineInd} = thetaROIWork(locPosWorkIndUL);
                    minPsi2RefInplane_omega{lineInd} = omegaROIWork(locPosWorkIndUL);
                    minPsi2RefInplane_sig{lineInd} = sigROIWork(locPosWorkIndUL);
                    minPsi2RefDegInplaneStdSeg(lineInd) = psi2RefDegInplaneStdSeg{lineInd}(workInd);
                end
                if omegaROIMedSeg{lineInd}(workInd) > maxOmegaROIMedSeg(lineInd)
                    maxOmega_interPoint{lineInd} = [interPointX{lineInd}(workInd) interPointY{lineInd}(workInd)];
                    maxOmega_locPos{lineInd} = [locPosWorkX(locPosWorkIndUL) locPosWorkY(locPosWorkIndUL)];
                    maxOmega_psi2RefInplane{lineInd} = psi2RefDegInplaneWork(locPosWorkIndUL);
                    maxOmega_psi2RefInplaneS{lineInd} = psi2RefDegInplaneSWork(locPosWorkIndUL);
                    maxOmega_phi{lineInd} = phiROIWork(locPosWorkIndUL);
                    maxOmega_theta{lineInd} = thetaROIWork(locPosWorkIndUL);
                    maxOmega_omega{lineInd} = omegaROIWork(locPosWorkIndUL);
                    maxOmega_sig{lineInd} = sigROIWork(locPosWorkIndUL);
                    maxOmegaROIMedSeg(lineInd) = omegaROIMedSeg{lineInd}(workInd);
                end
                if omegaROIMedSeg{lineInd}(workInd) < minOmegaROIMedSeg(lineInd)
                    minOmega_interPoint{lineInd} = [interPointX{lineInd}(workInd) interPointY{lineInd}(workInd)];
                    minOmega_locPos{lineInd} = [locPosWorkX(locPosWorkIndUL) locPosWorkY(locPosWorkIndUL)];
                    minOmega_psi2RefInplane{lineInd} = psi2RefDegInplaneWork(locPosWorkIndUL);
                    minOmega_psi2RefInplaneS{lineInd} = psi2RefDegInplaneSWork(locPosWorkIndUL);
                    minOmega_phi{lineInd} = phiROIWork(locPosWorkIndUL);
                    minOmega_theta{lineInd} = thetaROIWork(locPosWorkIndUL);
                    minOmega_omega{lineInd} = omegaROIWork(locPosWorkIndUL);
                    minOmega_sig{lineInd} = sigROIWork(locPosWorkIndUL);
                    minOmegaROIMedSeg(lineInd) = omegaROIMedSeg{lineInd}(workInd);
                end
            end
        end
    end
end

%% for phi deviation
[~,rgbPlotPsiDegMean,~,~,colorBarPsiDegMean,~,...
    ~,~,minParaMeanPsiDeg,maxParaMeanPsiDeg,~,~] = ...
    reconstVLocNumHAddPara_v2(locPosROI(~invalidLocROI,1),locPosROI(~invalidLocROI,2),  psiDeg(~invalidLocROI)  ,[],[],binSize,...
    Xedges,Yedges,max(max(histValues)),indN,colorMapAng);
[~,rgbPlotPsi2RefDegMean,~,~,colorBarPsi2RefDegMean,~,...
    ~,~,minParaMeanPsi2RefDeg,maxParaMeanPsi2RefDeg,~,~] = ...
    reconstVLocNumHAddPara_v2(locPosROI(~invalidLocROI,1),locPosROI(~invalidLocROI,2),  psi2RefDeg(~invalidLocROI)  ,[],[],binSize,...
    Xedges,Yedges,max(max(histValues)),indN,colorMapAng);

% for inplane phi deviation
[~,rgbPlotPsiDegInplaneMean,~,~,colorBarPsiDegInplaneMean,~,...
    ~,~,minParaMeanPsiDegInplane,maxParaMeanPsiDegInplane,~,~] = ...
    reconstVLocNumHAddPara_v2(locPosROI(~invalidLocROI,1),locPosROI(~invalidLocROI,2),  psiDegInplane(~invalidLocROI)  ,[],[],binSize,...
    Xedges,Yedges,max(max(histValues)),indN,colorMapAng);
[~,rgbPlotPsi2RefDegInplaneMean,~,~,colorBarPsi2RefDegInplaneMean,~,...
    ~,~,minParaMeanPsi2RefDegInplane,maxParaMeanPsi2RefDegInplane,~,~] = ...
    reconstVLocNumHAddPara_v2(locPosROI(~invalidLocROI,1),locPosROI(~invalidLocROI,2),  psi2RefDegInplane(~invalidLocROI)  ,[],[],binSize,...
    Xedges,Yedges,max(max(histValues)),indN,colorMapAng);

% for omega deviation
[~,rgbPlotOmegaMean,rgbPlotOmegaStd,~,colorBarOmegaMean,colorBarOmegaStd,...
    ~,~,minParaMeanOmega,maxParaMeanOmega,minParaStdOmega,maxParaStdOmega] = ...
    reconstVLocNumHAddPara_v2(locPosROI(~invalidLocROI,1),locPosROI(~invalidLocROI,2),  omegaROI(~invalidLocROI)  ,[],[],binSize,...
    Xedges,Yedges,max(max(histValues)),indN,colorMapAng);
[~,rgbPlotAbsDevOmegaMean,~,~,colorBarAbsDevOmegaMean,~,...
    ~,~,minParaMeanAbsDevOmega,maxParaMeanAbsDevOmega,~,~] = ...
    reconstVLocNumHAddPara_v2(locPosROI(~invalidLocROI,1),locPosROI(~invalidLocROI,2),  abs(devOmega(~invalidLocROI))  ,[],[],binSize,...
    Xedges,Yedges,max(max(histValues)),indN,colorMapAng);

% for tau_on
[~,rgbPlotTauMean,~,~,colorBarTauMean,~,...
    ~,~,minParaMeanTau,maxParaMeanTau,~,~] = ...
    reconstVLocNumHAddPara_v2(locPosGroup(:,1),locPosGroup(:,2),  onTimeCorrAll  ,[],[0.1],binSize,...
    Xedges,Yedges,max(max(histValuesG)),indNG,colorMapAng);

% for brightness per molecule
[~,rgbPlotPhotonMean,~,~,colorBarPhotonMean,~,...
    ~,~,minParaMeanPhoton,maxParaMeanPhoton,~,~] = ...
    reconstVLocNumHAddPara_v2(locPosGroup(:,1),locPosGroup(:,2),  photonCorrAll  ,[],[],binSize,...
    Xedges,Yedges,max(max(histValuesG)),indNG,colorMapAng);

% for brightness HZ
[~,rgbPlotHzMean,~,~,colorBarHzMean,~,...
    ~,~,minParaMeanHz,maxParaMeanHz,~,~] = ...
    reconstVLocNumHAddPara_v2(locPosGroup(:,1),locPosGroup(:,2),  brightnessHzAll  ,[],[],binSize,...
    Xedges,Yedges,max(max(histValuesG)),indNG,colorMapAng);

imgSize = size(rgbPlotPsiDegMean);
%% visualization
subplotTight = @(m,n,p)subtightplot(m, n, p, [0.001 0.001], [0.00 0.1], [0 0]);
%% phi plot
figure('Position',P1,'Visible','off');
ax1 = subplotTight(3,4,[1 5]);
imagesc(rgbPlotPsiDegMean(floor(imgSize(1)/2-imgSize(1)*shrinkFOVY/2)+1:floor(imgSize(1)/2+imgSize(1)*shrinkFOVY/2),...
    floor(imgSize(2)/2-imgSize(2)*shrinkFOVX/2)+1:floor(imgSize(2)/2+imgSize(2)*shrinkFOVX/2),: ))
axis image ij
axis off
title('mean of \psi map','FontSize',20)
co = colorbar('southoutside');
ax5 = subplotTight(3,4,9);
imagesc(rot90(colorBarPsiDegMean,2))
axis ij;
ax = gca;
ax.YTick = [ax.YLim(1) ax.YLim(2)];
ax.XTick = [ax.XLim(1) (ax.XLim(2)+ax.XLim(1))/2 ax.XLim(2)];
ax.YTickLabel= ({['\geq' num2str(indN)],0});
ax.XTickLabel = [round(minParaMeanPsiDeg*10)/10 round((minParaMeanPsiDeg+maxParaMeanPsiDeg)/2*10)/10 round(maxParaMeanPsiDeg*10)/10];
ax.Position(1) = ax.Position(1)+ax.Position(3)/8;
ax.Position(3) = ax.Position(3)*3/4;
ax.Position(2) = co.Position(2) - 4*co.Position(4);
ax.Position(4) = 2*co.Position(4);

ax2 = subplotTight(3,4,[2 6]);
hold on
imagesc(rgbPlotPsi2RefDegMean( floor(imgSize(1)/2-imgSize(1)*shrinkFOVY/2)+1:floor(imgSize(1)/2+imgSize(1)*shrinkFOVY/2),...
    floor(imgSize(2)/2-imgSize(2)*shrinkFOVX/2)+1:floor(imgSize(2)/2+imgSize(2)*shrinkFOVX/2),: ))
axis image ij
ax = gca;
rectangle('Position',[ax.XLim(2)-scaleBar/binSize-(ax.XLim(2)-ax.XLim(1))/20,...
    ax.YLim(2)-1-(ax.YLim(2)-ax.YLim(1))/20,...
    scaleBar/binSize,1],'FaceColor','w','EdgeColor','w');
axis off
title('mean of \psi_{ref} map','FontSize',20)
hold off
ax6 = subplotTight(3,4,10);
imagesc(rot90(colorBarPsi2RefDegMean,2))
axis ij;
ax = gca;
ax.YTick = [ax.YLim(1) ax.YLim(2)];
ax.XTick = [ax.XLim(1) (ax.XLim(2)+ax.XLim(1))/2 ax.XLim(2)];
ax.YTickLabel= ({['\geq' num2str(indN)],0});
ax.XTickLabel = [round(minParaMeanPsi2RefDeg*10)/10 round((minParaMeanPsi2RefDeg+maxParaMeanPsi2RefDeg)/2*10)/10 round(maxParaMeanPsi2RefDeg*10)/10];
ax.Position(1) = ax.Position(1)+ax.Position(3)/8;
ax.Position(3) = ax.Position(3)*3/4;
ax.Position(2) = co.Position(2) - 4*co.Position(4);
ax.Position(4) = 2*co.Position(4);

ax3 = subplotTight(3,4,[3 7]);
hold on
for n = 1:numOfChunk
    plot((xROI{n}-0.5).*binSize-(((sizeROI-1)/2+0.5)*pixelSize),(yROI{n}-0.5).*binSize-(((sizeROI-1)/2+0.5)*pixelSize))
end
nearNeiIndTemp = 1:size(locPosROI,1);
nearNeiIndTemp(nearNeiInd{plotLocInd}) = [];
scatter(locPosROI(nearNeiIndTemp,1),locPosROI(nearNeiIndTemp,2),'bo')
scatter(locPosROI(nearNeiInd{plotLocInd},1),locPosROI(nearNeiInd{plotLocInd},2),'r+')
scatter(locPosROI(plotLocInd,1),locPosROI(plotLocInd,2),'k*')
grid on
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'Box', 'on');
axis image ij
axis([Xedges(floor(length(Xedges)/2-length(Xedges)*shrinkFOVX/2)+1) Xedges(floor(length(Xedges)/2+length(Xedges)*shrinkFOVX/2))...
    Yedges(floor(length(Yedges)/2-length(Yedges)*shrinkFOVY/2)+1) Yedges(floor(length(Yedges)/2+length(Yedges)*shrinkFOVY/2))])
title({['Nearby localizations, locInd=' num2str(plotLocInd)]; ['red +: localizations nearby, blue o: otherwise']},'FontSize',16)
hold off

ax4 = subplotTight(3,4,[4 8]);
hold on
for n = 1:numOfChunk
    plot((xROI{n}-0.5).*binSize-(((sizeROI-1)/2+0.5)*pixelSize),(yROI{n}-0.5).*binSize-(((sizeROI-1)/2+0.5)*pixelSize))
end
% nearNeiIndTemp = 1:size(locPosROI,1);
% nearNeiIndTemp(nearNeiInd{plotLocInd}) = [];
scatter(locPosROI(nearNeiIndTemp,1),locPosROI(nearNeiIndTemp,2),'bo')
scatter(locPosROI(nearNeiInd{plotLocInd},1),locPosROI(nearNeiInd{plotLocInd},2),'r+')
scatter(locPosROI(plotLocInd,1),locPosROI(plotLocInd,2),'k*')
grid on
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'Box', 'on');
axis image ij
axis([(locPosROI(plotLocInd,1)-200) (locPosROI(plotLocInd,1)+200)...
    (locPosROI(plotLocInd,2)-200) (locPosROI(plotLocInd,2)+200)])
ax = gca;
rec = rectangle('Position',[ax.XLim(2)-scaleBarMin-(ax.XLim(2)-ax.XLim(1))/20,...
    ax.YLim(2)-1-(ax.YLim(2)-ax.YLim(1))/20,...
    scaleBarMin,3],'FaceColor','k','EdgeColor','k');
title({['Nearby localizations, locInd=' num2str(plotLocInd)]; ['red +: localizations nearby, blue o: otherwise']},'FontSize',16)
hold off

colorbar(ax1,'off');
h2 = gcf;
export_fig(h2,[dirNameNewNew filesep saveNameNewNew ' h2, rgb plot, phi deviation' saveFormat],'-transparent')

%% phi scatter plot
figure('Position',P1,'Visible','off');
ax1 = subplotTight(2,4,[1 2 5 6]);
hold on
scatter(locPosROI(~invalidLocROI,1),locPosROI(~invalidLocROI,2),10,psiDeg(~invalidLocROI),'filled','MarkerFaceAlpha',0.3)
grid on
axis image ij
axis off
caxis([round(minParaMeanPsiDeg*10)/10 round(maxParaMeanPsiDeg*10)/10])
colormap(flipud(colormapC));
colorbar('south')
axis([-sizeROI/2*pixelSize*shrinkFOVX sizeROI/2*pixelSize*shrinkFOVX -sizeROI/2*pixelSize*shrinkFOVY sizeROI/2*pixelSize*shrinkFOVY])
title('\psi map','FontSize',20)
hold off

ax2 = subplotTight(2,4,[3 4 7 8]);
hold on
scatter(locPosROI(~invalidLocROI,1),locPosROI(~invalidLocROI,2),10,psi2RefDeg(~invalidLocROI),'filled','MarkerFaceAlpha',0.3)
grid on
axis image ij
axis off
caxis([round(minParaMeanPsi2RefDeg*10)/10 round(maxParaMeanPsi2RefDeg*10)/10])
colormap(flipud(colormapC));
colorbar('south')
axis([-sizeROI/2*pixelSize*shrinkFOVX sizeROI/2*pixelSize*shrinkFOVX -sizeROI/2*pixelSize*shrinkFOVY sizeROI/2*pixelSize*shrinkFOVY])
ax = gca;
rectangle('Position',[ax.XLim(2)-scaleBar-(ax.XLim(2)-ax.XLim(1))/20,...
    ax.YLim(2)-1-(ax.YLim(2)-ax.YLim(1))/20,...
    scaleBar,3],'FaceColor','k','EdgeColor','k');
title('\psi_{ref} map','FontSize',20)
hold off
h2_1 = gcf;
export_fig(h2_1,[dirNameNewNew filesep saveNameNewNew ' h2_1, scatter plot, phi deviation' saveFormat],'-transparent')

% sorted version
locPosROIValid = locPosROI(~invalidLocROI,:);
psiDegValid = psiDeg(~invalidLocROI);
psi2RefDegValid = psi2RefDeg(~invalidLocROI);
[psiDegValidSort,psiDegValidSortI] = sort(psiDegValid);
[psi2RefDegValidSort,psi2RefDegValidSortI] = sort(psi2RefDegValid);

figure('Position',P1,'Visible','off');
ax1 = subplotTight(2,4,[1 2 5 6]);
hold on
scatter(locPosROIValid(psiDegValidSortI,1),locPosROIValid(psiDegValidSortI,2),10,psiDegValidSort,'filled','MarkerFaceAlpha',0.3)
grid on
axis image ij
axis off
caxis([round(minParaMeanPsiDeg*10)/10 round(maxParaMeanPsiDeg*10)/10])
colormap(flipud(colormapC));
colorbar('south')
axis([-sizeROI/2*pixelSize*shrinkFOVX sizeROI/2*pixelSize*shrinkFOVX -sizeROI/2*pixelSize*shrinkFOVY sizeROI/2*pixelSize*shrinkFOVY])
title('\psi map, sorted','FontSize',20)
hold off

ax2 = subplotTight(2,4,[3 4 7 8]);
hold on
scatter(locPosROIValid(psi2RefDegValidSortI,1),locPosROIValid(psi2RefDegValidSortI,2),10,psi2RefDegValidSort,'filled','MarkerFaceAlpha',0.3)
grid on
axis image ij
axis off
caxis([round(minParaMeanPsi2RefDeg*10)/10 round(maxParaMeanPsi2RefDeg*10)/10])
colormap(flipud(colormapC));
colorbar('south')
axis([-sizeROI/2*pixelSize*shrinkFOVX sizeROI/2*pixelSize*shrinkFOVX -sizeROI/2*pixelSize*shrinkFOVY sizeROI/2*pixelSize*shrinkFOVY])
ax = gca;
rectangle('Position',[ax.XLim(2)-scaleBar-(ax.XLim(2)-ax.XLim(1))/20,...
    ax.YLim(2)-1-(ax.YLim(2)-ax.YLim(1))/20,...
    scaleBar,3],'FaceColor','k','EdgeColor','k');
title('\psi_{ref} map, sorted','FontSize',20)
hold off
h2_1_1 = gcf;
export_fig(h2_1_1,[dirNameNewNew filesep saveNameNewNew ' h2_1_1, scatter plot, phi deviation, sorted' saveFormat],'-transparent')

% %% Segmented plot
% figure('Position',P1,'Visible','on')
% ax1 = subplotTight(2,4,[1 2 5 6]);
% hold on
% for lineInd = 1:length(lineCoor)
%     plot(interPointX{lineInd},interPointY{lineInd},'o','MarkerSize',5)
%     quiver(interPointX{lineInd},interPointY{lineInd},arrowScale*cos(interPointSlope{lineInd}),arrowScale*sin(interPointSlope{lineInd}),...
%     'AutoScale','off','LineWidth',3,'Color','r','ShowArrowHead','on');
%     plot(lineCoor{lineInd}(1,:),lineCoor{lineInd}(2,:),'+','MarkerSize',10)
% end
%
% scatter(locPosROIValid(psiDegValidSortI,1),locPosROIValid(psiDegValidSortI,2),10,psiDegValidSort,'filled','MarkerFaceAlpha',0.3)
% grid on
% axis image ij
% axis off
% caxis([round(minParaMeanPsiDeg*10)/10 round(maxParaMeanPsiDeg*10)/10])
% colormap(flipud(colormapC));
% colorbar('south')
% axis([-sizeROI/2*pixelSize*shrinkFOVX sizeROI/2*pixelSize*shrinkFOVX -sizeROI/2*pixelSize*shrinkFOVY sizeROI/2*pixelSize*shrinkFOVY])
% title('\psi map, sorted','FontSize',20)
% hold off
%
% ax2 = subplotTight(2,4,[3 4 7 8]);
% hold on
% scatter(locPosROIValid(psi2RefDegValidSortI,1),locPosROIValid(psi2RefDegValidSortI,2),10,psi2RefDegValidSort,'filled','MarkerFaceAlpha',0.3)
% grid on
% axis image ij
% axis off
% caxis([round(minParaMeanPsi2RefDeg*10)/10 round(maxParaMeanPsi2RefDeg*10)/10])
% colormap(flipud(colormapC));
% colorbar('south')
% axis([-sizeROI/2*pixelSize*shrinkFOVX sizeROI/2*pixelSize*shrinkFOVX -sizeROI/2*pixelSize*shrinkFOVY sizeROI/2*pixelSize*shrinkFOVY])
% ax = gca;
% rectangle('Position',[ax.XLim(2)-scaleBar-(ax.XLim(2)-ax.XLim(1))/20,...
%     ax.YLim(2)-1-(ax.YLim(2)-ax.YLim(1))/20,...
%     scaleBar,3],'FaceColor','k','EdgeColor','k');
% title('\psi_{ref} map, sorted','FontSize',20)
% hold off
% h2_1_2 = gcf;
% export_fig(h2_1_2,[dirNameNewNew filesep saveNameNewNew ' h2_1_1, scatter plot, phi deviation, sorted' saveFormat],'-transparent')
%

%% inplane phi plot
figure('Position',P1,'Visible','off');
ax1 = subplotTight(3,4,[1 5]);
imagesc(rgbPlotPsiDegInplaneMean( floor(imgSize(1)/2-imgSize(1)*shrinkFOVY/2)+1:floor(imgSize(1)/2+imgSize(1)*shrinkFOVY/2),...
    floor(imgSize(2)/2-imgSize(2)*shrinkFOVX/2)+1:floor(imgSize(2)/2+imgSize(2)*shrinkFOVX/2),: ))
axis image ij
axis off
title('mean of \psi^{inplane} map','FontSize',20)
co = colorbar('southoutside');
ax5 = subplotTight(3,4,9);
imagesc(rot90(colorBarPsiDegInplaneMean,2))
axis ij;
ax = gca;
ax.YTick = [ax.YLim(1) ax.YLim(2)];
ax.XTick = [ax.XLim(1) (ax.XLim(2)+ax.XLim(1))/2 ax.XLim(2)];
ax.YTickLabel= ({['\geq' num2str(indN)],0});
ax.XTickLabel = [round(minParaMeanPsiDegInplane*10)/10 round((minParaMeanPsiDegInplane+maxParaMeanPsiDegInplane)/2*10)/10 round(maxParaMeanPsiDegInplane*10)/10];
ax.Position(1) = ax.Position(1)+ax.Position(3)/8;
ax.Position(3) = ax.Position(3)*3/4;
ax.Position(2) = co.Position(2) - 4*co.Position(4);
ax.Position(4) = 2*co.Position(4);

ax2 = subplotTight(3,4,[2 6]);
hold on
imagesc(rgbPlotPsi2RefDegInplaneMean( floor(imgSize(1)/2-imgSize(1)*shrinkFOVY/2)+1:floor(imgSize(1)/2+imgSize(1)*shrinkFOVY/2),...
    floor(imgSize(2)/2-imgSize(2)*shrinkFOVX/2)+1:floor(imgSize(2)/2+imgSize(2)*shrinkFOVX/2),: ))
axis image ij
ax = gca;
rec = rectangle('Position',[ax.XLim(2)-scaleBar/binSize-(ax.XLim(2)-ax.XLim(1))/20,...
    ax.YLim(2)-1-(ax.YLim(2)-ax.YLim(1))/20,...
    scaleBar/binSize,1],'FaceColor','w','EdgeColor','w');
axis off
title('mean of \psi^{inplane}_{ref} map','FontSize',20)
hold off
ax6 = subplotTight(3,4,10);
imagesc(rot90(colorBarPsi2RefDegInplaneMean,2))
axis ij;
ax = gca;
ax.YTick = [ax.YLim(1) ax.YLim(2)];
ax.XTick = [ax.XLim(1) (ax.XLim(2)+ax.XLim(1))/2 ax.XLim(2)];
ax.YTickLabel= ({['\geq' num2str(indN)],0});
ax.XTickLabel = [round(minParaMeanPsi2RefDegInplane*10)/10 round((minParaMeanPsi2RefDegInplane+maxParaMeanPsi2RefDegInplane)/2*10)/10 round(maxParaMeanPsi2RefDegInplane*10)/10];
ax.Position(1) = ax.Position(1)+ax.Position(3)/8;
ax.Position(3) = ax.Position(3)*3/4;
ax.Position(2) = co.Position(2) - 4*co.Position(4);
ax.Position(4) = 2*co.Position(4);

ax3 = subplotTight(3,4,[3 7]);
hold on
for n = 1:numOfChunk
    plot((xROI{n}-0.5).*binSize-(((sizeROI-1)/2+0.5)*pixelSize),(yROI{n}-0.5).*binSize-(((sizeROI-1)/2+0.5)*pixelSize))
end
% nearNeiIndTemp = 1:size(locPosROI,1);
% nearNeiIndTemp(nearNeiInd{plotLocInd}) = [];
scatter(locPosROI(nearNeiIndTemp,1),locPosROI(nearNeiIndTemp,2),'bo')
scatter(locPosROI(nearNeiInd{plotLocInd},1),locPosROI(nearNeiInd{plotLocInd},2),'r+')
scatter(locPosROI(plotLocInd,1),locPosROI(plotLocInd,2),'k*')
grid on
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'Box', 'on');
axis image ij
axis([Xedges(floor(length(Xedges)/2-length(Xedges)*shrinkFOVX/2)+1) Xedges(floor(length(Xedges)/2+length(Xedges)*shrinkFOVX/2))...
    Yedges(floor(length(Yedges)/2-length(Yedges)*shrinkFOVY/2)+1) Yedges(floor(length(Yedges)/2+length(Yedges)*shrinkFOVY/2))])
title({['Nearby localizations, locInd=' num2str(plotLocInd)]; ['red +: localizations nearby, blue o: otherwise']},'FontSize',16)
hold off

ax4 = subplotTight(3,4,[4 8]);
hold on
for n = 1:numOfChunk
    plot((xROI{n}-0.5).*binSize-(((sizeROI-1)/2+0.5)*pixelSize),(yROI{n}-0.5).*binSize-(((sizeROI-1)/2+0.5)*pixelSize))
end
% nearNeiIndTemp = 1:size(locPosROI,1);
% nearNeiIndTemp(nearNeiInd{plotLocInd}) = [];
scatter(locPosROI(nearNeiIndTemp,1),locPosROI(nearNeiIndTemp,2),'bo')
scatter(locPosROI(nearNeiInd{plotLocInd},1),locPosROI(nearNeiInd{plotLocInd},2),'r+')
scatter(locPosROI(plotLocInd,1),locPosROI(plotLocInd,2),'k*')
grid on
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'Box', 'on');
axis image ij
axis([(locPosROI(plotLocInd,1)-200) (locPosROI(plotLocInd,1)+200)...
    (locPosROI(plotLocInd,2)-200) (locPosROI(plotLocInd,2)+200)])
ax = gca;
rec = rectangle('Position',[ax.XLim(2)-scaleBarMin-(ax.XLim(2)-ax.XLim(1))/20,...
    ax.YLim(2)-1-(ax.YLim(2)-ax.YLim(1))/20,...
    scaleBarMin,3],'FaceColor','k','EdgeColor','k');
title({['Nearby localizations, locInd=' num2str(plotLocInd)]; ['red +: localizations nearby, blue o: otherwise']},'FontSize',16)
hold off
colorbar(ax1,'off');
h2_2 = gcf;

export_fig(h2_2,[dirNameNewNew filesep saveNameNewNew ' h2_2, rgb plot, inplane phi deviation' saveFormat],'-transparent')

%% inplane phi scatter plot
figure('Position',P1,'Visible','off');
ax1 = subplotTight(2,4,[1 2 5 6]);
hold on
scatter(locPosROI(~invalidLocROI,1),locPosROI(~invalidLocROI,2),10,psiDegInplane(~invalidLocROI),'filled','MarkerFaceAlpha',0.3)
grid on
axis image ij
axis off
caxis([round(minParaMeanPsiDegInplane*10)/10 round(maxParaMeanPsiDegInplane*10)/10])
colormap(flipud(colormapC));
colorbar('south')
axis([-sizeROI/2*pixelSize*shrinkFOVX sizeROI/2*pixelSize*shrinkFOVX -sizeROI/2*pixelSize*shrinkFOVY sizeROI/2*pixelSize*shrinkFOVY])
title('\psi^{inplane} map','FontSize',20)
hold off

ax2 = subplotTight(2,4,[3 4 7 8]);
hold on
scatter(locPosROI(~invalidLocROI,1),locPosROI(~invalidLocROI,2),10,psi2RefDegInplane(~invalidLocROI),'filled','MarkerFaceAlpha',0.3)
grid on
axis image ij
axis off
caxis([round(minParaMeanPsi2RefDegInplane*10)/10 round(maxParaMeanPsi2RefDegInplane*10)/10])
colormap(flipud(colormapC));
colorbar('south')
axis([-sizeROI/2*pixelSize*shrinkFOVX sizeROI/2*pixelSize*shrinkFOVX -sizeROI/2*pixelSize*shrinkFOVY sizeROI/2*pixelSize*shrinkFOVY])
ax = gca;
rec = rectangle('Position',[ax.XLim(2)-scaleBar-(ax.XLim(2)-ax.XLim(1))/20,...
    ax.YLim(2)-1-(ax.YLim(2)-ax.YLim(1))/20,...
    scaleBar,3],'FaceColor','k','EdgeColor','k');
title('\psi^{inplane}_{ref} map','FontSize',20)
hold off
h2_3 = gcf;
export_fig(h2_3,[dirNameNewNew filesep saveNameNewNew ' h2_3, scatter plot, inplane phi deviation' saveFormat],'-transparent')

% sorted version
psiDegInplaneValid = psiDegInplane(~invalidLocROI);
psi2RefDegInplaneValid = psi2RefDegInplane(~invalidLocROI);
[psiDegInplaneValidSort,psiDegInplaneValidSortI] = sort(psiDegInplaneValid);
[psi2RefDegInplaneValidSort,psi2RefDegInplaneValidSortI] = sort(psi2RefDegInplaneValid);

figure('Position',P1,'Visible','on');
ax1 = subplotTight(2,4,[1 2 5 6]);
hold on
scatter(locPosROIValid(psiDegInplaneValidSortI,1),locPosROIValid(psiDegInplaneValidSortI,2),10,psiDegInplaneValidSort,'filled','MarkerFaceAlpha',0.3)
grid on
axis image ij
axis off
caxis([round(minParaMeanPsiDegInplane*10)/10 round(maxParaMeanPsiDegInplane*10)/10])
colormap(flipud(colormapC));
colorbar('south')
axis([-sizeROI/2*pixelSize*shrinkFOVX sizeROI/2*pixelSize*shrinkFOVX -sizeROI/2*pixelSize*shrinkFOVY sizeROI/2*pixelSize*shrinkFOVY])
title('\psi^{inplane} map, sorted','FontSize',20)
hold off

ax2 = subplotTight(2,4,[3 4 7 8]);
hold on
scatter(locPosROIValid(psi2RefDegInplaneValidSortI,1),locPosROIValid(psi2RefDegInplaneValidSortI,2),10,psi2RefDegInplaneValidSort,'filled','MarkerFaceAlpha',0.3)
grid on
axis image ij
axis off
caxis([round(minParaMeanPsi2RefDegInplane*10)/10 round(maxParaMeanPsi2RefDegInplane*10)/10])
colormap(flipud(colormapC));
colorbar('south')
axis([-sizeROI/2*pixelSize*shrinkFOVX sizeROI/2*pixelSize*shrinkFOVX -sizeROI/2*pixelSize*shrinkFOVY sizeROI/2*pixelSize*shrinkFOVY])
ax = gca;
rec = rectangle('Position',[ax.XLim(2)-scaleBar-(ax.XLim(2)-ax.XLim(1))/20,...
    ax.YLim(2)-1-(ax.YLim(2)-ax.YLim(1))/20,...
    scaleBar,3],'FaceColor','k','EdgeColor','k');
title('\psi^{inplane}_{ref} map, sorted','FontSize',20)
hold off
h2_3_1 = gcf;
export_fig(h2_3_1,[dirNameNewNew filesep saveNameNewNew ' h2_3_1, scatter plot, inplane phi deviation, sorted' saveFormat],'-transparent')


figure('Position',P1,'Visible','on');
subplot(2,1,1);
hold on
locPosLineValidX = locPosLineX{visLineInd(1)};
locPosLineValidY = locPosLineY{visLineInd(1)};
psi2RefDegInplaneLineValid = psi2RefDegInplaneLine{visLineInd(1)};
[psi2RefDegInplaneLineValidSort,psi2RefDegInplaneLineValidSortI] = sort(psi2RefDegInplaneLineValid);
scatter(locPosLineValidX(psi2RefDegInplaneLineValidSortI),locPosLineValidY(psi2RefDegInplaneLineValidSortI),15,psi2RefDegInplaneLineValidSort,'filled','MarkerFaceAlpha',0.5)
% plot(interPointX{visLineInd(1)},interPointY{visLineInd(1)},'+','MarkerSize',5,'Color',[0.5 0.5 0.5])
% hold on
% visRot = rad2deg(mean(interPointSlope{visLineInd(1)}));
% psi2RefInplanePlot = psiDegInplaneLine{visLineInd(1)};
% psi2RefInplanePlot(psi2RefInplanePlot > maxDegPlot) = maxDegPlot;
% psi2RefInplaneColor = round(psi2RefInplanePlot/maxDegPlot*(size(colormapArrow,1)-1))+1;
% for j = 1:length(phiLine{visLineInd(1)})
% quiver(locPosLineX{visLineInd(1)}(j),locPosLineY{visLineInd(1)}(j),arrowScaleNoHead*cos(phiLine{visLineInd(1)}(j)),arrowScaleNoHead*sin(phiLine{visLineInd(1)}(j)),...
%     'AutoScale','off','LineWidth',2,'Color',colormapArrow(psi2RefInplaneColor(j),:),'ShowArrowHead','off');
% end
axis image ij
axis off
% if isempty(visRot)
visRot = rad2deg(mean(interPointSlope{visLineInd(1)}));
% end
camroll(visRot)
camzoom(visZoom)
caxis([round(minParaMeanPsi2RefDegInplane*10)/10 round(maxParaMeanPsi2RefDegInplane*10)/10])
colormap(flipud(colormapC));
colorbar('south')
% axis([-sizeROI/2*pixelSize*shrinkFOVX sizeROI/2*pixelSize*shrinkFOVX -sizeROI/2*pixelSize*shrinkFOVY sizeROI/2*pixelSize*shrinkFOVY])
% title(['\psi^{inplane}_{ref} map, BB=' num2str(visLineInd(1))],'FontSize',20)
hold off

subplot(2,1,2);
hold on
locPosLineValidX = locPosLineX{visLineInd(2)};
locPosLineValidY = locPosLineY{visLineInd(2)};
psi2RefDegInplaneLineValid = psi2RefDegInplaneLine{visLineInd(2)};
[psi2RefDegInplaneLineValidSort,psi2RefDegInplaneLineValidSortI] = sort(psi2RefDegInplaneLineValid);
scatter(locPosLineValidX(psi2RefDegInplaneLineValidSortI),locPosLineValidY(psi2RefDegInplaneLineValidSortI),15,psi2RefDegInplaneLineValidSort,'filled','MarkerFaceAlpha',0.5)
% plot(interPointX{visLineInd(2)},interPointY{visLineInd(2)},'+','MarkerSize',5,'Color',[0.5 0.5 0.5])
% hold on
% psi2RefInplanePlot = psiDegInplaneLine{visLineInd(2)};
% psi2RefInplanePlot(psi2RefInplanePlot > maxDegPlot) = maxDegPlot;
% psi2RefInplaneColor = round(psi2RefInplanePlot/maxDegPlot*(size(colormapArrow,1)-1))+1;
% for j = 1:length(phiLine{visLineInd(2)})
% quiver(locPosLineX{visLineInd(2)}(j),locPosLineY{visLineInd(2)}(j),arrowScaleNoHead*cos(phiLine{visLineInd(2)}(j)),arrowScaleNoHead*sin(phiLine{visLineInd(2)}(j)),...
%     'AutoScale','off','LineWidth',2,'Color',colormapArrow(psi2RefInplaneColor(j),:),'ShowArrowHead','off');
% end
axis image ij
axis off
% if isempty(visRot)
visRot = rad2deg(mean(interPointSlope{visLineInd(2)}));
% end
camroll(visRot)
camzoom(visZoom)
caxis([round(minParaMeanPsi2RefDegInplane*10)/10 round(maxParaMeanPsi2RefDegInplane*10)/10])
colormap(flipud(colormapC));
colorbar('south')
% axis([-sizeROI/2*pixelSize*shrinkFOVX sizeROI/2*pixelSize*shrinkFOVX -sizeROI/2*pixelSize*shrinkFOVY sizeROI/2*pixelSize*shrinkFOVY])
ax = gca;
rec = rectangle('Position',[ax.XLim(2)-scaleBar-(ax.XLim(2)-ax.XLim(1))/20,...
    ax.YLim(2)-1-(ax.YLim(2)-ax.YLim(1))/20,...
    scaleBar,3],'FaceColor','k','EdgeColor','k');
% title(['\psi^{inplane}_{ref} map, BB=' num2str(visLineInd(2))],'FontSize',20)
hold off
h2_3_2 = gcf;
export_fig(h2_3_2,[dirNameNewNew filesep saveNameNewNew ' h2_3_2, scatter plot, inplane phi deviation, 2 BBs' saveFormat],'-transparent')

% figure('Position',P1,'Visible','on');
% subplot(2,1,1);
% hold on
% locPosLineValidX = locPosLineX{visLineInd(1)};
% locPosLineValidY = locPosLineY{visLineInd(1)};
% phiLineValid = phiLine{visLineInd(1)};
% psi2RefDegInplaneLineValid = psi2RefDegInplaneLine{visLineInd(1)};
% [psi2RefDegInplaneLineValidSort,psi2RefDegInplaneLineValidSortI] = sort(psi2RefDegInplaneLineValid);
% % scatter(locPosLineValidX(psi2RefDegInplaneLineValidSortI),locPosLineValidY(psi2RefDegInplaneLineValidSortI),15,psi2RefDegInplaneLineValidSort,'filled','MarkerFaceAlpha',0.5)
% plot(interPointX{visLineInd(1)},interPointY{visLineInd(1)},'+','MarkerSize',5,'Color',[0.5 0.5 0.5])
% % hold on
% % psi2RefInplanePlot = psi2RefDegInplaneLine{visLineInd(1)};
% locPosLinePlotX = locPosLineValidX(psi2RefDegInplaneLineValidSortI);
% locPosLinePlotY = locPosLineValidY(psi2RefDegInplaneLineValidSortI);
% phiLinePlot = phiLineValid(psi2RefDegInplaneLineValidSortI);
% psi2RefInplanePlot = psi2RefDegInplaneLineValidSort;
% psi2RefInplanePlot(psi2RefInplanePlot > maxDegPlot) = maxDegPlot;
% psi2RefInplaneColor = round(psi2RefInplanePlot/maxDegPlot*(size(colormapArrow,1)-1))+1;
% for j = 1:length(phiLinePlot)
%     quiver(locPosLinePlotX(j),locPosLinePlotY(j),arrowScaleNoHead*cos(phiLinePlot(j)),arrowScaleNoHead*sin(phiLinePlot(j)),...
%         'AutoScale','off','LineWidth',2,'Color',colormapArrow(psi2RefInplaneColor(j),:),'ShowArrowHead','off');
% end
% axis image ij
% axis off
% % if isempty(visRot)
% visRot = rad2deg(mean(interPointSlope{visLineInd(1)}));
% % end
% camroll(visRot)
% camzoom(visZoom)
% % caxis([round(minParaMeanPsi2RefDegInplane*10)/10 round(maxParaMeanPsi2RefDegInplane*10)/10])
% % colormap(flipud(colormapC));
% % colorbar('south')
% % axis([-sizeROI/2*pixelSize*shrinkFOVX sizeROI/2*pixelSize*shrinkFOVX -sizeROI/2*pixelSize*shrinkFOVY sizeROI/2*pixelSize*shrinkFOVY])
% % title(['\psi^{inplane}_{ref} map, BB=' num2str(visLineInd(1))],'FontSize',20)
% hold off
%
% subplot(2,1,2);
% hold on
% locPosLineValidX = locPosLineX{visLineInd(2)};
% locPosLineValidY = locPosLineY{visLineInd(2)};
% phiLineValid = phiLine{visLineInd(2)};
% psi2RefDegInplaneLineValid = psi2RefDegInplaneLine{visLineInd(2)};
% [psi2RefDegInplaneLineValidSort,psi2RefDegInplaneLineValidSortI] = sort(psi2RefDegInplaneLineValid);
% % scatter(locPosLineValidX(psi2RefDegInplaneLineValidSortI),locPosLineValidY(psi2RefDegInplaneLineValidSortI),15,psi2RefDegInplaneLineValidSort,'filled','MarkerFaceAlpha',0.5)
% plot(interPointX{visLineInd(2)},interPointY{visLineInd(2)},'+','MarkerSize',5,'Color',[0.5 0.5 0.5])
% % hold on
% % psi2RefInplanePlot = psi2RefDegInplaneLine{visLineInd(1)};
% locPosLinePlotX = locPosLineValidX(psi2RefDegInplaneLineValidSortI);
% locPosLinePlotY = locPosLineValidY(psi2RefDegInplaneLineValidSortI);
% phiLinePlot = phiLineValid(psi2RefDegInplaneLineValidSortI);
% psi2RefInplanePlot = psi2RefDegInplaneLineValidSort;
% psi2RefInplanePlot(psi2RefInplanePlot > maxDegPlot) = maxDegPlot;
% psi2RefInplaneColor = round(psi2RefInplanePlot/maxDegPlot*(size(colormapArrow,1)-1))+1;
% for j = 1:length(phiLinePlot)
%     quiver(locPosLinePlotX(j),locPosLinePlotY(j),arrowScaleNoHead*cos(phiLinePlot(j)),arrowScaleNoHead*sin(phiLinePlot(j)),...
%         'AutoScale','off','LineWidth',2,'Color',colormapArrow(psi2RefInplaneColor(j),:),'ShowArrowHead','off');
% end
% axis image ij
% axis off
% % if isempty(visRot)
% visRot = rad2deg(mean(interPointSlope{visLineInd(2)}));
% % end
% camroll(visRot)
% camzoom(visZoom)
% caxis([round(minParaMeanPsi2RefDegInplane*10)/10 round(maxParaMeanPsi2RefDegInplane*10)/10])
% % colormap(flipud(colormapC));
% % colorbar('south')
% % axis([-sizeROI/2*pixelSize*shrinkFOVX sizeROI/2*pixelSize*shrinkFOVX -sizeROI/2*pixelSize*shrinkFOVY sizeROI/2*pixelSize*shrinkFOVY])
% ax = gca;
% rec = rectangle('Position',[ax.XLim(2)-scaleBar-(ax.XLim(2)-ax.XLim(1))/20,...
%     ax.YLim(2)-1-(ax.YLim(2)-ax.YLim(1))/20,...
%     scaleBar,3],'FaceColor','k','EdgeColor','k');
% % title(['\psi^{inplane}_{ref} map, BB=' num2str(visLineInd(2))],'FontSize',20)
% hold off
% h2_3_3 = gcf;
% export_fig(h2_3_3,[dirNameNewNew filesep saveNameNewNew ' h2_3_3, quiver plot, inplane phi deviation, 2 BBs' saveFormat],'-transparent')

% %% validation of assigned backbone orientation
% figure('Position',P1,'Visible','on');
% subplot(2,1,1);
% hold on
% locPosLineValidX = locPosLineX{visLineInd(1)};
% locPosLineValidY = locPosLineY{visLineInd(1)};
% locPosROIRefMuLineValid = locPosROIRefMuLine{visLineInd(1)};
% plot(interPointX{visLineInd(1)},interPointY{visLineInd(1)},'+','MarkerSize',5,'Color',[0.5 0.5 0.5])
% quiver(locPosLineValidX,locPosLineValidY,arrowScaleNoHead*locPosROIRefMuLineValid(:,1),arrowScaleNoHead*locPosROIRefMuLineValid(:,2),...
%     'AutoScale','off','LineWidth',3,'Color','r','ShowArrowHead','off');
% axis image ij
% axis off
% % if isempty(visRot)
% visRot = rad2deg(mean(interPointSlope{visLineInd(1)}));
% % end
% camroll(visRot)
% camzoom(visZoom)
% hold off
%
% subplot(2,1,2);
% hold on
% locPosLineValidX = locPosLineX{visLineInd(2)};
% locPosLineValidY = locPosLineY{visLineInd(2)};
% locPosROIRefMuLineValid = locPosROIRefMuLine{visLineInd(2)};
% plot(interPointX{visLineInd(2)},interPointY{visLineInd(2)},'+','MarkerSize',5,'Color',[0.5 0.5 0.5])
% quiver(locPosLineValidX,locPosLineValidY,arrowScaleNoHead*locPosROIRefMuLineValid(:,1),arrowScaleNoHead*locPosROIRefMuLineValid(:,2),...
%     'AutoScale','off','LineWidth',3,'Color','r','ShowArrowHead','off');
% axis image ij
% axis off
% % if isempty(visRot)
% visRot = rad2deg(mean(interPointSlope{visLineInd(2)}));
% % end
% camroll(visRot)
% camzoom(visZoom)
% ax = gca;
% rec = rectangle('Position',[ax.XLim(2)-scaleBar-(ax.XLim(2)-ax.XLim(1))/20,...
%     ax.YLim(2)-1-(ax.YLim(2)-ax.YLim(1))/20,...
%     scaleBar,3],'FaceColor','k','EdgeColor','k');
% % title(['\psi^{inplane}_{ref} map, BB=' num2str(visLineInd(2))],'FontSize',20)
% hold off
% h2_3_4 = gcf;
% export_fig(h2_3_4,[dirNameNewNew filesep saveNameNewNew ' h2_3_4, quiver plot, locPosROIRef, 2 BBs' saveFormat],'-transparent')

%% phi/inplane phi arrow plot
figure('Position',P1,'Visible','off');
hold on
quiver(locPosROI(~invalidLocROI,1),locPosROI(~invalidLocROI,2),arrowScale*cos(phiROI(~invalidLocROI)),arrowScale*sin(phiROI(~invalidLocROI)),...
    'AutoScale','off','LineWidth',2,'Color','b','ShowArrowHead','off');
quiver(lineCoorAll(1,:),lineCoorAll(2,:),arrowScale*cos(slopeLineAll),arrowScale*sin(slopeLineAll),...
    'AutoScale','off','LineWidth',3,'Color','r','ShowArrowHead','on');
grid on
axis image ij
axis off
% caxis([round(minParaMeanCosPhi*10)/10 round(maxParaMeanCosPhi*10)/10])
% colormap(colormapC);
% colorbar('south')
axis([-sizeROI/2*pixelSize*shrinkFOVX sizeROI/2*pixelSize*shrinkFOVX -sizeROI/2*pixelSize*shrinkFOVY sizeROI/2*pixelSize*shrinkFOVY])
title('\phi arrow map','FontSize',20)
hold off
h2_4 = gcf;
export_fig(h2_4,[dirNameNewNew filesep saveNameNewNew ' h2_4, arrow plot, inplane phi' saveFormat],'-transparent')

%% Omega plot
figure('Position',P1,'Visible','off');
ax1 = subplotTight(3,3,[1 4]);
imagesc(rgbPlotOmegaMean( floor(imgSize(1)/2-imgSize(1)*shrinkFOVY/2)+1:floor(imgSize(1)/2+imgSize(1)*shrinkFOVY/2),...
    floor(imgSize(2)/2-imgSize(2)*shrinkFOVX/2)+1:floor(imgSize(2)/2+imgSize(2)*shrinkFOVX/2),: ))
axis image ij
axis off
title('mean of \Omega map')
co = colorbar('southoutside');
ax5 = subplotTight(3,3,7);
imagesc(rot90(colorBarOmegaMean,2))
axis ij;
ax = gca;
ax.YTick = [ax.YLim(1) ax.YLim(2)];
ax.XTick = [ax.XLim(1) (ax.XLim(2)+ax.XLim(1))/2 ax.XLim(2)];
ax.YTickLabel= ({['\geq' num2str(indN)],0});
ax.XTickLabel = [round(minParaMeanOmega*10)/10 round((maxParaMeanOmega+minParaMeanOmega)/2*10)/10 round(maxParaMeanOmega*10)/10];
ax.Position(1) = ax.Position(1)+ax.Position(3)/8;
ax.Position(3) = ax.Position(3)*3/4;
ax.Position(2) = co.Position(2) - 4*co.Position(4);
ax.Position(4) = 2*co.Position(4);

ax2 = subplotTight(3,3,[2 5]);
imagesc(rgbPlotOmegaStd( floor(imgSize(1)/2-imgSize(1)*shrinkFOVY/2)+1:floor(imgSize(1)/2+imgSize(1)*shrinkFOVY/2),...
    floor(imgSize(2)/2-imgSize(2)*shrinkFOVX/2)+1:floor(imgSize(2)/2+imgSize(2)*shrinkFOVX/2),: ))
axis image ij
axis off
title('std \Omega map')
ax6 = subplotTight(3,3,8);
imagesc(rot90(colorBarOmegaStd,2))
axis ij;
ax = gca;
ax.YTick = [ax.YLim(1) ax.YLim(2)];
ax.XTick = [ax.XLim(1) (ax.XLim(2)+ax.XLim(1))/2 ax.XLim(2)];
ax.YTickLabel= ({['\geq' num2str(indN)],0});
ax.XTickLabel = [round(minParaStdOmega*10)/10 round((maxParaStdOmega+minParaStdOmega)/2*10)/10 round(maxParaStdOmega*10)/10];
ax.Position(1) = ax.Position(1)+ax.Position(3)/8;
ax.Position(3) = ax.Position(3)*3/4;
ax.Position(2) = co.Position(2) - 4*co.Position(4);
ax.Position(4) = 2*co.Position(4);

ax3 = subplotTight(3,3,[3 6]);
hold on
imagesc(rgbPlotAbsDevOmegaMean( floor(imgSize(1)/2-imgSize(1)*shrinkFOVY/2)+1:floor(imgSize(1)/2+imgSize(1)*shrinkFOVY/2),...
    floor(imgSize(2)/2-imgSize(2)*shrinkFOVX/2)+1:floor(imgSize(2)/2+imgSize(2)*shrinkFOVX/2),: ))
axis image ij
axis off
ax = gca;
rec = rectangle('Position',[ax.XLim(2)-scaleBar/binSize-(ax.XLim(2)-ax.XLim(1))/20,...
    ax.YLim(2)-1-(ax.YLim(2)-ax.YLim(1))/20,...
    scaleBar/binSize,1],'FaceColor','w','EdgeColor','w');
title('mean of |\Omega-mean(\Omega)| map')
hold off
ax7 = subplotTight(3,3,9);
imagesc(rot90(colorBarAbsDevOmegaMean,2))
axis ij;
ax = gca;
ax.YTick = [ax.YLim(1) ax.YLim(2)];
ax.XTick = [ax.XLim(1) (ax.XLim(2)+ax.XLim(1))/2 ax.XLim(2)];
ax.YTickLabel= ({['\geq' num2str(indN)],0});
ax.XTickLabel = [round(minParaMeanAbsDevOmega*10)/10 round((maxParaMeanAbsDevOmega+minParaMeanAbsDevOmega)/2*10)/10 round(maxParaMeanAbsDevOmega*10)/10];
ax.Position(1) = ax.Position(1)+ax.Position(3)/8;
ax.Position(3) = ax.Position(3)*3/4;
ax.Position(2) = co.Position(2) - 4*co.Position(4);
ax.Position(4) = 2*co.Position(4);

colorbar(ax1,'off');
h5 = gcf;
export_fig(h5,[dirNameNewNew filesep saveNameNewNew ' h5, rgb plot, Omega deviation' saveFormat],'-transparent')

%% omega scatter plot
figure('Position',P1,'Visible','off');
ax1 = subplotTight(2,4,[1 2 5 6]);
hold on
scatter(locPosROI(~invalidLocROI,1),locPosROI(~invalidLocROI,2),10,omegaROI(~invalidLocROI),'filled','MarkerFaceAlpha',0.3)
grid on
axis image ij
axis off
caxis([round(minParaMeanOmega*10)/10 round(maxParaMeanOmega*10)/10])
colormap(flipud(colormapC));
colorbar('south')
axis([-sizeROI/2*pixelSize*shrinkFOVX sizeROI/2*pixelSize*shrinkFOVX -sizeROI/2*pixelSize*shrinkFOVY sizeROI/2*pixelSize*shrinkFOVY])
title('\Omega map','FontSize',20)
hold off

ax2 = subplotTight(2,4,[3 4 7 8]);
hold on
scatter(locPosROI(~invalidLocROI,1),locPosROI(~invalidLocROI,2),10,abs(devOmega(~invalidLocROI)),'filled','MarkerFaceAlpha',0.3)
grid on
axis image ij
axis off
caxis([round(minParaMeanAbsDevOmega*10)/10 round(maxParaMeanAbsDevOmega*10)/10])
colormap(flipud(colormapC));
colorbar('south')
axis([-sizeROI/2*pixelSize*shrinkFOVX sizeROI/2*pixelSize*shrinkFOVX -sizeROI/2*pixelSize*shrinkFOVY sizeROI/2*pixelSize*shrinkFOVY])
ax = gca;
rec = rectangle('Position',[ax.XLim(2)-scaleBar-(ax.XLim(2)-ax.XLim(1))/20,...
    ax.YLim(2)-1-(ax.YLim(2)-ax.YLim(1))/20,...
    scaleBar,3],'FaceColor','k','EdgeColor','k');
title('|\Omega-mean(\Omega)| map','FontSize',20)
hold off
h5_1 = gcf;
export_fig(h5_1,[dirNameNewNew filesep saveNameNewNew ' h5_1, scatter plot, omega deviation' saveFormat],'-transparent')

% sorted version
omegaROIValid = omegaROI(~invalidLocROI);
absDevOmegaValid = abs(devOmega(~invalidLocROI));
[omegaROIValidSort,omegaROIValidSortI] = sort(omegaROIValid);
[absDevOmegaValidSort,absDevOmegaValidSortI] = sort(absDevOmegaValid);

figure('Position',P1);
ax1 = subplotTight(2,4,[1 2 5 6]);
hold on
scatter(locPosROIValid(omegaROIValidSortI,1),locPosROIValid(omegaROIValidSortI,2),10,omegaROIValidSort,'filled','MarkerFaceAlpha',0.3)
grid on
axis image ij
axis off
caxis([round(minParaMeanOmega*10)/10 round(maxParaMeanOmega*10)/10])
colormap(flipud(colormapC));
colorbar('south')
axis([-sizeROI/2*pixelSize*shrinkFOVX sizeROI/2*pixelSize*shrinkFOVX -sizeROI/2*pixelSize*shrinkFOVY sizeROI/2*pixelSize*shrinkFOVY])
title('\Omega map, sorted','FontSize',20)
hold off

ax2 = subplotTight(2,4,[3 4 7 8]);
hold on
scatter(locPosROIValid(absDevOmegaValidSortI,1),locPosROIValid(absDevOmegaValidSortI,2),10,absDevOmegaValidSort,'filled','MarkerFaceAlpha',0.3)
grid on
axis image ij
axis off
caxis([round(minParaMeanAbsDevOmega*10)/10 round(maxParaMeanAbsDevOmega*10)/10])
colormap(flipud(colormapC));
colorbar('south')
axis([-sizeROI/2*pixelSize*shrinkFOVX sizeROI/2*pixelSize*shrinkFOVX -sizeROI/2*pixelSize*shrinkFOVY sizeROI/2*pixelSize*shrinkFOVY])
ax = gca;
rec = rectangle('Position',[ax.XLim(2)-scaleBar-(ax.XLim(2)-ax.XLim(1))/20,...
    ax.YLim(2)-1-(ax.YLim(2)-ax.YLim(1))/20,...
    scaleBar,3],'FaceColor','k','EdgeColor','k');
title('|\Omega-mean(\Omega)| map, sorted','FontSize',20)
hold off
h5_1_1 = gcf;
export_fig(h5_1_1,[dirNameNewNew filesep saveNameNewNew ' h5_1_1, scatter plot, omega deviation, sorted' saveFormat],'-transparent')

figure('Position',P1,'Visible','on');
subplot(2,1,1);
hold on
locPosLineValidX = locPosLineX{visLineInd(1)};
locPosLineValidY = locPosLineY{visLineInd(1)};
omegaLineValid = omegaLine{visLineInd(1)};
[omegaLineValidSort,omegaLineValidSortI] = sort(omegaLineValid);
scatter(locPosLineValidX(omegaLineValidSortI),locPosLineValidY(omegaLineValidSortI),15,omegaLineValidSort,'filled','MarkerFaceAlpha',0.5)
grid on
axis image ij
axis off
% if isempty(visRot)
visRot = rad2deg(mean(interPointSlope{visLineInd(1)}));
% end
camroll(visRot)
camzoom(visZoom)
caxis([round(minParaMeanOmega*10)/10 round(maxParaMeanOmega*10)/10])
colormap(flipud(colormapC));
colorbar('south')
% axis([-sizeROI/2*pixelSize*shrinkFOVX sizeROI/2*pixelSize*shrinkFOVX -sizeROI/2*pixelSize*shrinkFOVY sizeROI/2*pixelSize*shrinkFOVY])
% title(['\psi^{inplane}_{ref} map, BB=' num2str(visLineInd(1))],'FontSize',20)
hold off

subplot(2,1,2);
hold on
locPosLineValidX = locPosLineX{visLineInd(2)};
locPosLineValidY = locPosLineY{visLineInd(2)};
omegaLineValid = omegaLine{visLineInd(2)};
[omegaLineValidSort,omegaLineValidSortI] = sort(omegaLineValid);
scatter(locPosLineValidX(omegaLineValidSortI),locPosLineValidY(omegaLineValidSortI),15,omegaLineValidSort,'filled','MarkerFaceAlpha',0.5)
grid on
axis image ij
axis off
% if isempty(visRot)
visRot = rad2deg(mean(interPointSlope{visLineInd(2)}));
% end
camroll(visRot)
camzoom(visZoom)
caxis([round(minParaMeanOmega*10)/10 round(maxParaMeanOmega*10)/10])
colormap(flipud(colormapC));
colorbar('south')
% axis([-sizeROI/2*pixelSize*shrinkFOVX sizeROI/2*pixelSize*shrinkFOVX -sizeROI/2*pixelSize*shrinkFOVY sizeROI/2*pixelSize*shrinkFOVY])
ax = gca;
rec = rectangle('Position',[ax.XLim(2)-scaleBar-(ax.XLim(2)-ax.XLim(1))/20,...
    ax.YLim(2)-1-(ax.YLim(2)-ax.YLim(1))/20,...
    scaleBar,3],'FaceColor','k','EdgeColor','k');
% title(['\psi^{inplane}_{ref} map, BB=' num2str(visLineInd(2))],'FontSize',20)
hold off
h5_1_2 = gcf;
export_fig(h5_1_2,[dirNameNewNew filesep saveNameNewNew ' h5_1_2, scatter plot, omega deviation, 2 BBs' saveFormat],'-transparent')

%% tau_on scatter plot
figure('Position',P1,'Visible','off');
ax1 = subplotTight(2,6,[1 2 7 8]);
hold on
scatter(locPosGroup(:,1),locPosGroup(:,2),10,onTimeCorrAll,'filled','MarkerFaceAlpha',0.3)
grid on
axis image ij
axis off
caxis([round(minParaMeanTau*1000)/1000 round(maxParaMeanTau*1000)/1000])
colormap(flipud(colormapC));
colorbar('south')
axis([-sizeROI/2*pixelSize*shrinkFOVX sizeROI/2*pixelSize*shrinkFOVX -sizeROI/2*pixelSize*shrinkFOVY sizeROI/2*pixelSize*shrinkFOVY])
title('\tau_{on} map','FontSize',20)
hold off

ax2 = subplotTight(2,6,[3 4 9 10]);
hold on
scatter(locPosGroup(:,1),locPosGroup(:,2),10,photonCorrAll,'filled','MarkerFaceAlpha',0.3)
grid on
axis image ij
axis off
caxis([photonLowThreshold median(photonCorrAll)*3])
colormap(flipud(colormapC));
colorbar('south')
axis([-sizeROI/2*pixelSize*shrinkFOVX sizeROI/2*pixelSize*shrinkFOVX -sizeROI/2*pixelSize*shrinkFOVY sizeROI/2*pixelSize*shrinkFOVY])
ax = gca;
rec = rectangle('Position',[ax.XLim(2)-scaleBar-(ax.XLim(2)-ax.XLim(1))/20,...
    ax.YLim(2)-1-(ax.YLim(2)-ax.YLim(1))/20,...
    scaleBar,3],'FaceColor','k','EdgeColor','k');
title('photon per burst map','FontSize',20)
hold off

ax3 = subplotTight(2,6,[5 6 11 12]);
hold on
scatter(locPosGroup(:,1),locPosGroup(:,2),10,brightnessHzAll,'filled','MarkerFaceAlpha',0.3)
grid on
axis image ij
axis off
caxis([round(minParaMeanHz) round(maxParaMeanHz)])
colormap(flipud(colormapC));
colorbar('south')
axis([-sizeROI/2*pixelSize*shrinkFOVX sizeROI/2*pixelSize*shrinkFOVX -sizeROI/2*pixelSize*shrinkFOVY sizeROI/2*pixelSize*shrinkFOVY])
ax = gca;
rec = rectangle('Position',[ax.XLim(2)-scaleBar-(ax.XLim(2)-ax.XLim(1))/20,...
    ax.YLim(2)-1-(ax.YLim(2)-ax.YLim(1))/20,...
    scaleBar,3],'FaceColor','k','EdgeColor','k');
title('brightness/\tau_{on} map [kHz]','FontSize',20)
hold off
h6 = gcf;
export_fig(h6,[dirNameNewNew filesep saveNameNewNew ' h6, scatter plot, tau and brightness' saveFormat],'-transparent')

%% Omega vs psi
figure('Position',P1,'Visible','on');
ax1 = subplot(2,4,[1 2 5 6]);
hold on
scatter(psi2RefDegValid,omegaROIValid,10,'filled','MarkerFaceAlpha',0.3)
grid on
% axis image ij
title('\Omega vs \psi_{ref}')
hold off
ax2 = subplot(2,4,[3 4 7 8]);
hold on
scatter(psi2RefDegInplaneValid,omegaROIValid,10,'filled','MarkerFaceAlpha',0.3)
grid on
% axis image ij
title('\Omega vs \psi^{inplane}_{ref}')
hold off
h7 = gcf;
export_fig(h7,[dirNameNewNew filesep saveNameNewNew ' h7, scatter plot, omega vs ps in ROI' saveFormat],'-transparent')

%% Segmented plots
% trajectories of stds of psi2RefDegInplane and omega on the specified
% backbones
figure('Position',P1,'Visible','on');
ax1 = subplot(2,2,[1]);
xx = [0:length(psi2RefDegInplaneStdSeg{visLineInd(1)})-1]*interPDist;
plot(xx,psi2RefDegInplaneStdSeg{visLineInd(1)})
xlim([xx(1) xx(end)])
% legend(['med \psi^{inplane}_{ref}=' num2str(round(psiDegInplaneLineMed(visLineInd(1)),3,'significant'))])
% xlabel('distance')
set(gca, 'XTickLabel', []);
ylabel('std [deg]')
title(['std of \psi^{inplane}_{ref}, backbone ' num2str(visLineInd(1)) ', ' num2str(interPDist*2) 'nm window'])

ax2 = subplot(2,2,[3]);
xx = [0:length(omegaROIMedSeg{visLineInd(1)})-1]*interPDist;
plot(xx,omegaROIMedSeg{visLineInd(1)})
xlim([xx(1) xx(end)])
% legend(['med \Omega=' num2str(round(omegaLineMed(visLineInd(1)),3,'significant'))])
xlabel('distance')
ylabel('Median [sr]')
title(['Median of \Omega, backbone ' num2str(visLineInd(1)) ', ' num2str(interPDist*2) 'nm window'])

ax3 = subplot(2,2,[2]);
xx = [0:length(psi2RefDegInplaneStdSeg{visLineInd(2)})-1]*interPDist;
plot(xx,psi2RefDegInplaneStdSeg{visLineInd(2)})
xlim([xx(1) xx(end)])
% legend(['med \psi^{inplane}_{ref}=' num2str(round(psiDegInplaneLineMed(visLineInd(2)),3,'significant'))])
% xlabel('distance')
set(gca, 'XTickLabel', []);
ylabel('std [deg]')
title(['std of \psi^{inplane}_{ref}, backbone ' num2str(visLineInd(2)) ', ' num2str(interPDist*2) 'nm window'])

ax4 = subplot(2,2,[4]);
xx = [0:length(omegaROIMedSeg{visLineInd(2)})-1]*interPDist;
plot(xx,omegaROIMedSeg{visLineInd(2)})
xlim([xx(1) xx(end)])
% legend(['med \Omega=' num2str(round(omegaLineMed(visLineInd(2)),3,'significant'))])
xlabel('distance')
ylabel('Median [sr]')
title(['Median of \Omega, backbone ' num2str(visLineInd(2)) ', ' num2str(interPDist*2) 'nm window'])
h8 = gcf;
export_fig(h8,[dirNameNewNew filesep saveNameNewNew ' h8, trajectory of std of psi2RefInplane vs med omega on BBs' saveFormat],'-transparent')

%% histogram of psi ref inplane and omega in two backbone regions
figure('Position',P1,'Visible','on');
ax1 = subplot(1,2,[1]);
hold on
hist1 = histogram(psi2RefDegInplaneLine{visLineInd(1)},'Normalization','probability');
hist2 = histogram(psi2RefDegInplaneLine{visLineInd(2)},'Normalization','probability');
binEdges = hist1.BinEdges;
hist2.BinEdges = binEdges;
legend(['BB' num2str(visLineInd(1)) ',med=' num2str(round(psi2RefDegInplaneLineMed(visLineInd(1)),3,'significant'))],...
    ['BB' num2str(visLineInd(2)) ',med=' num2str(round(psi2RefDegInplaneLineMed(visLineInd(2)),3,'significant'))])
xlabel('\psi^{inplane}_{ref} [deg]')
% set(gca, 'XTickLabel', []);
% ylabel('std [deg]')
title(['Inplane deviation'])
hold off

ax2 = subplot(1,2,[2]);
hold on
hist1 = histogram(omegaLine{visLineInd(1)},'Normalization','probability');
hist2 = histogram(omegaLine{visLineInd(2)},'Normalization','probability');
binEdges = hist1.BinEdges;
hist2.BinEdges = binEdges;
legend(['BB' num2str(visLineInd(1)) ',med=' num2str(round(omegaLineMed(visLineInd(1)),3,'significant'))],...
    ['BB' num2str(visLineInd(2)) ',med=' num2str(round(omegaLineMed(visLineInd(2)),3,'significant'))])
xlabel('\Omega [sr]')
title(['Wobbling area'])
hold off

h8_1 = gcf;
export_fig(h8_1,[dirNameNewNew filesep saveNameNewNew ' h8_1, hist, psi2RefInplane vs med omega on BBs' saveFormat],'-transparent')

%% Visualizaton of psi2RefDegInplane and omega distribution within the region
% with the largest psi2RefDegInplaneStdSeg on the longest backbone
figure('Position',P1,'Visible','on');
ax1 = subplot(4,4,[1 2 5 6]);
plot(interPointX{visLineInd(1)},interPointY{visLineInd(1)},'+','MarkerSize',5,'Color',[0.5 0.5 0.5])
hold on
maxPsi2RefInplanePlot = maxPsi2RefInplane_psi2RefInplane{visLineInd(1)};
maxPsi2RefInplanePlot(maxPsi2RefInplanePlot > maxDegPlot) = maxDegPlot;
psi2RefInplaneColor = round(maxPsi2RefInplanePlot/maxDegPlot*(size(colormapArrow,1)-1))+1;
for j = 1:length(maxPsi2RefInplane_phi{visLineInd(1)})
    %     quiver(maxPsi2RefInplane_locPos{visLineInd(1)}(j,1),maxPsi2RefInplane_locPos{visLineInd(1)}(j,2),...
    %         arrowScaleNoHead*cos(maxPsi2RefInplane_phi{visLineInd(1)}(j)),arrowScaleNoHead*sin(maxPsi2RefInplane_phi{visLineInd(1)}(j)),...
    %         'AutoScale','off','LineWidth',2,'Color',colormapArrow(psi2RefInplaneColor(j),:),'ShowArrowHead','off');
    line([maxPsi2RefInplane_locPos{visLineInd(1)}(j,1)-arrowScaleNoHead/2*cos(maxPsi2RefInplane_phi{visLineInd(1)}(j))...
        maxPsi2RefInplane_locPos{visLineInd(1)}(j,1)+arrowScaleNoHead/2*cos(maxPsi2RefInplane_phi{visLineInd(1)}(j))],...
        [maxPsi2RefInplane_locPos{visLineInd(1)}(j,2)-arrowScaleNoHead/2*sin(maxPsi2RefInplane_phi{visLineInd(1)}(j))...
        maxPsi2RefInplane_locPos{visLineInd(1)}(j,2)+arrowScaleNoHead/2*sin(maxPsi2RefInplane_phi{visLineInd(1)}(j))],...
        'LineWidth',2,'Color',colormapArrow(psi2RefInplaneColor(j),:));
end
xlim([maxPsi2RefInplane_interPoint{visLineInd(1)}(1)-visMargin maxPsi2RefInplane_interPoint{visLineInd(1)}(1)+visMargin])
ylim([maxPsi2RefInplane_interPoint{visLineInd(1)}(2)-visMargin maxPsi2RefInplane_interPoint{visLineInd(1)}(2)+visMargin])
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'Box', 'on');
axis square ij
title('\psi^{inplane}_{ref} in the largest \psi^{inplane}_{ref} std region')
hold off
ax2 = subplot(4,4,[3 4 7 8]);
plot(interPointX{visLineInd(1)},interPointY{visLineInd(1)},'+','MarkerSize',5,'Color',[0.5 0.5 0.5])
hold on
minPsi2RefInplanePlot = minPsi2RefInplane_psi2RefInplane{visLineInd(1)};
minPsi2RefInplanePlot(minPsi2RefInplanePlot > maxDegPlot) = maxDegPlot;
psi2RefInplaneColor = round(minPsi2RefInplanePlot/maxDegPlot*(size(colormapArrow,1)-1))+1;
for j = 1:length(minPsi2RefInplane_phi{visLineInd(1)})
    %     quiver(minPsi2RefInplane_locPos{visLineInd(1)}(j,1),minPsi2RefInplane_locPos{visLineInd(1)}(j,2),...
    %         arrowScaleNoHead*cos(minPsi2RefInplane_phi{visLineInd(1)}(j)),arrowScaleNoHead*sin(minPsi2RefInplane_phi{visLineInd(1)}(j)),...
    %         'AutoScale','off','LineWidth',2,'Color',colormapArrow(psi2RefInplaneColor(j),:),'ShowArrowHead','off');
    line([minPsi2RefInplane_locPos{visLineInd(1)}(j,1)-arrowScaleNoHead/2*cos(minPsi2RefInplane_phi{visLineInd(1)}(j))...
        minPsi2RefInplane_locPos{visLineInd(1)}(j,1)+arrowScaleNoHead/2*cos(minPsi2RefInplane_phi{visLineInd(1)}(j))],...
        [minPsi2RefInplane_locPos{visLineInd(1)}(j,2)-arrowScaleNoHead/2*sin(minPsi2RefInplane_phi{visLineInd(1)}(j))...
        minPsi2RefInplane_locPos{visLineInd(1)}(j,2)+arrowScaleNoHead/2*sin(minPsi2RefInplane_phi{visLineInd(1)}(j))],...
        'LineWidth',2,'Color',colormapArrow(psi2RefInplaneColor(j),:));
end
xlim([minPsi2RefInplane_interPoint{visLineInd(1)}(1)-visMargin minPsi2RefInplane_interPoint{visLineInd(1)}(1)+visMargin])
ylim([minPsi2RefInplane_interPoint{visLineInd(1)}(2)-visMargin minPsi2RefInplane_interPoint{visLineInd(1)}(2)+visMargin])
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'Box', 'on');
axis square ij
title('\psi^{inplane}_{ref} in the smallest \psi^{inplane}_{ref} std region')
hold off
ax3 = subplot(4,4,[9 10 13 14]);
plot(interPointX{visLineInd(1)},interPointY{visLineInd(1)},'+','MarkerSize',5,'Color',[0.5 0.5 0.5])
hold on
scatter(maxOmega_locPos{visLineInd(1)}(:,1),maxOmega_locPos{visLineInd(1)}(:,2),20,maxOmega_omega{visLineInd(1)},'filled')
caxis([round(minParaMeanOmega*10)/10 round(maxParaMeanOmega*10)/10])
colormap(flipud(colormapC));
xlim([maxOmega_interPoint{visLineInd(1)}(1)-visMargin maxOmega_interPoint{visLineInd(1)}(1)+visMargin])
ylim([maxOmega_interPoint{visLineInd(1)}(2)-visMargin maxOmega_interPoint{visLineInd(1)}(2)+visMargin])
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'Box', 'on');
axis square ij
title('\Omega in the largest \Omega med region')
hold off
ax1 = subplot(4,4,[11 12 15 16]);
plot(interPointX{visLineInd(1)},interPointY{visLineInd(1)},'+','MarkerSize',5,'Color',[0.5 0.5 0.5])
hold on
scatter(minOmega_locPos{visLineInd(1)}(:,1),minOmega_locPos{visLineInd(1)}(:,2),20,minOmega_omega{visLineInd(1)},'filled')
caxis([round(minParaMeanOmega*10)/10 round(maxParaMeanOmega*10)/10])
colormap(flipud(colormapC));
xlim([minOmega_interPoint{visLineInd(1)}(1)-visMargin minOmega_interPoint{visLineInd(1)}(1)+visMargin])
ylim([minOmega_interPoint{visLineInd(1)}(2)-visMargin minOmega_interPoint{visLineInd(1)}(2)+visMargin])
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'Box', 'on');
axis square ij
title('\Omega in the smallest \Omega med region')
hold off
h9 = gcf;
export_fig(h9,[dirNameNewNew filesep saveNameNewNew ' h9, trajectory of std of psi2RefInplane vs omega on longest BB' saveFormat],'-transparent')

clear h* Bx By allVectorDis brightness_* nearNeiInd rgbPlot*

save([dirNameNewNew filesep saveNameNewNew ' all analyzed results' '.mat'])

frmAll = [frmROI; frmOFF];
locPosAll = [locPosROI; locPosOFF];
sigAll = [sigROI; sigOFF];
MAll = [MROI; MOFF];
thetaAll = [thetaROI;thetaOFF];
phiAll = [phiROI;phiOFF];
omegaAll = [omegaROI;omegaOFF];

loc_data_ana = [frmAll locPosAll sigAll MAll thetaAll phiAll omegaAll/2];
tau_on_fibrils = [posCorrAll onTimeCorrAll];

save([dirNameNewNew filesep saveNameNewNew 'loc and tau_on data.mat'],'loc_data_ana','tau_on_fibrils')