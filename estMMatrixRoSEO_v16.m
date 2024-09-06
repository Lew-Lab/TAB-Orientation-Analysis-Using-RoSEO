% 210403 - v16 TD - Remove localizations too close to FOV edges.

% 210128 - TD - Updated PSF crop format, that is locDXYCha removal
% depending on locDInd. Also, use old RoSEO instead of RoSEO_rcnstImg.

% 201213 - TD - Update zernikecoefficient input style. Read
% zernikecoefficient mat file instead of analyzed-result mat file.

% 201210 - v14_1 TD - Use RoSEO_rcnstImg to reconstruct noise less image.
% Save the reconstructed and raw images.

% 201210 - v14 TD - Add estimated zernikecoefficient

% 200409 - TD - Changed crop pixel center and added email sending for work
% station analysis

% 200319 - v13 TD - Read ThunderSTORM histogram instead of csv file for
% defining working FOV quickly. Change tiffread2 to Tiff class within
% tformGenSM. Reading an image sequence with Tiff is much faster than
% tiffread2. See resutls in \Detector calibration, offset subtraction, data
% transfer\Analysis\200319_readImgSpeed_v1.

% 200316 - v12 TD - Turn back to refractiveIndxSam for classical Nanoscope.
% Use reconstruction instead of scatter plots for defining FOV. Introduce
% parameters representing background estimator prefix and numApt.

% 190718 - TD - Changed refractiveIndxSam to sampleRefractiveIndx for new
% Nanoscope

% 190607,09,13 - TD - Changed h2 visualization. Changed crop region for
% h1_1. Test pixel shift. 

% 190604 - v11 TD - Use RoSEO instead of vIniCrossTZero. Do not start
% estimation from zeros for the cross terms. Use original SVD. Remove
% projection and visualization parts (moved to readRidgeDetCalDisorder_v6.m
% or later)

% 190512 - v10 TD - Performe secondM2SymmCone projection every anaLoop 
% frame group. Obtain average localized PSF in both channel.

% 190507 - v9 TD - heterogeneity analysis using inner product and
% gamma/alpha deviation. Turn off PCA.

% 190501 - v8 TD - Analyze second moments directly using PCA.

% 190412 - v7 TD - Added FOV checking loop and fibril ROI selection. Also
% included tform generation/identification function

% 190411 - v6 TD - Changed the saving format for cropped images (.mat ->
% .h5)

% 190411 - TD - Added colorbar in the h1 reconstruction and another
% visualization for estimated signals (h6). Added 'phasemaskpara' as an
% input for effectively change phase mask. Added 'refractiveIndxSam' as an
% input parameter and put this and 'emitWaveL' as input to 'Nanoscope'.

% 190404 - v5 TD - use FI to weight projection from M vector to angles.

% 190319 - v4_2 TD - remove weighted FPSFs and use RoSEO_vIniCrossTZero for
% start estimation of cross terms (mu_xy,mu_xz,mu_yz) from zeros instead of
% SVD results.

% 190301 - v4_1 TD - test weighted FPSFs

% 190225 - v4 TD - Full estimation including mux, muy, rotational
% constraint. And some code tuning for compatibility with workstation
% (ll03) analysis

% 190114 - v3 TD - Calculate background using
% colRawGaussFit_waveletFit_backgroundCorr_v1

% 190114 - v2 TD - Calculate uniform background for RoSEO estimation

% 190111 Tianben Ding
% Estimate M matrix of single molecules using RoSE-O. See startHere.mlx for
% more details.

% Input raw data and estimated background data should be in tif sequence.
% Offset data should be in multi-tif format...

%% Analysis configuration
clear; close all;clc;
tic;

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
%% Input parameters - fixed parameters, you do not need to change these parameters
pixelSize = 58.5; % nm
conversionF = 0.29; %0.49;% conversion factor (electron/count)

sampleRefractiveIndx=1.334;%1.518; % refractive index of the sample, water: 1.334, matched: 1.518, air: 1
% sampleRefractiveIndx=1.334;
numApt =  1.4; % 1.5; %numeric aperture


binSize = 20; % nm, reconstruction bin size
satPop = 0.99; % 0-1, saturation rate in the final reconstruction

% only for two channel registration using data sets captured with the standard PSF
regPara.boundaryXPosSM = 1024;%1034;%1060; % boudary pixel of two channels in x axis for single molecule imaging (pixel)
regPara.imageSizeVSM = 200; % vertical image size(pixel)
regPara.imageSizeHSM = 2048; % horizontal image size(pixel)
regPara.sigmaLowThreshold = 50;% threshold for removing localization results with very small sigma (nm)
regPara.sigmaHighThreshold = 150;% threshold for removing localization results with very large sigma (nm)
regPara.photonLowThresholdReg = 0; % photon threshold for two channel registration
regPara.squareS = 7; % square size for intensity integration (pixel, odd)
regPara.preThre = binSize/2; % precision thresholding, only keep localizations with very high precision
regPara.bmin=0;% nm, 'BinLimits' in a visualization of tformGenSM_v1 or later
regPara.bmax=800; %nm, 'BinLimits' in a visualization of tformGenSM_v1 or later
regPara.biasThre = 15;% bias tolerance
regPara.ensFactor = 2;% factor of most common ensemble

photonLowThreshold = 100; % threshold for removing localization results with very dim intensity (photon number)

colorMap = parula_hdr;

scaleBar = 300; % nm, scalebar size in the SR reconstruction

threROI = 0;

saveFormat = '.png';
%% variable parameters depending on data sets
workStatInd = 0; % index for identify if the code is working on workstation or not, 0: local computer, 1: workstation

sizeROI = 177;% pixels, odd number
centerROI =  [593    98];%center pixel coordinate of the ROI in the left channel
shrinkFOVX = 0.98; % shrink rate of final FOV visualization
shrinkFOVY = 0.8;

emitWaveL =  676; %nm, emission wave length, center at BP
% emitWaveL = 593;
% channelRatio = 0.8873; % for ThT
channelRatio = 1.1400; % for NB, based on Jin's lipid data
% channelRatio = 1.1449; % intensity ratio y cha/x cha, from Jin's NR on lipid individual experiment (NR perpendicular to coverslip) with Di03-488/561, 523/610 bandpass
% channelRatio = 1.3118;%1.1398;%1.5350; %2.0741;% Di03-R488/561, 593/46
% channelRatio = 1.1613; % intensity ratio y cha/x cha, from alignment checking slide with Di03-488/561, 523/610 bandpass
alpha_mean= [0 0];%[0.17 0.15];%[x_channel y_channel] or [right channel left channel] zeroth order leakage factor; set to 0 if no zeroth order leakage

frameN = 10000; % total frame number to be analyzed

% penalty of localization. In Jin's NR on lipid experiment, 0.2 for non
% diffusive bright molecules, 0.5 for diffusive dimmer molecules. In
% Hesam's RoSE (not RoSEO) analysis, 0.21 for 190109 Data9 amyloid data
% (standard PSF) 
reg = 0.25;%0.5;%

dataDate = '190729';
standardPSFData = 'Data12'; % standard PSF data name, should be the same FOV as your orientaton data
oriPSFData = 'Data12'; % orientation PSF data name
offsetData = 'Data18_stack'; % offset data name, tif file in the folder should tif-stack instead of tif-sequence

avePSFSize = 11;% crop pixel size for PSF averaging

% mask information, parameters of the phase mask mounted on SLM, the
% structure names should be in the lower case due to opt2struct function
phasemaskpara=struct('maskname','flatMask_0_PupilRadius_80',...
    'pupilradius',80,...%84.375,...%
    'x_shift_center_position',0,...
    'y_shift_center_position',0,...
    'maskrotation',0);
% mask options:
% tri-spot
% mBisectedXY, pRamp_3pi, pRadius_80
% mBisectedYYZZ_matched
% flatMask_0_PupilRadius_80
% trispot_3pi_r80_17-Oct-2018_aperture_80
% radial
% radialNoAperN90deg_abrComp_xyCha_201123Aper
% radialNoAperN90deg_abrComp_xyCha_210115R635,676_37Aper

remPix = 5; % margin for removing localizations too close to FOV edges

zernikeFileFolder = []; % folder containing retrieved Zernike polynomial, leave it empty if you do not consider aberration

% geometric transform information for channel registration (rough
% registration by bright nanostructures)
tformDataBead = ['190110_controlPointsMatching_r5' filesep '190110 Red beads,Di03-R488-561,523_610,vS200 geometric transformation data.mat'];

%-- use anaLoop = 1000 for workstation
anaLoop = 200; % partitioned frame number for preventing RAM overflow
%--

%-- save raw frame with localized positions, no meaning if working on
% workstation
saveSFrame = 1;
saveFrame = 50;
%--

%-- change format for workstation analysis
% data path
dataFolderPath = '/home/tding/Raw data/190729_2158BBHPTara laser irradiation tests';% for linux workstation
% dataFolderPath = 'G:\RawData\190729_2158BBHPTara laser irradiation tests'; % for local windows computer
BGEstSuf = 'lineFitWavelet';
%--

saveName = [dataDate ' ' oriPSFData ' reg' num2str(reg) ' ' BGEstSuf];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get information for saving
t=datetime('today');

dirName=datestr(t,'yymmdd');
dirName=[dirName '_' mfilename];

if exist(dirName,'dir') ~= 7
    mkdir(dirName);
end

save([dirName filesep saveName ' initial setting' '.mat'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Working ROI and two channel registration
% suppress a warning message due to unknown field provided by
% HCImageLive software
warningID = 'imageio:tiffmexutils:libtiffWarning';
warning('off',warningID)
warningID = 'MATLAB:imagesci:tiffmexutils:libtiffWarning';
warning('off',warningID)

% generate or load registration function constructed by standard PSF
tformFolder = [folderPath filesep 'tform'];
tformFile = [dataDate ' ' standardPSFData ' tform.mat'];
tformInfo = dir([tformFolder filesep tformFile]);
if ~isempty(tformInfo)
    load([tformFolder filesep tformFile])
    
%     % Read localized data
%     num = csvread([dataFolderPath filesep 'offsetSubtracted' filesep standardPSFData '.csv'],1,0);
%     num(:,5:end) = [];
%     
%     % Remove very small sigma localization
%     index = num(:,4) < regPara.sigmaLowThreshold;
%     num(index,:)=[];
%     
%     % Remove very large sigma localization
%     index = num(:,4) > regPara.sigmaHighThreshold;
%     num(index,:)=[];
        
    t = Tiff([dataFolderPath filesep 'offsetSubtracted' filesep standardPSFData '_Histogram.tif']);
    recHist = t.read();
    t.close();
    % define working region and shrink rate for final visulaization
    [centerROI,sizeROI,shrinkFOVX,shrinkFOVY,fig] =...
        checkFOV_v3(centerROI,sizeROI,shrinkFOVX,shrinkFOVY,recHist,pixelSize,regPara,...
        tformSM_L2R,P1,colorMap,satPop,scaleBar);
    
    export_fig(fig.h0,[dirName filesep saveName ' h0, scatter, localizations in both cha, large FOV' saveFormat],'-transparent')
    export_fig(fig.h0_1,[dirName filesep saveName ' h0_1, scatter, localizations in both cha, zoom' saveFormat],'-transparent')
    close(fig.h0,fig.h0_1)
    clear fig       
else
    % Load geometric transform information
    load(tformDataBead)
       
    t = Tiff([dataFolderPath filesep 'offsetSubtracted' filesep standardPSFData '_Histogram.tif']);
    recHist = t.read();
    t.close();
    
    % define working region and shrink rate for final visulaization
    [centerROI,sizeROI,shrinkFOVX,shrinkFOVY,fig] =...
        checkFOV_v3(centerROI,sizeROI,shrinkFOVX,shrinkFOVY,recHist,pixelSize,regPara,...
        tform,P1,colorMap,satPop,scaleBar);
    
    export_fig(fig.h0,[dirName filesep saveName ' h0, scatter, localizations in both cha, large FOV' saveFormat],'-transparent')
    export_fig(fig.h0_1,[dirName filesep saveName ' h0_1, scatter, localizations in both cha, zoom' saveFormat],'-transparent')
    close(fig.h0,fig.h0_1)
    clear fig
    
    centerROIFlip = centerROI;
    centerROIFlip(1) = -centerROI(1) + regPara.boundaryXPosSM;
    
    zoomRegion = [centerROIFlip(1)-(sizeROI-1)/2 centerROIFlip(2)-(sizeROI-1)/2;...
        centerROIFlip(1)+(sizeROI-1)/2 centerROIFlip(2)+(sizeROI-1)/2];
    zoomRegionRight = transformPointsInverse(tform,zoomRegion.*pixelSize);
    zoomRegionRight = [floor(zoomRegionRight(1,1)/pixelSize)+1 floor(zoomRegionRight(2,1)/pixelSize)+1 floor(zoomRegionRight(1,2)/pixelSize)+1 floor(zoomRegionRight(2,2)/pixelSize)+1];
    zoomRegion = [zoomRegion(1,1) zoomRegion(2,1) zoomRegion(1,2) zoomRegion(2,2)];
    
    % Read localized data
    num = csvread([dataFolderPath filesep 'offsetSubtracted' filesep standardPSFData '.csv'],1,0);
    
    % Remove very small sigma localization
    index = num(:,4) < regPara.sigmaLowThreshold;
    num(index,:)=[];
    
    % Remove very large sigma localization
    index = num(:,4) > regPara.sigmaHighThreshold;
    num(index,:)=[];
    
    % Clear the localized emitters which are not in the working FOV (right channel)
    indexXL = num(:,2) > (zoomRegionRight(1) + regPara.boundaryXPosSM)*pixelSize+regPara.squareS*pixelSize;
    indexXU = num(:,2) < (zoomRegionRight(2) + regPara.boundaryXPosSM)*pixelSize-regPara.squareS*pixelSize;
    indexYL = num(:,3) > zoomRegionRight(3)*pixelSize+regPara.squareS*pixelSize;
    indexYU = num(:,3) < zoomRegionRight(4)*pixelSize-regPara.squareS*pixelSize;
    indexAll = indexXL + indexXU + indexYL + indexYU;
    
    numRight=num(indexAll > 3,:);
    
    % Clear the localized emitters which are not in the working FOV (left channel)
    indexXL = num(:,2) > (-zoomRegion(2) + regPara.boundaryXPosSM)*pixelSize+regPara.squareS*pixelSize;
    indexXU = num(:,2) < (-zoomRegion(1) + regPara.boundaryXPosSM)*pixelSize-regPara.squareS*pixelSize;
    indexYL = num(:,3) > zoomRegion(3)*pixelSize+regPara.squareS*pixelSize;
    indexYU = num(:,3) < zoomRegion(4)*pixelSize-regPara.squareS*pixelSize;
    indexAll = indexXL + indexXU + indexYL + indexYU;
    
    numLeft=num(indexAll > 3,:);
    
    % Only keep the emitters in the working FOV
    num = [numLeft; numRight];
    
    %-- generate tform function using localization of standard PSF
    if workStatInd == 0
        [tformSM_L2R,tformSM_R2L,fig] =...
            tformGenSM_v2(num,pixelSize,conversionF,channelRatio,tform,tformRight2Left,zoomRegion*pixelSize,...
            dataFolderPath,standardPSFData,binSize,regPara,P1,1);
        
        export_fig(fig.h1,[dirName filesep saveName ' h0_2, scatter, localizations with high prec, without pair and corr' saveFormat],'-transparent')
        export_fig(fig.h2,[dirName filesep saveName ' h0_3, hist, prec distribution after thresholding' saveFormat],'-transparent')
        export_fig(fig.h3,[dirName filesep saveName ' h0_4, hist, all possible dist after prec thresholding before pairing' saveFormat],'-transparent')
        
        if isfield(fig,'h4')
            export_fig(fig.h4,[dirName filesep saveName ' h0_5, hist, all possible CosSin after prec, dist thresholding' saveFormat],'-transparent')
            export_fig(fig.h5,[dirName filesep saveName ' h0_6, scatter, localization after pairing, without corr' saveFormat],'-transparent')
            export_fig(fig.h6,[dirName filesep saveName ' h0_7, hist, prec of localization for pairing' saveFormat],'-transparent')
            export_fig(fig.h7,[dirName filesep saveName ' h0_8, hist, x,y bias for paired localizations, without corr' saveFormat],'-transparent')
            export_fig(fig.h8,[dirName filesep saveName ' h0_9, hist, dist between paired localizations, without corr' saveFormat],'-transparent')
            
            export_fig(fig.h9,[dirName filesep saveName ' h0_10, scatter, localization after pairing, with corr' saveFormat],'-transparent')
            export_fig(fig.h10,[dirName filesep saveName ' h0_11, hist, x,y bias for paired localizations, with corr' saveFormat],'-transparent')
            export_fig(fig.h11,[dirName filesep saveName ' h0_12, hist, dist between paired localizations, with corr' saveFormat],'-transparent')
            close(fig.h4,fig.h5,fig.h6,fig.h7,fig.h8,fig.h9,fig.h10,fig.h11)
        end
        close(fig.h1,fig.h2,fig.h3)
        clear fig
    elseif workStatInd == 1
        [tformSM_L2R,tformSM_R2L,~] =...
            tformGenSM_v2(num,pixelSize,conversionF,channelRatio,tform,tformRight2Left,zoomRegion*pixelSize,...
            dataFolderPath,standardPSFData,binSize,regPara,P1,0);
    end
    %--
    
    save([tformFolder filesep tformFile],'tformSM_L2R','tformSM_R2L')
end

clear num

%% Construct or read raw and background image data
% option of raw image input for ROSE_O
optionsImg={'TransformScale',pixelSize,'registershifty',0,'registershiftx',0,...
    'nointeractive',true,...
    'sizeROI',sizeROI,...
    'centerroi',centerROI,...
    'datapath',[dataFolderPath filesep oriPSFData],...
    'offsetpath',[dataFolderPath filesep offsetData],...
    'tform',tformSM_L2R};

% option for background image imput for ROSE_O
optionsBack={'TransformScale',pixelSize,'registershifty',0,'registershiftx',0,...
    'nointeractive',true,'nooffset',true,...
    'sizeROI',sizeROI,...
    'centerroi',centerROI,...
    'datapath',[dataFolderPath filesep 'offsetSubtracted' filesep oriPSFData filesep 'BGEst_' BGEstSuf],...%'BGEst_lfw'],...
    'offsetpath',[dataFolderPath filesep offsetData],...
    'tform',tformSM_L2R};
%%
n_default=Nanoscope('pixelSize',pixelSize,'imageSize',sizeROI,'emissWavelength',emitWaveL,...
    'sampleRefractiveIndx',sampleRefractiveIndx,'phasemaskpara',phasemaskpara,'ADcount',conversionF,'numApt',numApt);

optionsImgStruct=opt2struct(optionsImg);
optionsBackStruct=opt2struct(optionsBack);

if exist([optionsImgStruct.datapath filesep '4RoSEO'],'dir') ~= 7
    mkdir([optionsImgStruct.datapath filesep '4RoSEO']);
end
imgMatInfo = dir([optionsImgStruct.datapath filesep '4RoSEO' filesep '*' '.mat']);
if exist([optionsBackStruct.datapath filesep '4RoSEO'],'dir') ~= 7
    mkdir([optionsBackStruct.datapath filesep '4RoSEO']);
end
backMatInfo = dir([optionsBackStruct.datapath filesep '4RoSEO' filesep '*' '.mat']);

if isempty(imgMatInfo)
    [SMLM_img,~,~,~]=n_default.loadImgInteractive(optionsImg{:});
    SMLM_img = reshape(SMLM_img,[size(SMLM_img,1),size(SMLM_img,2),size(SMLM_img,4)]);
    save([optionsImgStruct.datapath filesep '4RoSEO' filesep 'croppedOffsetSubtractedData4RoSEO.mat'],'n_default','optionsImg','saveName')
    h5create([optionsImgStruct.datapath filesep '4RoSEO' filesep 'croppedOffsetSubtractedData4RoSEO' '.h5'],'/SMLM_img',size(SMLM_img))
    h5write([optionsImgStruct.datapath filesep '4RoSEO' filesep 'croppedOffsetSubtractedData4RoSEO' '.h5'],'/SMLM_img',SMLM_img)
else
    SMLM_img = h5read([optionsImgStruct.datapath filesep '4RoSEO' filesep 'croppedOffsetSubtractedData4RoSEO' '.h5'],'/SMLM_img');
end

if isempty(backMatInfo)
    [backg,~,~,~] = n_default.loadImgInteractive(optionsBack{:});
    backg = reshape(backg,[size(backg,1),size(backg,2),size(backg,4)]);
    save([optionsBackStruct.datapath filesep '4RoSEO' filesep 'croppedData4RoSEO.mat'],'n_default','optionsBack','saveName')
    h5create([optionsBackStruct.datapath filesep '4RoSEO' filesep 'croppedData4RoSEO' '.h5'],'/backg',size(backg))
    h5write([optionsBackStruct.datapath filesep '4RoSEO' filesep 'croppedData4RoSEO' '.h5'],'/backg',backg)
else
    backg = h5read([optionsBackStruct.datapath filesep '4RoSEO' filesep 'croppedData4RoSEO' '.h5'],'/backg');
end

%% RoSEO analysis
% load zernikecoefficient
% dataFolderInfo = dir([zernikeFileFolder filesep '* analyzed results.mat']);
dataFolderInfo = dir([zernikeFileFolder filesep '* zernikecoefficient.mat']);
% ExpRecov_xAll = [];
% ExpRecov_yAll = [];
% for dataInd = 1:length(dataFolderInfo)
%     load([zernikeFileFolder filesep dataFolderInfo(dataInd).name],'ExpRecov_x','ExpRecov_y')
%     % remove xy tilt
%     ExpRecov_x(1:2,:) = 0;
%     ExpRecov_y(1:2,:) = 0;
%     % remove defocus
%     ExpRecov_x(4,:) = 0;
%     ExpRecov_y(4,:) = 0;
%         
%     ExpRecov_xAll = [ExpRecov_xAll ExpRecov_x];
%     ExpRecov_yAll = [ExpRecov_yAll ExpRecov_y];
% end
% ExpRecov_x = ExpRecov_xAll;
% ExpRecov_y = ExpRecov_yAll;
ExpRecov_x = [];
ExpRecov_y = [];
for dataInd = 1:length(dataFolderInfo)
    load([zernikeFileFolder filesep dataFolderInfo(dataInd).name],'zernikecoefficient')
    % remove xy tilt
    zernikecoefficient.x(1:2,:) = 0;
    zernikecoefficient.y(1:2,:) = 0;
    % remove defocus
    zernikecoefficient.x(4,:) = 0;
    zernikecoefficient.y(4,:) = 0;
        
    ExpRecov_x = [ExpRecov_x zernikecoefficient.x];
    ExpRecov_y = [ExpRecov_y zernikecoefficient.y];
end
clear zernikecoefficient
zernikecoefficient.x = median(ExpRecov_x,2);
zernikecoefficient.y = median(ExpRecov_y,2);

% zernikecoefficient.x = [];
% zernikecoefficient.y = [];

phasemaskpara.zeroorder=alpha_mean; % add zeroorder leakage into phasemaskpara

n1=Nanoscope('pixelSize',pixelSize,'imageSize',sizeROI,'ADcount',1,'emissWavelength',emitWaveL,...
    'sampleRefractiveIndx',sampleRefractiveIndx,'phasemaskpara',phasemaskpara,'numApt',numApt,...
    'ZernikeCoefficient',zernikecoefficient);

% create PSF matrix accounting for channel transmission ratio
[FPSFx,FPSFy]=n1.createPSFstruct(n1,'ytoxchanneltransratio',channelRatio);

if isempty(frameN)
    frameN = size(SMLM_img,3);
end

% devide image stack into sub-stack for parfor analysis
subLoopNum = frameN/anaLoop;

% handle for recovery
RoSEO_h=@(img,background)RoSEO(n1,img,background,FPSFx,FPSFy,'regVal',reg);
% RoSEO_h=@(img,background)RoSEO_rcnstImg(n1,img,background,FPSFx,FPSFy,'regVal',reg);

loc_data=cell(1,frameN);
% rcnstImg = nan(sizeROI,2*sizeROI,frameN);
p=gcp(); %current pool

cropInd = 0;
cropXAll = zeros(avePSFSize,avePSFSize);
cropYAll = zeros(avePSFSize,avePSFSize);

for l = 1:subLoopNum
    disp('l=')
    disp(l)
    
    SMLM_imgTemp = SMLM_img(:,:, ((l-1)*anaLoop+1) : (l*anaLoop) );
    backgTemp = backg(:,:, ((l-1)*anaLoop+1) : (l*anaLoop) );
    
    for i=1:anaLoop
        F(i)=parfeval(p,RoSEO_h,3,SMLM_imgTemp(:,:,i),backgTemp(:,:,i));
%         F(i)=parfeval(p,RoSEO_h,4,SMLM_imgTemp(:,:,i),backgTemp(:,:,i));
    end
    
    for i=1:anaLoop
        % fetchNext blocks until next results are available.
        [completedIndx,~,~,locD]=fetchNext(F);
%         [completedIndx,~,~,locD,rcnstImgTemp]=fetchNext(F);
        
        loc_data{(l-1)*anaLoop+completedIndx}=locD;
%         rcnstImg(:,:,(l-1)*anaLoop+completedIndx) = rcnstImgTemp;
        
        if ~isempty(locD)
            %             locDX = round( locD(:,2)/pixelSize + (sizeROI+1)/2); % 1/2 is added because sizeROI is odd
            %             locDY = round( locD(:,3)/pixelSize + (sizeROI+1)/2);
%             locDX = round( locD(:,2)/pixelSize + (sizeROI-1)/2-1/2); % changed at 190609
%             locDY = round( locD(:,3)/pixelSize + (sizeROI-1)/2-1/2);
            locDX = floor( (locD(:,2)+sizeROI*pixelSize/2)/pixelSize-1-1/2)+1; % changed at 200409
            locDXYCha = floor( (locD(:,2)+3*sizeROI*pixelSize/2)/pixelSize-1-1/2)+1;
            locDY = floor( (locD(:,3)+sizeROI*pixelSize/2)/pixelSize-1-1/2)+1;
            
            %             locDInd = (  locDX - (avePSFSize-1)/2  ) < 1;
            %             locD(locDInd,:) = [];
            %             locDX(locDInd) = [];
            %             locDY(locDInd) = [];
            %             locDInd = (  locDY - (avePSFSize-1)/2  ) < 1;
            %             locD(locDInd,:) = [];
            %             locDX(locDInd) = [];
            %             locDY(locDInd) = [];
            %             locDInd = (  locDX + (avePSFSize-1)/2  ) > sizeROI;
            %             locD(locDInd,:) = [];
            %             locDX(locDInd) = [];
            %             locDY(locDInd) = [];
            %             locDInd = (  locDY + (avePSFSize-1)/2  ) > sizeROI;
            %             locD(locDInd,:) = [];
            %             locDX(locDInd) = [];
            %             locDY(locDInd) = [];
            locDInd = (  locDX - (avePSFSize-1)/2 - 1/2 ) < 1;% changed at 190613
            locD(locDInd,:) = [];
            locDX(locDInd) = [];
            locDXYCha(locDInd) = [];
            locDY(locDInd) = [];
            locDInd = (  locDY - (avePSFSize-1)/2 - 1/2 ) < 1;
            locD(locDInd,:) = [];
            locDX(locDInd) = [];
            locDXYCha(locDInd) = [];
            locDY(locDInd) = [];
%             locDInd = (  locDX + (avePSFSize-1)/2 - 1/2 ) > sizeROI;
%             locD(locDInd,:) = [];
%             locDX(locDInd) = [];
%             locDY(locDInd) = [];
%             locDInd = (  locDY + (avePSFSize-1)/2 - 1/2 ) > sizeROI;
%             locD(locDInd,:) = [];
%             locDX(locDInd) = [];
%             locDY(locDInd) = [];
            locDInd = (  locDX + (avePSFSize-1)/2 + 1/2 ) > sizeROI; % changed at 200409
            locD(locDInd,:) = [];
            locDX(locDInd) = [];
            locDXYCha(locDInd) = [];
            locDY(locDInd) = [];
            locDInd = (  locDY + (avePSFSize-1)/2 + 1/2 ) > sizeROI;
            locD(locDInd,:) = [];
            locDX(locDInd) = [];
            locDXYCha(locDInd) = [];
            locDY(locDInd) = [];
            
            if ~isempty(locD)
                for ii = 1:size(locD,1)
                    %                     cropXTemp = SMLM_imgTemp( (locDY(ii)-(avePSFSize-1)/2) : (locDY(ii)+(avePSFSize-1)/2),...
                    %                         (locDX(ii)-(avePSFSize-1)/2) : (locDX(ii)+(avePSFSize-1)/2),...
                    %                         i);
                    %                     cropYTemp = SMLM_imgTemp( (locDY(ii)-(avePSFSize-1)/2) : (locDY(ii)+(avePSFSize-1)/2),...
                    %                         (sizeROI+locDX(ii)-(avePSFSize-1)/2) : (sizeROI+locDX(ii)+(avePSFSize-1)/2),...
                    %                         i);
                    cropXTemp = SMLM_imgTemp( (locDY(ii)-(avePSFSize-1)/2) : (locDY(ii)+(avePSFSize-1)/2),...% changed at 190613
                        (locDX(ii)-(avePSFSize-1)/2) : (locDX(ii)+(avePSFSize-1)/2),...
                        i);
%                     cropYTemp = SMLM_imgTemp( (locDY(ii)-(avePSFSize-1)/2) : (locDY(ii)+(avePSFSize-1)/2),...
%                         (sizeROI+locDX(ii)-(avePSFSize-1)/2) : (sizeROI+locDX(ii)+(avePSFSize-1)/2),...
%                         i);
                    cropYTemp = SMLM_imgTemp( (locDY(ii)-(avePSFSize-1)/2) : (locDY(ii)+(avePSFSize-1)/2),...% changed at 200409
                        (locDXYCha(ii)-(avePSFSize-1)/2) : (locDXYCha(ii)+(avePSFSize-1)/2),...
                        i);
                    cropXAll = cropXAll + cropXTemp;
                    cropYAll = cropYAll + cropYTemp;
                    cropInd = cropInd + 1;
                end
            end
        end
    end
    
    clear F   
end

cropXAll = cropXAll./cropInd;
cropYAll = cropYAll./cropInd;
cropCAxis = [min(min(min(cropXAll)),min(min(cropYAll))) max(max(max(cropXAll)),max(max(cropYAll)))];

%% final analysis and SR recunstruction
%--
if workStatInd == 0
    SMLM_save_img = reshape( SMLM_img(:,:,saveSFrame:(saveSFrame+saveFrame-1)), [size(SMLM_img,1),size(SMLM_img,2),saveFrame] );
end
%--
% rawImg = SMLM_img(:,:,1:frameN) - backg(:,:,1:frameN);
% save([dirName filesep saveName ' rawImg,rcnstImg' '.mat'],'rawImg','rcnstImg')
clear SMLM_img SMLM_imgTemp backgTemp %rawImg rcnstImg rcnstImgTemp
save([dirName filesep saveName ' lock_data' '.mat'],'loc_data')

% accumulating estimated M vector
M= [];
sig = [];
frm = [];
locPos = [];
for i = 1:frameN
    if ~isempty(loc_data{i})
        ind = loc_data{i}(:,4) < photonLowThreshold;
        loc_data{i}(ind,:) = [];
        % remove localizations too close to edges
        indL = loc_data{i}(:,2) <= (-((sizeROI-1)/2+0.5)*pixelSize) + remPix*pixelSize;
        indR = loc_data{i}(:,2) >= (((sizeROI-1)/2+0.5)*pixelSize) - remPix*pixelSize;
        indU = loc_data{i}(:,3) <= (-((sizeROI-1)/2+0.5)*pixelSize) + remPix*pixelSize;
        indB = loc_data{i}(:,3) >= (((sizeROI-1)/2+0.5)*pixelSize) - remPix*pixelSize;
        indAll = (indL + indR + indU + indB) > 0;
        loc_data{i}(indAll,:) = [];
        
        M = [M; loc_data{i}(:,5:10)];
        sig = [sig; loc_data{i}(:,4)];
        frm = [frm; i*ones(size(loc_data{i},1),1)];
        locPos = [locPos; loc_data{i}(:,2:3)];
    end
end

% center pixel ind: (size-1)/2+1, shift 0.5 to the negative direction for
% avoiding axis mismatch
Xedges = ( (-((sizeROI-1)/2+0.5)*pixelSize):binSize:(((sizeROI-1)/2+0.5)*pixelSize) );
Yedges = ( (-((sizeROI-1)/2+0.5)*pixelSize):binSize:(((sizeROI-1)/2+0.5)*pixelSize) );

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
title('SR image of RoSEO','FontSize',24)
hold off
h1 = gcf;
export_fig(h1,[dirName filesep saveName ' h1 reconstruction' saveFormat],'-transparent')

figure('Position',P1,'Visible','off');
subplot(1,2,1)
imagesc(cropYAll)
colormap(gray)
axis image ij
caxis(cropCAxis);
colorbar;
title('average Y channel')
subplot(1,2,2)
imagesc(cropXAll)
colormap(gray)
axis image ij
caxis(cropCAxis);
colorbar;
title('average X channel')
h1_1 = gcf;
export_fig(h1_1,[dirName filesep saveName ' h1_1 average PSFs' saveFormat],-'transparent')
close(h1_1)

reconstructed = hist.Values;
reconstructed = reconstructed.';

% ROI selection
binRec = reconstructed > threROI;
binRecFill = imfill(binRec,'holes');
CC = bwconncomp(binRecFill);
numPixels = cellfun(@numel,CC.PixelIdxList);
[~,idx] = max(numPixels);
fibROI = zeros(size(reconstructed));
fibROI(CC.PixelIdxList{idx}) = 1;
fibROIAlpha = fibROI.*0.3;
fibROI = logical(fibROI);
fibROIBW =~fibROI;
pixIdx = 1;
while fibROIBW(pixIdx) == 1
    pixIdx = pixIdx + 1;
end
pixIdx = pixIdx - 1;
if mod(pixIdx,size(fibROIBW,1)) == 0
    pixR = size(fibROIBW,1);
else
    pixR = mod(pixIdx,size(fibROIBW,1));
end
pixC = floor(pixIdx/size(fibROIBW,1))+1;
fibROIBound = bwtraceboundary(fibROIBW,[pixR pixC],'E',4,Inf,'clockwise');
xROI = fibROIBound(:,2);
yROI = fibROIBound(:,1);
recSat = reconstructed;
recSat(reconstructed>indN) = indN;
recSat = round(recSat./indN.*length(colorMap));
I = ind2rgb(recSat,colorMap);
green = cat(3, zeros(size(recSat)),ones(size(recSat)), zeros(size(recSat)));
figure('Position',P1,'Visible','off');
hold on
imagesc(I)
h = imagesc(green);
set(h,'AlphaData',fibROIAlpha);
plot(xROI,yROI,'w','LineWidth',1)
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'Box', 'off');
axis image ij
title({['Selected ROI']; ['transparent green: ROI, white line: boundary']},'FontSize',24)
hold off
hROI1 = gcf;

in = inpolygon(locPos(:,1),locPos(:,2),(xROI-0.5).*binSize - (((sizeROI-1)/2+0.5)*pixelSize),(yROI-0.5).*binSize- (((sizeROI-1)/2+0.5)*pixelSize));

export_fig(hROI1,[dirName filesep saveName ' hROI1, 2Dhist, check ROI vs reconstruction, image' saveFormat],'-transparent')
close(hROI1)
clear hROI1

% save the reconstructed image as a tif file
reconstructed = uint16(reconstructed);
imwrite(reconstructed,[dirName filesep saveName ' reconstructed' '.tif'],'Compression','none')

[reconstructedIn,~,~] = histcounts2(locPos(in,1),locPos(in,2),Xedges,Yedges);
reconstructedIn = reconstructedIn.';
reconstructedIn(reconstructedIn>indN) = indN;
reconstructedIn = imgaussfilt(reconstructedIn,1);

reconstructedIn = reconstructedIn./max(max(reconstructedIn))*255;
reconstructedIn = uint8(reconstructedIn);
imwrite(reconstructedIn,[dirName filesep saveName ' reconstructedIn' '.tif'],'Compression','none')

%--
if workStatInd == 0
    minPho = min(min(min(SMLM_save_img)));
    maxPho = max(max(max(SMLM_save_img)));
    
    figure('Position',P1,'Visible','off');
    for frameInd = 1:saveFrame
        imagesc(SMLM_save_img(:,:,frameInd))
        hold on
        axis image ij;
        caxis([minPho maxPho])
        colorbar
        colormap gray;
        if ~isempty(loc_data{saveSFrame+frameInd-1})
            %             scatter( (loc_data{saveSFrame+frameInd-1}(:,2)/pixelSize+(sizeROI+1)/2 ),...
            %                 (loc_data{saveSFrame+frameInd-1}(:,3)/pixelSize+(sizeROI+1)/2 ),...
            %                 'r+')
            scatter( (loc_data{saveSFrame+frameInd-1}(:,2)/pixelSize+(sizeROI-1)/2-1/2 ),...% changed at 190609 
                (loc_data{saveSFrame+frameInd-1}(:,3)/pixelSize+(sizeROI-1)/2-1/2 ),...
                'r+')
        end
        xlabel('x (pix)')
        ylabel('y (pix)')
        hold off
        h2 = gcf;
        export_fig(h2,[dirName filesep saveName ' h2 pho con raw data with loc, frame=' num2str(saveSFrame+frameInd-1) saveFormat],'-transparent')
    end
    clear h2
end
%--

mLabel = ['Mxx';'Myy';'Mzz';'Mxy';'Mxz';'Myz'];

figure('Position',P1,'Visible','off');
for mInd = 1:6
    subplot(2,3,mInd)
    mHist = histogram(M(:,mInd));
    if mInd < 4
        mHist.BinEdges = 0:0.05:1;
    else
        mHist.BinEdges = -0.5:0.05:0.5;
    end
    xlabel('Estimated value')
    title(mLabel(mInd,:))
end
h3 = gcf;
export_fig(h3,[dirName filesep saveName ' h3 distribution of estimated M'])
close(h3)

figure('Position',P1,'Visible','off');
for mInd = 1:6
    subplot(2,3,mInd)
    mHist = histogram(M(in,mInd));
    if mInd < 4
        mHist.BinEdges = 0:0.05:1;
    else
        mHist.BinEdges = -0.5:0.05:0.5;
    end
    xlabel('Estimated value')
    title(['ROI' mLabel(mInd,:)])
end
h3_1 = gcf;
export_fig(h3_1,[dirName filesep saveName ' h3_1 distribution of estimated M in ROI'])
clear h1 h1_1 h3 h3_1

% calculate Hadamard products of basis for FIM calculation in LS projection
bx.XX = n1.XXxBasis;
bx.YY = n1.YYxBasis;
bx.ZZ = n1.ZZxBasis;
bx.XY = n1.XYxBasis;
bx.XZ = n1.XZxBasis;
bx.YZ = n1.YZxBasis;

by.XX = n1.XXyBasis;
by.YY = n1.YYyBasis;
by.ZZ = n1.ZZyBasis;
by.XY = n1.XYyBasis;
by.XZ = n1.XZyBasis;
by.YZ = n1.YZyBasis;

% Hadamard products, for x channel
Bx.aa = (bx.XX).*(bx.XX);
Bx.ab = (bx.XX).*(bx.YY);
Bx.ac = (bx.XX).*(bx.ZZ);
Bx.ad = (bx.XX).*(bx.XY);
Bx.ae = (bx.XX).*(bx.XZ);
Bx.af = (bx.XX).*(bx.YZ);

Bx.bb = (bx.YY).*(bx.YY);
Bx.bc = (bx.YY).*(bx.ZZ);
Bx.bd = (bx.YY).*(bx.XY);
Bx.be = (bx.YY).*(bx.XZ);
Bx.bf = (bx.YY).*(bx.YZ);

Bx.cc = (bx.ZZ).*(bx.ZZ);
Bx.cd = (bx.ZZ).*(bx.XY);
Bx.ce = (bx.ZZ).*(bx.XZ);
Bx.cf = (bx.ZZ).*(bx.YZ);

Bx.dd = (bx.XY).*(bx.XY);
Bx.de = (bx.XY).*(bx.XZ);
Bx.df = (bx.XY).*(bx.YZ);

Bx.ee = (bx.XZ).*(bx.XZ);
Bx.ef = (bx.XZ).*(bx.YZ);

Bx.ff = (bx.YZ).*(bx.YZ);

% for y channel
By.aa = (by.XX).*(by.XX);
By.ab = (by.XX).*(by.YY);
By.ac = (by.XX).*(by.ZZ);
By.ad = (by.XX).*(by.XY);
By.ae = (by.XX).*(by.XZ);
By.af = (by.XX).*(by.YZ);

By.bb = (by.YY).*(by.YY);
By.bc = (by.YY).*(by.ZZ);
By.bd = (by.YY).*(by.XY);
By.be = (by.YY).*(by.XZ);
By.bf = (by.YY).*(by.YZ);

By.cc = (by.ZZ).*(by.ZZ);
By.cd = (by.ZZ).*(by.XY);
By.ce = (by.ZZ).*(by.XZ);
By.cf = (by.ZZ).*(by.YZ);

By.dd = (by.XY).*(by.XY);
By.de = (by.XY).*(by.XZ);
By.df = (by.XY).*(by.YZ);

By.ee = (by.XZ).*(by.XZ);
By.ef = (by.XZ).*(by.YZ);

By.ff = (by.YZ).*(by.YZ);

sizeVBackg = size(backg,1);
backgMeanX = sum(sum(sum(backg(:,1:sizeVBackg,:),3),2),1)/numel(backg(:,1:sizeVBackg,:));
backgMeanY = sum(sum(sum(backg(:,(sizeVBackg+1):(2*sizeVBackg),:),3),2),1)/numel(backg(:,(sizeVBackg+1):(2*sizeVBackg),:));
backgMean = backgMeanX+backgMeanY;

figure('Position',P1,'Visible','off');
histogram(sig);
xlabel('estimated emitted photons (maximum collectable photons if molecules oriented in-plane)','FontSize',18)
title(  [ 'Photon number, median=' num2str( round(median(sig),3,'significant') ) ', mean=' num2str( round(mean(sig),3,'significant') )...
    ', mean of backg (X+Y)=' num2str( round(backgMean,2,'significant'))]  );
h4 = gcf;
export_fig(h4,[dirName filesep saveName ' h4 histogram estimated emitted photons' saveFormat],'-transparent')
close(h4)
clear h4

clear backg n1 n_default index

anaTime = toc;

save([dirName filesep saveName ' all analyzed results' '.mat'])

if workStatInd == 1
    exit;
end