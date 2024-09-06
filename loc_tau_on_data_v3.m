% 200830 - TD - Also collect ROI only data and inplane angle deviation from
% fibril backbones


% 200723 Tianben Ding
% Construct loc and tau_on data.mat

%% Analysis configuration
clear; close all;clc;

% Get output figure size (for large figures)
P0 = get(0,'ScreenSize');
P1(3:4) = P0(3:4)-150;

% Default figure setting
set(0, 'DefaultFigureColor','White',...
    'DefaultTextFontName','Arial','DefaultAxesFontName','Arial',...
    'DefaultAxesFontSize', 18,'DefaultTextFontSize',18,...
    'DefaultTextColor','k',...
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
%% Input parameters - fixed parameters

% data path --
dataFolderPath = '210709_readRidgeDetCalDisorder_v22';
fileName = '190729 Data12 reg0.25 lineFitWavelet,200 photonThre, nearNeiDis=30, thThre=70';

dirNamePre = '';
saveNameFinal = ['190729 Data12 loc and tau_on data'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get information for saving
t=datetime('today');

dirNameFinal=datestr(t,'yymmdd');
dirNameFinal=[dirNamePre dirNameFinal '_' mfilename];

if exist(dirNameFinal,'dir') ~= 7
    mkdir(dirNameFinal);
end

save([dirNameFinal filesep saveNameFinal ' initial setting' '.mat'])
dirNameFinal1 = dirNameFinal;
saveNameFinal1 = saveNameFinal;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([dataFolderPath filesep fileName ' all analyzed results.mat'])

if length(frmROI) > sum(validThetaID)
    frmROI = frmROI(validThetaID);
end
if length(frmOFF) > sum(validThetaIDOFF)
    frmOFF = frmOFF(validThetaIDOFF);
end
frm = [frmROI(~invalidLocROI); frmOFF];
locPos = [locPosROI(~invalidLocROI,:); locPosOFF];
sig = [sigROI(~invalidLocROI); sigOFF];
M = [MROI(~invalidLocROI,:); MOFF];
theta = [thetaROI(~invalidLocROI); thetaOFF];
phi = [phiROI(~invalidLocROI); phiOFF];
omega = [omegaROI(~invalidLocROI); omegaOFF]/2;

loc_data_ana = [frm locPos sig M theta phi omega];
tau_on_fibrils = [locPosGroup onTimeCorrAll];

loc_data_ana_ROI = [frmROI(~invalidLocROI) locPosROI(~invalidLocROI,:) sigROI(~invalidLocROI) ...
    MROI(~invalidLocROI,:) thetaROI(~invalidLocROI) phiROI(~invalidLocROI) deg2rad(psi2RefDegInplaneS(~invalidLocROI)) omegaROI(~invalidLocROI)/2];
%% clear fig* SMLM_imgTemp allVectorDis brightness_* Bx By index green I n1 nearNeiInd
save([dirNameFinal1 filesep saveNameFinal1 '.mat'],'loc_data_ana','tau_on_fibrils','loc_data_ana_ROI','sizeROI','pixelSize')