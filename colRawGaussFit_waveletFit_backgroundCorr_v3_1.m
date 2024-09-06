% 190604 - v3_1 TD - Changed saving format from uint16 to single. Also, use
% Tiff class for saving the floating point numbers.

% 190522 - v3 TD - Added "margin" in the fitting procedure

% 190226 - v2 TD - Analysis on multiple stacks. Some code tuning for
% compatibility with workstation (ll03) analysis

% Tianben Ding 190212
% This program estimate non-uniform background by fitting individual
% columns and raws by 1D 2 Gaussian. A wavelet filter is further applied to
% the sum.

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

% Input parameters
% Image address
imageFile = ...
    'G:\RawData\190729_2158BBHPTara laser irradiation tests\offsetSubtracted';
imageFilePrefix = 'Data'; % prefix of your data name
imageFileInd = [12]; % data index you want to analyze. it accepts multiple data set.
% File format
fFormat = '.tif';

% Vertical image size
imageSizeV = 200;% pixel
% Horizontal image size
imageSizeH = 2048;% pixel
% Horisonzal boundary
imageBoundH = 1024;% pixel

% holizontal margin during fitting procedure, from the left and right in
% each channel
marginH = 250;
marginV = 0;

aveRate = 200; % frame

sigmaMin = 80; % minimum sigma to be fit, pixel

% wavelet filtering para
sorh = 's';

wname = 'bior6.8';
level = 6;

% Input part end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% tic
xL = 1:imageBoundH;
xLFit = (1+marginH):(imageBoundH-marginH);
yL = 1:imageSizeV;
yLFit = (1+marginV):(imageSizeV-marginV);

xR = imageBoundH+1:imageSizeH;
xRFit = (imageBoundH+1+marginH):(imageSizeH-marginH);
yR = yL;
yRFit = yLFit;

% tags of TIFF class for writing single estimated background into tiff
% files
tagstruct.ImageLength = imageSizeV;
tagstruct.ImageWidth = imageSizeH;
tagstruct.Compression = Tiff.Compression.None;
tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample = 32;
tagstruct.SamplesPerPixel = 1; % save individual estimated background in single tifs
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;


for stackInd = 1:length(imageFileInd)
    saveBGFile = [imageFile filesep imageFilePrefix num2str(imageFileInd(stackInd)) filesep 'BGEst_lineFitWavelet'];
    
    if exist(saveBGFile,'dir') ~= 7
        mkdir(saveBGFile);
    end
    
    % Read the images
    dataFolderInfo = dir([imageFile filesep imageFilePrefix num2str(imageFileInd(stackInd)) filesep '*' fFormat]); % read folder information of raw data
    numFiles = length(dataFolderInfo); % store the frame number of the folder
    
    totSubStack = numFiles/aveRate;
    
    for hSS = 1:totSubStack
        colFitL = nan(length(yLFit),length(xLFit));
        colFitR = nan(length(yRFit),length(xRFit));
        rowFitL = nan(length(yLFit),length(xLFit));
        rowFitR = nan(length(yRFit),length(xRFit));
        capturedImageSum = zeros(imageSizeV,imageSizeH);
        for h = ( (hSS-1)*aveRate+1 ):(hSS*aveRate)
            capturedImage = tiffread2([imageFile filesep imageFilePrefix num2str(imageFileInd(stackInd)) filesep dataFolderInfo(h).name],1,1);
            capturedImage = double(capturedImage.data);
            capturedImageSum = capturedImageSum + capturedImage;
        end
        aveImage = capturedImageSum/aveRate;
        aveImageLFit = aveImage(yLFit,xLFit);
        aveImageRFit = aveImage(yRFit,xRFit);
        
        parfor i = 1:length(yLFit) % scan from row#1 to row#end, left channel
            options = fitoptions('gauss2', 'Lower', [0 -Inf sigmaMin 0 -Inf sigmaMin]);
            f = fit(xLFit.',aveImageLFit(i,:).','gauss2',options);
            a1 = f.a1;
            b1 = f.b1;
            c1 = f.c1;
            a2 = f.a2;
            b2 = f.b2;
            c2 = f.c2;
            %         if c1 < sigmaMin
            %             a1 = 0;
            %         end
            %         if c2 < sigmaMin
            %             a2 = 0;
            %         end
            rowFitL(i,:) =  a1.*exp( -((xLFit-b1)./c1).^2 )  +  a2.*exp( -((xLFit-b2)./c2).^2 );
        end
        parfor i = 1:length(yRFit) % scan from row#1 to row#end, right channel
            options = fitoptions('gauss2', 'Lower', [0 -Inf sigmaMin 0 -Inf sigmaMin]);
            f = fit(xRFit.',aveImageRFit(i,:).','gauss2',options);
            a1 = f.a1;
            b1 = f.b1;
            c1 = f.c1;
            a2 = f.a2;
            b2 = f.b2;
            c2 = f.c2;
            %         if c1 < sigmaMin
            %             a1 = 0;
            %         end
            %         if c2 < sigmaMin
            %             a2 = 0;
            %         end
            rowFitR(i,:) =  a1.*exp( -((xRFit-b1)./c1).^2 )  +  a2.*exp( -((xRFit-b2)./c2).^2 );
        end
        parfor i = 1:length(xLFit) % scan from column#1 to column#end, left channel
            options = fitoptions('gauss2', 'Lower', [0 -Inf sigmaMin 0 -Inf sigmaMin]);
            f = fit(yLFit.',aveImageLFit(:,i),'gauss2',options);
            a1 = f.a1;
            b1 = f.b1;
            c1 = f.c1;
            a2 = f.a2;
            b2 = f.b2;
            c2 = f.c2;
            %         if c1 < sigmaMin
            %             a1 = 0;
            %         end
            %         if c2 < sigmaMin
            %             a2 = 0;
            %         end
            colFitL(:,i) =  a1.*exp( -((yLFit-b1)./c1).^2 )  +  a2.*exp( -((yLFit-b2)./c2).^2 );
        end
        parfor i = 1:length(xRFit) % scan from column#1 to column#end, right channel
            options = fitoptions('gauss2', 'Lower', [0 -Inf sigmaMin 0 -Inf sigmaMin]);
            f = fit(yRFit.',aveImageRFit(:,i),'gauss2',options);
            a1 = f.a1;
            b1 = f.b1;
            c1 = f.c1;
            a2 = f.a2;
            b2 = f.b2;
            c2 = f.c2;
            %         if c1 < sigmaMin
            %             a1 = 0;
            %         end
            %         if c2 < sigmaMin
            %             a2 = 0;
            %         end
            colFitR(:,i) =  a1.*exp( -((yRFit-b1)./c1).^2 )  +  a2.*exp( -((yRFit-b2)./c2).^2 );
        end
        ave1DFitL = (rowFitL + colFitL)./2;
        ave1DFitR = (rowFitR + colFitR)./2;
        
        %     backgEstL = Wavelet_backg_est_noParfor(ave1DFit(:,xL),thresh_bg,wavelet_level,'bior3.5',num_iter);
        %     backgEstR = Wavelet_backg_est_noParfor(ave1DFit(:,xR),thresh_bg,wavelet_level,'bior3.5',num_iter);
        %     backgEstBoth = [backgEstL backgEstR];
        %     backgroundEstAll = uint16(backgEstBoth);
        
        [CL,SL] = wavedec2(ave1DFitL,level,wname);
        thrL = wthrmngr('dw2ddenoLVL','sqtwolog',CL,SL,'one');
        [XDENL,~,~] = wdencmp('lvd',CL,SL,wname,level,thrL,sorh);
        
        [CR,SR] = wavedec2(ave1DFitR,level,wname);
        thrR = wthrmngr('dw2ddenoLVL','sqtwolog',CR,SR,'one');
        [XDENR,~,~] = wdencmp('lvd',CR,SR,wname,level,thrR,sorh);
        
        %         [thr,sorh,keepapp] = ddencmp('den','wv',ave1DFitL);
        %         XDENL = wdencmp('gbl',ave1DFitL,'sym4',2,thr,sorh,keepapp);
        %         [thr,sorh,keepapp] = ddencmp('den','wv',ave1DFitR);
        %         XDENR = wdencmp('gbl',ave1DFitR,'sym4',2,thr,sorh,keepapp);
        
        backgroundEstAll = zeros(imageSizeV,imageSizeH);
        backgroundEstAll(yLFit,xLFit) = XDENL;
        backgroundEstAll(yRFit,xRFit) = XDENR;
        %         backgroundEstAll = uint16(backgroundEstAll);
        backgroundEstAll = single(backgroundEstAll);
        
        
        %             figure;subplot(3,1,1);
        %             imagesc(ave1DFit); colormap gray;colorbar; axis image ij;
        %             title('Noisy Image');
        %             subplot(3,1,2);
        %             imagesc(XDEN); colormap gray;colorbar; axis image ij;
        %             title('Denoised Image');
        %             subplot(3,1,3);
        %             imagesc(backgroundEstAll); colormap gray;colorbar; axis image ij;
        %             title('Rounded Image');
        
        parfor h = ( (hSS-1)*aveRate+1 ):(hSS*aveRate)
            %             imwrite(backgroundEstAll,[saveBGFile filesep dataFolderInfo(h).name],'Compression','none');
            t = Tiff([saveBGFile filesep dataFolderInfo(h).name],'w');
            t.setTag(tagstruct);
            t.write(backgroundEstAll);
            t.close();
        end
    end
end
% toc