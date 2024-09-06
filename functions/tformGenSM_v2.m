% 200319 - v2 TD - Change tiffread2 to Tiff class. Reading an image
% sequence with Tiff is much faster than tiffread2. See resutls in
% \Detector calibration, offset subtraction, data
% transfer\Analysis\200319_readImgSpeed_v1.

% 190718 - TD - Changed the bin size definition in distance histogram of
% pared localizations

% 190415 Tianben Ding
% tform generator based on SM localizations
% This code is based on SMRegisterReconstLargeEnsRegine_v22.m.

function  [tformSML2R,tformSMR2L,fig] =...
    tformGenSM_v2(num,pixelSize,conversionF,channelRatio,tformL2R,tformR2L,zoomRegion,...
    dataFolderPath,standardPSFData,binSize,regPara,visuWinSize,visInd)

boundaryXPosSM = regPara.boundaryXPosSM;
squareS = regPara.squareS;
photonLowThresholdReg =regPara.photonLowThresholdReg;
preThre = regPara.preThre;
bmin = regPara.bmin;
bmax = regPara.bmax;
biasThre = regPara.biasThre;
ensFactor = regPara.ensFactor;

dataFolderInfo = dir([dataFolderPath filesep 'offsetSubtracted' filesep standardPSFData filesep '*.tif' ]);

% Initialize parameters for the following process
xPosRightPreRegAll = [];
yPosRightPreRegAll = [];
xPosLeftPreRegAll = [];
yPosLeftPreRegAll = [];

precRightAll = [];
precLeftAll = [];

allVectorDisAll = [];
allVectorRadAll = [];
allVectorPosAll = [];
allVectorPrecAll = [];

for h = 1:length(dataFolderInfo)
    id = transpose(num(num(:,1)==h,1));
    xPos = transpose(num(num(:,1)==h,2));
    yPos = transpose(num(num(:,1)==h,3));
    sigma = transpose(num(num(:,1)==h,4));
    
    if ~isempty(id)
        % Read the images
        t = Tiff([dataFolderPath filesep 'offsetSubtracted' filesep standardPSFData filesep dataFolderInfo(h).name]);
        capturedImage = t.read();
        t.close();
        capturedImage = double(capturedImage);
%         capturedImage = tiffread2([dataFolderPath filesep 'offsetSubtracted' filesep standardPSFData filesep dataFolderInfo(h).name],1,1);
%         capturedImage = double(capturedImage.data);
        
        xPosCopy = xPos;
        yPosCopy = yPos;
        for i = 1:size(xPos,1)
            distMat = sqrt((xPos(:)-xPos(i)).^2+(yPos(:)-yPos(i)).^2);
            ind = distMat < max(squareS,squareS)*pixelSize;
            if sum(ind) > 1 % if there is too close emitters
                % remove too close emitters
                xPosCopy(ind) = nan;
                yPosCopy(ind) = nan;
                sigma(ind) = nan;
            end
        end
        xPos = xPosCopy;
        yPos = yPosCopy;
        
        xPos(isnan(xPos)) = [];
        yPos(isnan(yPos)) = [];
        sigma(isnan(sigma)) = [];
        
        % Subtraction of background
        % insideMatrix: contains all emitter candidates' PSF information
        insideMatrix = ...
            capturedImage(reshape(repmat(floor(yPos./pixelSize)+1,[squareS,1]),[1,(squareS)*length(yPos)])+repmat(-(squareS-1)/2:(squareS-1)/2,[1,length(yPos)]),...
            reshape(repmat(floor(xPos./pixelSize)+1,[squareS,1]),[1,(squareS)*length(xPos)])+repmat(-(squareS-1)/2:(squareS-1)/2,[1,length(xPos)])...
            );
        insideIndex = kron(eye(length(xPos)),ones(squareS,squareS));
        insideMatrix(~logical(insideIndex)) = 0;
        insideVector = sum(insideMatrix,1);
        insideVector = reshape(insideVector,[squareS,length(xPos)]);
        insideVector = sum(insideVector,1);
        
        % largeMatrix: contains all emitter candidates' PSF information (larger region)
        largeMatrix = ...
            capturedImage(reshape(repmat(floor(yPos./pixelSize)+1,[squareS+2,1]),[1,(squareS+2)*length(yPos)])+repmat(-(squareS+2-1)/2:(squareS+2-1)/2,[1,length(yPos)]),...
            reshape(repmat(floor(xPos./pixelSize)+1,[squareS+2,1]),[1,(squareS+2)*length(xPos)])+repmat(-(squareS+2-1)/2:(squareS+2-1)/2,[1,length(xPos)])...
            );
        largeIndex = kron(eye(length(xPos)),ones(squareS+2,squareS+2));
        largeMatrix(~logical(largeIndex)) = 0;
        largeVector = sum(largeMatrix,1);
        largeVector = reshape(largeVector,[squareS+2,length(xPos)]);
        largeVector = sum(largeVector,1);
        
        backgroundVector = largeVector - insideVector;
        backgroundVector = backgroundVector./(squareS*2+(squareS+2)*2);
        
        % backgroundMatrix: contains background information at all emitters
        % candidates' positions
        backgroundMatrix = reshape(repmat(backgroundVector,[length(xPos)*squareS*squareS,1]),[length(xPos)*squareS,length(xPos)*squareS]);
        backgroundMatrix(~logical(insideIndex)) = 0;
        
        % Subtract background from inside box
        % intensMatrix: contains all emitter candidates' PSF information after
        % background subtraction
        intensMatrix = insideMatrix - backgroundMatrix;
        
        intensAllVector = sum(intensMatrix,1);
        intensAllVector = reshape(intensAllVector,[squareS,length(xPos)]); % not rounded yet
        intensAllVector = sum(intensAllVector,1).*conversionF;
        backgroundAllVector = backgroundVector.*conversionF;
        
        % Remove localizations with negative intensities
        xPos(intensAllVector<=0) = [];
        yPos(intensAllVector<=0) = [];
        sigma(intensAllVector<=0) = [];
        backgroundAllVector(intensAllVector<=0) = [];
        intensAllVector(intensAllVector<=0) = [];
                
        % Distinguish two channels
        rightInd = xPos > boundaryXPosSM*pixelSize;
        leftInd = xPos <= boundaryXPosSM*pixelSize;
        
        % Store the right channel information
        xPosRight = xPos(rightInd);
        yPosRight = yPos(rightInd);
        sigmaRight = sigma(rightInd);
        intensRight = intensAllVector(rightInd);
        backgroundRight = backgroundAllVector(rightInd);
        if channelRatio > 1
            intensRight = intensRight.*channelRatio;
            backgroundRight = backgroundRight.*channelRatio;
        end
        photonRight = round(intensRight);
        backgroundRight = round(backgroundRight);
        
        % Remove localization if its photon number is smaller than
        % threshold (registration)
        ind = photonRight <= photonLowThresholdReg;
        xPosRight(ind) = [];
        yPosRight(ind) = [];
        sigmaRight(ind) = [];
        photonRight(ind) = [];
        backgroundRight(ind) = [];
        
        % Store the left channel information
        xPosLeft = xPos(leftInd);
        yPosLeft = yPos(leftInd);
        sigmaLeft = sigma(leftInd);
        intensLeft = intensAllVector(leftInd);
        backgroundLeft = backgroundAllVector(leftInd);
        if channelRatio <= 1
            intensLeft = intensLeft./channelRatio;
            backgroundLeft = backgroundLeft./channelRatio;
        end
        photonLeft = round(intensLeft);
        backgroundLeft = round(backgroundLeft);
        
        % Remove localization if its photon number is smaller than
        % threshold (registration)
        ind = photonLeft <= photonLowThresholdReg;
        xPosLeft(ind) = [];
        yPosLeft(ind) = [];
        sigmaLeft(ind) = [];
        photonLeft(ind) = [];
        backgroundLeft(ind) = [];
        
        % Calculate localization precision based on captured photon number
        % (Least square)
        tauRight = 2*pi*backgroundRight.*(sigmaRight.^2+pixelSize^2/12)./photonRight./(pixelSize^2);
        precRight = (sigmaRight.^2+pixelSize^2/12)./photonRight.*(16/9+4*tauRight);
        precRight = sqrt(precRight);
        
        tauLeft = 2*pi*backgroundLeft.*(sigmaLeft.^2+pixelSize^2/12)./photonLeft./(pixelSize^2);
        precLeft = (sigmaLeft.^2+pixelSize^2/12)./photonLeft.*(16/9+4*tauLeft);
        precLeft = sqrt(precLeft);
        
        % thresholding using the calculated precision
        ind = precRight > preThre;
        xPosRight(ind) = [];
        yPosRight(ind) = [];
        precRight(ind) = [];
        precRightAll = [precRightAll precRight];
        
        ind = precLeft > preThre;
        xPosLeft(ind) = [];
        yPosLeft(ind) = [];
        precLeft(ind) = [];
        precLeftAll = [precLeftAll precLeft];
        
        if ~isempty(xPosRight) && ~isempty(xPosLeft)
            xPosRight = xPosRight-pixelSize*boundaryXPosSM;
            xPosLeft = -xPosLeft+pixelSize*boundaryXPosSM;
            
            xPosRightPreRegAll = [xPosRightPreRegAll xPosRight];
            yPosRightPreRegAll = [yPosRightPreRegAll yPosRight];
            xPosLeftPreRegAll = [xPosLeftPreRegAll xPosLeft];
            yPosLeftPreRegAll = [yPosLeftPreRegAll yPosLeft];
            
            [xPosLeftTrans,yPosLeftTrans] = transformPointsInverse(tformL2R,xPosLeft.',yPosLeft.');
            xPosLeftTrans = xPosLeftTrans.';
            yPosLeftTrans = yPosLeftTrans.';
            
            rightMeshX = repmat(xPosRight,length(xPosLeftTrans),1);
            rightMeshY = repmat(yPosRight,length(yPosLeftTrans),1);
            leftMeshX = repmat(xPosLeft.',1,length(xPosRight));
            leftMeshY = repmat(yPosLeft.',1,length(yPosRight));
            leftTransMeshX = repmat(xPosLeftTrans.',1,length(xPosRight));
            leftTransMeshY = repmat(yPosLeftTrans.',1,length(yPosRight));
            precRightMesh = repmat(precRight,length(precLeft),1);
            precLeftMesh = repmat(precLeft.',1,length(precRight));
            
            % calculate Euclidean distance between all SMs in two channels
            allVectorDis = hypot(rightMeshX-leftTransMeshX,rightMeshY-leftTransMeshY);
            allVectorCos = (rightMeshX-leftTransMeshX)./allVectorDis;
            allVectorSin = (rightMeshY-leftTransMeshY)./allVectorDis;
            allVectorRad = angle(allVectorCos + 1j*allVectorSin);
            allVectorDisAll = [allVectorDisAll reshape(allVectorDis,[1,size(allVectorDis,1)*size(allVectorDis,2)])];
            allVectorRadAll = [allVectorRadAll reshape(allVectorRad,[1,size(allVectorRad,1)*size(allVectorRad,2)])];
            allVectorPosAll = [allVectorPosAll ...
                [reshape(rightMeshX,[1,size(rightMeshX,1)*size(rightMeshX,2)]);...
                reshape(rightMeshY,[1,size(rightMeshY,1)*size(rightMeshY,2)]);...
                reshape(leftMeshX,[1,size(leftMeshX,1)*size(leftMeshX,2)]);...
                reshape(leftMeshY,[1,size(leftMeshY,1)*size(leftMeshY,2)])]];
            allVectorPrecAll = [allVectorPrecAll ...
                [reshape(precRightMesh,[1,size(precRightMesh,1)*size(precRightMesh,2)]);...
                reshape(precLeftMesh,[1,size(precLeftMesh,1)*size(precLeftMesh,2)])]];
        end
    end
end

[xPosLeftPreRegAllTrans, yPosLeftPreRegAllTrans] = transformPointsInverse(tformL2R,xPosLeftPreRegAll.',yPosLeftPreRegAll.');
xPosLeftPreRegAllTrans = xPosLeftPreRegAllTrans.';
yPosLeftPreRegAllTrans = yPosLeftPreRegAllTrans.';

if visInd == 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % all localization used for registration check (before pairing, with tform, only ones with high precision)
    figure('Position',visuWinSize);
    hold on;
    scatter(xPosRightPreRegAll,yPosRightPreRegAll,'LineWidth',1)
    scatter(xPosLeftPreRegAllTrans, yPosLeftPreRegAllTrans,'LineWidth',1)
    axis image ij
    axis(zoomRegion)
    grid on
    grid minor
    legend('rightPoints','leftPointsTrans')
    xlabel('x position [nm]')
    ylabel('y position [nm]')
    title({'Geometric transform of SM images';...
        ['before pairing, localizations with precision < ' num2str(preThre) ' nm']})
    hold off
    fig.h1 = gcf;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % localization precision of the all localization (histogram, before pairing, only ones with high precision)
    figure('Position',visuWinSize);
    hold on;
    hPR = histogram(precRightAll);
    hPL = histogram(precLeftAll);
    hPL.BinWidth = hPR.BinWidth;
    legend(['right cha, median = ' num2str(median(precRightAll)) ' nm'],['left cha, median = ' num2str(median(precLeftAll)) ' nm'])
    xlabel('localization precision [nm]')
    title({'Localization precision distribution of emitters';...
        ['after ' num2str(photonLowThresholdReg) ' photon and ' num2str(preThre) ' nm precision thresholding']})
    hold off
    fig.h2 = gcf;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% all possible pairing distance (histogram, only ones with high precision)
figure('Position',visuWinSize);
hold on;
if ~isempty(bmin) && ~isempty(bmax)
    histDist = histogram(allVectorDisAll,'BinLimits',[bmin,bmax]);
else
    histDist = histogram(allVectorDisAll);
end
[~,histDistMaxInd] = max(histDist.Values);
histDist.BinEdges = 0:histDist.BinEdges(histDistMaxInd+1)/15:histDist.BinEdges(histDistMaxInd+1)*2;
%histDist.BinEdges = 0:binSize/2:histDist.BinEdges(histDistMaxInd+1)*2; %190718
% Gaussian equation for one dimentional fitting
gaussEqn = 'a*exp(-(x-b)^2/(2*c^2))+d';
x =histDist.BinEdges+histDist.BinWidth/2;
x = x(1:end-1);
[maxi,maxiInd] = max(histDist.Values);
startPoints = [maxi x(maxiInd) binSize 0];
f = fit(x.',(histDist.Values).',gaussEqn,'Start', startPoints);
plot(f,x,(histDist.Values));
xlabel('Pairing distance [nm]')
title({'Distribution of all possible pairing distance';...
    ['after ' num2str(photonLowThresholdReg) ' photon and ' num2str(preThre) ' nm precision thresholding']})
hold off
fig.h3 = gcf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disPeak = f.b;
disSigma = (f.c)/ensFactor;

if visInd == 0
    close(fig.h3)
end

if disPeak < hypot(biasThre,biasThre)
    % If there is little bias on the registration, go forward with the
    % calibratyed transform function
    tformSML2R = tformL2R;
    tformSMR2L = tformR2L;
else
    % If there is any obvious linear bias, compensate it and re-generate a
    % new transform function
    indL = allVectorDisAll > (disPeak - disSigma);
    indU = allVectorDisAll < (disPeak + disSigma);
    indLU = (indL+indU) > 1;
    
    % Radian filtering
    workingRad = allVectorRadAll(indLU);
    resVec = sum(exp(1j*workingRad));
    meanRad = angle(resVec);
    xRadMin = meanRad - pi/2;
    if xRadMin < -pi
        workingRad((pi-abs(-pi-xRadMin)) < workingRad) = workingRad((pi-abs(-pi-xRadMin)) < workingRad) - 2*pi;
    end
    xRadMax = meanRad + pi/2;
    if xRadMax > pi
        workingRad(workingRad < (-pi+abs(xRadMax-pi))) = workingRad(workingRad < (-pi+abs(xRadMax-pi))) + 2*pi;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Aangles between all possible pairing localization, after distance thresholding (histogram)
    figure('Position',visuWinSize);
    hold on;
    histRad = histogram(workingRad);
    histRad.BinEdges = xRadMin:0.05:xRadMax;
    % Gaussian equation for one dimentional fitting
    gaussEqn = 'a*exp(-(x-b)^2/(2*c^2))+d';
    x =  histRad.BinEdges+histRad.BinWidth/2;
    x = x(1:end-1);
    [maxRad,maxRadInd] = max(histRad.Values);
    startPointsRad = [maxRad x(maxRadInd) pi/4 0];
    fRad = fit(x.',(histRad.Values).',gaussEqn,'Start',startPointsRad);
    plot(fRad,x,(histRad.Values));
    xlabel('Radian')
    title({['Angles of all possible pairs, after ' num2str(photonLowThresholdReg) ' photon,'];...
        [num2str(preThre) ' nm precision, and distance thresholding']})
    hold off
    fig.h4 = gcf;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    radPeak = fRad.b;
    radSigma = (fRad.c)/ensFactor;
    indRadU = allVectorRadAll < (radPeak + radSigma);
    indRadL = allVectorRadAll > (radPeak - radSigma);
    indLUAll = (indL+indU+indRadL+indRadU) > 3;
    
    pairedPos = allVectorPosAll(:,indLUAll);
    pairedPrec = allVectorPrecAll(:,indLUAll);
    
    % Fit the old geometric transformation to the control point pairs
    movingPoints = pairedPos(1:2,:).';
    fixedPoints = pairedPos(3:4,:).';
    fixedPointsTransUncorr = transformPointsInverse(tformL2R,fixedPoints);
    medXDir = median(movingPoints(:,1)-fixedPointsTransUncorr(:,1));
    medYDir = median(movingPoints(:,2)-fixedPointsTransUncorr(:,2));
    
    if visInd == 0
        close(fig.h4)
    end
    if visInd == 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % test figure transform without correction (paired localizations)
        figure('Position',visuWinSize);
        hold on;
        scatter(movingPoints(:,1),movingPoints(:,2),'LineWidth',1)
        scatter(fixedPointsTransUncorr(:,1),fixedPointsTransUncorr(:,2),'LineWidth',1)
        for k = 1:size(movingPoints,1)
            plot([movingPoints(k,1) fixedPointsTransUncorr(k,1)],[movingPoints(k,2) fixedPointsTransUncorr(k,2)],'-g');
        end
        axis image ij
        axis(zoomRegion)
        % axis(workingRegion)
        grid on
        grid minor
        legend('rightPoints','leftPointsTrans')
        xlabel('x position [nm]')
        ylabel('y position [nm]')
        title({'Geometric transform of SM images';...
            ['without corrected transform, Paired # = ' num2str(size(movingPoints,1))]})
        hold off
        fig.h5 = gcf;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % localization precision used for pairing (histogram)
        figure('Position',visuWinSize);
        hold on;
        histPrec1 = histogram(pairedPrec(1,:));
        histPrec2 = histogram(pairedPrec(2,:));
        histPrec2.BinWidth = histPrec1.BinWidth;
        legend(['right cha, #' num2str(length(pairedPrec(1,:))) ' median = ' num2str(median(pairedPrec(1,:))) ' nm'],...
            ['left cha, #' num2str(length(pairedPrec(2,:))) ' median = ' num2str(median(pairedPrec(2,:))) ' nm'])
        xlabel('localization precision [nm]')
        % ylabel('Frequency')
        title('Localization precision used for pairing and re-calibration')
        hold off
        fig.h6 = gcf;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % test figure transform without correction (histogram)
        figure('Position',visuWinSize);
        hold on;
        histX = histogram(movingPoints(:,1)-fixedPointsTransUncorr(:,1));
        histY = histogram(movingPoints(:,2)-fixedPointsTransUncorr(:,2));
        histY.BinWidth = histX.BinWidth;
        legend(['xDir, median:' num2str(medXDir) 'nm'],...
            ['yDir, median:' num2str(medYDir) 'nm'])
        xlabel('Error [nm]')
        % ylabel('Frequency')
        title({'Error distribution of geometric transform';...
            'without corrected transform'})
        hold off
        fig.h7 = gcf;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        lengthReg = hypot((movingPoints(:,1)-fixedPointsTransUncorr(:,1)),(movingPoints(:,2)-fixedPointsTransUncorr(:,2)));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % test figure transform without correction (histogram of distance)
        figure('Position',visuWinSize);
        hold on;
        histogram(lengthReg);
        legend(['median= ' num2str(median(lengthReg)) ' nm']);
        xlabel('Error [nm]')
        % ylabel('Frequency')
        title({'Error distribution of geometric transform';...
            'without corrected transform (distance)'})
        hold off
        fig.h8 = gcf;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    tformSML2R = fitgeotrans(movingPoints,fixedPoints,'polynomial',2);
    tformSMR2L = fitgeotrans(fixedPoints,movingPoints,'polynomial',2);
    
    if visInd == 1
        % Check the quality of fitgeotrans
        fixedPointsTrans = transformPointsInverse(tformSML2R,fixedPoints);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % test figure transform with correction (paired localizations)
        figure('Position',visuWinSize);
        hold on;
        scatter(movingPoints(:,1),movingPoints(:,2),'LineWidth',1)
        scatter(fixedPointsTrans(:,1),fixedPointsTrans(:,2),'LineWidth',1)
        for k = 1:size(movingPoints,1)
            plot([movingPoints(k,1) fixedPointsTrans(k,1)],[movingPoints(k,2) fixedPointsTrans(k,2)],'-g');
        end
        axis image ij
        axis(zoomRegion)
        % axis(workingRegion)
        grid on
        grid minor
        legend('rightPoints','leftPointsTrans')
        xlabel('x position [nm]')
        ylabel('y position [nm]')
        title({'Geometric transform of SM images';...
            ['with corrected transform, Paired # = ' num2str(size(movingPoints,1))]})
        hold off
        fig.h9 = gcf;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % test figure transform with correction (histogram)
        figure('Position',visuWinSize);
        hold on;
        histTest1 = histogram(movingPoints(:,1)-fixedPointsTrans(:,1));
        histTest2 = histogram(movingPoints(:,2)-fixedPointsTrans(:,2));
        histTest2.BinWidth = histTest1.BinWidth;
        legend(['xDir, median:' num2str(median(movingPoints(:,1)-fixedPointsTrans(:,1))) 'nm'],...
            ['yDir, median:' num2str(median(movingPoints(:,2)-fixedPointsTrans(:,2))) 'nm'])
        xlabel('Error [nm]')
        % ylabel('Frequency')
        title({'Error distribution of geometric transform';...
            'with corrected transform (FRE)'})
        hold off
        fig.h10 = gcf;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        lengthRegCorr = hypot((movingPoints(:,1)-fixedPointsTrans(:,1)),(movingPoints(:,2)-fixedPointsTrans(:,2)));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % test figure transform with correction (histogram of distance)
        figure('Position',visuWinSize);
        hold on;
        histogram(lengthRegCorr);
        legend(['median= ' num2str(median(lengthRegCorr)) ' nm']);
        xlabel('Error [nm]')
        % ylabel('Frequency')
        title({'Error distribution of geometric transform';...
            'with corrected transform (distance, FRE)'})
        hold off
        fig.h11 = gcf;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end    
end