% 190508 - v2 TD - Added empty minPara and maxPara options. Also added a 
% mean output plot. Also added a std output plot

% Tianben Ding 190228
% Generate a SR reconstruction by encoding localization number and an
% additional parameter to value and hue in hsv colorscale respectively.

function [rgbPlotMedian,rgbPlotMean,rgbPlotStd,colorBarMedian,colorBarMean,colorBarStd,...
    minParaMedian,maxParaMedian,minParaMean,maxParaMean,minParaStd,maxParaStd] =...
    reconstVLocNumHAddPara_v2(x,y,para,minPara,maxPara,binSize,Xedges,Yedges,maxLocNum,satLocNum,colorMap)
% x: a column vector includes localized x positions
% y: a column vector includes localized y positions
% para: a column vector includes additional parameter values, should have
% minPara: minimum value for Para in the final visualization, only for
% median and mean plots
% maxPara: maximum value for Para in the final visualization, only for
% median and mean plots
% the same size as x and y
% Xedges,Yedges: edge information for 2D histogram
% maxLocNum: maximum localiztion number out of bins, you can get this value
% from a standard SR recontruction with the same edges
% satLocNum: represents the saturation threshold of localization number
% colorMap: colormap, e.g., parula

binNum=zeros(2,1);
binNum(1) = numel(Xedges)-1;
binNum(2) = numel(Yedges)-1;

indAng = zeros(binNum(2),binNum(1));
allAng = nan(binNum(2),binNum(1),maxLocNum);
for ii = 1:length(para)
    xInd = floor((x(ii)-Xedges(1))/binSize)+1;
    yInd = floor((y(ii)-Yedges(1))/binSize)+1;
    indAng(yInd,xInd) = indAng(yInd,xInd) + 1;
    allAng(yInd,xInd,indAng(yInd,xInd)) = para(ii);
end
medianAng = nanmedian(allAng,3);
meanAng = nanmean(allAng,3);
stdAng = nanstd(allAng,0,3);

if isempty(minPara)
    medianAngRef = medianAng;
    medianAngRef(indAng < satLocNum/2) = nan;
    minParaMedian = nanmin(nanmin(medianAngRef));
    meanAngRef = meanAng;
    meanAngRef(indAng < satLocNum/2) = nan;
    minParaMean = nanmin(nanmin(meanAngRef));
else
    minParaMedian = minPara;
    minParaMean = minPara;
end
if isempty(maxPara)
    medianAngRef = medianAng;
    medianAngRef(indAng < satLocNum/2) = nan;
    maxParaMedian = nanmax(nanmax(medianAngRef));
    meanAngRef = meanAng;
    meanAngRef(indAng < satLocNum/2) = nan;
    maxParaMean = nanmax(nanmax(meanAngRef));
else
    maxParaMedian = maxPara;
    maxParaMean = maxPara;        
end
stdAngRef = stdAng;
stdAngRef(indAng < satLocNum/2) = nan;
maxParaStd = nanmax(nanmax(stdAngRef));
minParaStd = nanmin(nanmin(stdAngRef));

medianAng(isnan(medianAng)) = 0; % these pixels will appeal as black pixels in final LD map due to zero "value"
meanAng(isnan(meanAng)) = 0;
stdAng(isnan(stdAng)) = 0;

% output for median
hueInd = medianAng - minParaMedian;
hueInd = round(hueInd*(size(colorMap,1)-1)/(maxParaMedian-minParaMedian)+1);
hueInd( hueInd < 1 ) = 1;
hueInd( hueInd > size(colorMap,1) ) = size(colorMap,1);
colorRGB=reshape(reshape(colorMap,1,size(colorMap,1)*3) , [1,size(colorMap,1),3]);
colorHSV=rgb2hsv(colorRGB);
colorHue=colorHSV(:,:,1);
colorSat=colorHSV(:,:,2);
colorVal=colorHSV(:,:,3);
huePlot = colorHue(hueInd);
saturationPlot = colorSat(hueInd);
valuePlot = colorVal(hueInd);
indAng(indAng > satLocNum) = satLocNum;
valuePlot = valuePlot.*(indAng/satLocNum);
hsvPlot = cat(3,huePlot,saturationPlot,valuePlot);
rgbPlotMedian = hsv2rgb(hsvPlot);

colorHSV1 = repmat(flipud(colorHSV(:,:,1).'),[1,satLocNum+1]);
colorHSV2 = repmat(flipud(colorHSV(:,:,2).'),[1,satLocNum+1]);
colorHSV3 = nan(size(colorHSV1));
for indn = 1:satLocNum+1
    colorHSV3(:,indn) = (indn-1)/satLocNum*flipud(colorHSV(:,:,3).');
end
colorHSVBar = cat(3,colorHSV1.',colorHSV2.',colorHSV3.');
colorBarMedian = hsv2rgb(colorHSVBar);

% output for mean
hueInd = meanAng - minParaMean;
hueInd = round(hueInd*(size(colorMap,1)-1)/(maxParaMean-minParaMean)+1);
hueInd( hueInd < 1 ) = 1;
hueInd( hueInd > size(colorMap,1) ) = size(colorMap,1);
colorRGB=reshape(reshape(colorMap,1,size(colorMap,1)*3) , [1,size(colorMap,1),3]);
colorHSV=rgb2hsv(colorRGB);
colorHue=colorHSV(:,:,1);
colorSat=colorHSV(:,:,2);
colorVal=colorHSV(:,:,3);
huePlot = colorHue(hueInd);
saturationPlot = colorSat(hueInd);
valuePlot = colorVal(hueInd);
indAng(indAng > satLocNum) = satLocNum;
valuePlot = valuePlot.*(indAng/satLocNum);
hsvPlot = cat(3,huePlot,saturationPlot,valuePlot);
rgbPlotMean = hsv2rgb(hsvPlot);

colorHSV1 = repmat(flipud(colorHSV(:,:,1).'),[1,satLocNum+1]);
colorHSV2 = repmat(flipud(colorHSV(:,:,2).'),[1,satLocNum+1]);
colorHSV3 = nan(size(colorHSV1));
for indn = 1:satLocNum+1
    colorHSV3(:,indn) = (indn-1)/satLocNum*flipud(colorHSV(:,:,3).');
end
colorHSVBar = cat(3,colorHSV1.',colorHSV2.',colorHSV3.');
colorBarMean = hsv2rgb(colorHSVBar);

% output for std
hueInd = stdAng - minParaStd;
hueInd = round(hueInd*(size(colorMap,1)-1)/(maxParaStd-minParaStd)+1);
hueInd( hueInd < 1 ) = 1;
hueInd( hueInd > size(colorMap,1) ) = size(colorMap,1);
colorRGB=reshape(reshape(colorMap,1,size(colorMap,1)*3) , [1,size(colorMap,1),3]);
colorHSV=rgb2hsv(colorRGB);
colorHue=colorHSV(:,:,1);
colorSat=colorHSV(:,:,2);
colorVal=colorHSV(:,:,3);
huePlot = colorHue(hueInd);
saturationPlot = colorSat(hueInd);
valuePlot = colorVal(hueInd);
indAng(indAng > satLocNum) = satLocNum;
valuePlot = valuePlot.*(indAng/satLocNum);
hsvPlot = cat(3,huePlot,saturationPlot,valuePlot);
rgbPlotStd = hsv2rgb(hsvPlot);

colorHSV1 = repmat(flipud(colorHSV(:,:,1).'),[1,satLocNum+1]);
colorHSV2 = repmat(flipud(colorHSV(:,:,2).'),[1,satLocNum+1]);
colorHSV3 = nan(size(colorHSV1));
for indn = 1:satLocNum+1
    colorHSV3(:,indn) = (indn-1)/satLocNum*flipud(colorHSV(:,:,3).');
end
colorHSVBar = cat(3,colorHSV1.',colorHSV2.',colorHSV3.');
colorBarStd = hsv2rgb(colorHSVBar);
end