% WZ 20211020 select ROI and reserve the locs of SMs within the ROI

% clear
% load('C:\Users\23116\OneDrive - Washington University in St. Louis\code\TDcode\TAB SMOLM using RoSEO TD\211012_estMMatrixRoSEO_v16\210928 Data_7 reg0.25 lineFitWavelet all analyzed results.mat')

newdatasetFlag = 0;
loc_data_all = nan(length(M),10);
frameloc = zeros(frameN,1);

for frameidx = 1:frameN
    locDataTmp = loc_data{frameidx};
    frameloc(frameidx) = size(loc_data{frameidx},1);
    if frameidx == 1
        framelocidx = 1:frameloc(frameidx);
    else
        framelocidx = sum(frameloc(1:(frameidx-1)))+(1:frameloc(frameidx));
    end
    for idx = 1:length(framelocidx)
        idxidx = framelocidx(idx);
        loc_data_all(idxidx,:) = locDataTmp(idx,:);
    end
end

loc_data_all = [loc_data_all mux muy muz rotMobil];

xq = loc_data_all(:,2)/pixelSize;
yq = loc_data_all(:,3)/pixelSize;
figure; 
plot(xq,-yq,'.'); % "-" is add to match the picture returned by RoSE-O.
axis image;

if newdatasetFlag
    roi = drawpolygon;    % interactive, draw polygon according to SR image
    % save the roi data e.g., save 210928_Data_10_ROI.mat roi
else
%     load 210928_Data_10_ROI.mat
end

%%
xv = roi.Position(:,1);
yv = -roi.Position(:,2);

[in,on] = inpolygon(xq,yq,xv,yv);

%%%%%% check the separation of SMs inside and outside of ROI
figure; plot(xv,yv); axis equal; hold on
plot(xq(in),yq(in),'r.'); plot(xq(~in),yq(~in),'b.'); hold off

figure; histogram(mux(in)); title('\mu_x ROI')
figure; histogram(muy(in)); title('\mu_y ROI')
figure; histogram(muz(in)); title('\mu_z ROI')

theta = nan(length(M),1); 
muzflag = (muz>0);
theta(muzflag) = acos(muz(muzflag));
theta(~muzflag) = acos(-muz(~muzflag));

phi = nan(length(M),1);
phi(muzflag) = atan2(muy(muzflag),mux(muzflag));
phi(~muzflag) = atan2(muy(~muzflag),mux(~muzflag));

alpha = acos( (  -1 + sqrt(1+8*rotMobil)  )/2);
omega = 2.*4.*pi.*sin(alpha/2).^2;
figure;
subplot(2,2,1)
histogram(rad2deg(alpha));
xlabel('[deg]','FontSize',18)
title('\alpha')
subplot(2,2,2)
histogram(rotMobil);
title('\gamma')
subplot(2,2,3)
histogram(rad2deg(theta));
xlabel('[deg]','FontSize',18)
title('\theta')
subplot(2,2,4)
histogram(rad2deg(phi));
xlabel('[deg]','FontSize',18)
title('\phi')
loc_data_roi = loc_data_all(in,:); 
% figure; plot(loc_data_roi(:,2),loc_data_roi(:,3),'.') % check ROI selection
thetaROI = theta(in);
alphaROI = alpha(in);
phiROI = phi(in);

figure;
h = binscatter(rad2deg(thetaROI),rad2deg(alphaROI),[90 90]);
xlabel('\theta [deg]')
ylabel('\alpha [deg]')
axis equal
xlim([0 90]); ylim([0 90])
colormap(gca,'parula')
% caxis([0,10])









