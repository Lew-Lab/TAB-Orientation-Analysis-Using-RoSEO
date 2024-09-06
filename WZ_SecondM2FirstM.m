% WZ 20211020

clear
codepath = 'D:\MATLAB_Codes\TAB SMOLM using RoSEO TD\';
addpath([codepath,'functions']);
addpath(genpath([codepath,'utils']));
addpath([codepath,'phasemask']);

backgFile = ['K:\WZ 1 data\210928_2158BBP LD_Nile Red_ultra pure H2O (NR in 1X PBS)\offsetSubtracted\Data_7\BGEst_lineFitWavelet\4RoSEO' '\croppedData4RoSEO.h5'];
% backgFile = ['K:\WZ 1 data\210902_2158BBP LL_Nile Red_ultra pure H2O (NR in 1X PBS)\offsetSubtracted\Data_6\BGEst_lineFitWavelet\4RoSEO' '\croppedData4RoSEO.h5'];
load('D:\OneDrive - Washington University in St. Louis\code\TDcode\TAB SMOLM using RoSEO TD\211012_estMMatrixRoSEO_v16\210928 Data_17 reg0.25 lineFitWavelet all analyzed results.mat')

%% estimation of physical angle from estimated M matrix
% projection to physical angles
mux = nan(length(M),1);
muy = nan(length(M),1);
rotMobil = nan(length(M),1);
muz = nan(length(M),1);

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


Bs.XX = [XXx XXy];
Bs.YY = [YYx YYy];
Bs.ZZ = [ZZx ZZy];
Bs.XY = [XYx XYy];
Bs.XZ = [XZx XZy];
Bs.YZ = [YZx YZy];
B_Ha.aa = [Bx.aa By.aa];
B_Ha.ab = [Bx.ab By.ab];
B_Ha.ac = [Bx.ac By.ac];
B_Ha.ad = [Bx.ad By.ad];
B_Ha.ae = [Bx.ae By.ae];
B_Ha.af = [Bx.ad By.af];
B_Ha.bb = [Bx.bb By.bb];
B_Ha.bc = [Bx.bc By.bc];
B_Ha.bd = [Bx.bd By.bd];
B_Ha.be = [Bx.be By.be];
B_Ha.bf = [Bx.bf By.bf];
B_Ha.cc = [Bx.cc By.cc];
B_Ha.cd = [Bx.cd By.cd];
B_Ha.ce = [Bx.ce By.ce];
B_Ha.cf = [Bx.cf By.cf];
B_Ha.dd = [Bx.dd By.dd];
B_Ha.de = [Bx.de By.de];
B_Ha.df = [Bx.df By.df];
B_Ha.ee = [Bx.ee By.ee];
B_Ha.ef = [Bx.ef By.ef];
B_Ha.ff = [Bx.ff By.ff];

backg = h5read(backgFile,'/backg');
frameloc = zeros(frameN,1);

for frameidx = 1:frameN
    locDataTmp = loc_data{frameidx};
    backgTemp = backg(:,:,frameidx);
    frameloc(frameidx) = size(loc_data{frameidx},1);
    if frameidx == 1
        framelocidx = 1:frameloc(frameidx);
    else
        framelocidx = sum(frameloc(1:(frameidx-1)))+(1:frameloc(frameidx));
    end
    for idx = 1:length(framelocidx)
        idxidx = framelocidx(idx);
        [mux(idxidx),muy(idxidx),muz(idxidx),rotMobil(idxidx)] = secondM2SymmConeWeighted(Bs,B_Ha,sumNorm,locDataTmp(idx,5:10),locDataTmp(idx,4),backgTemp);
    end
end

% check
MP = rotMobil.*[mux.^2 muy.^2 muz.^2 mux.*muy mux.*muz muy.*muz]+(1-rotMobil)/3.*[1 1 1 0 0 0];