% 190719 - TD - SVD only

% 190718 - v7 TD - Use calFIMSecondM_v3 instead of calFIMSecondM_v2 for
% fast calculation

% 190716 - v6 TD - Run fmincon twice with two different initial points

% 190716 - v5 TD - Use Nanoscope object directly. Use calFIMSecondM_v2.

% 190701 - TD - Excluded mux flipping from the code.

% 190701 - TD - changed LS algorithm in from interior-point to sqp

% 190628 - TD - Included the mux sign flipping into symmCone2secM(z)

% 190618 - v4 TD - Start estimation of mux, muy, muz and rotMobil from the
% nearest backbone orientations.

% 190617 - v3 TD - Flip the estimated mux for matching coordinate system of
% objective

% 190404 Tianben Ding
% Maps the second moments onto parameters of a symmetric cone model for a
% rotating molecule

% Weight LS estimation using FIM
function [mux,muy,muz,rotMobil,x0]= secondM2SymmConeWeighted_v7(bx,by,Bx,By,sumNorm,secM,signal,backg)
%% define a least square problem structure

Ix = secM(1).*(bx.XX) + secM(2).*(bx.YY) + secM(3).*(bx.ZZ) + ...
    secM(4).*(bx.XY) + secM(5).*(bx.XZ) + secM(6).*(bx.YZ);

Iy = secM(1).*(by.XX) + secM(2).*(by.YY) + secM(3).*(by.ZZ) + ...
    secM(4).*(by.XY) + secM(5).*(by.XZ) + secM(6).*(by.YZ);

Ix = Ix/sumNorm;
Iy = Iy/sumNorm;

Ix = signal.*Ix;
Iy = signal.*Iy;

FIM = calFIMSecondM_v3(Bx,By,Ix,Iy,signal,backg);

%objective function
%------------------------------------------------------------
objFcn=@(z)(symmCone2secM(z).'-secM)*FIM*(symmCone2secM(z)-secM.');

%  constraints in full 3D coordinate

lb=[-ones(3,1);zeros(1,1)];
ub=ones(4,1);

% options=optimoptions('fmincon','Display','none');
options = optimoptions('fmincon','Display','none','Algorithm','interior-point');

% % Initial value based on a prior knowledge
% x0Prior(1) = locPosROIRefMu(1);
% x0Prior(2) = locPosROIRefMu(2);
% x0Prior(3) = locPosROIRefMu(3);
% x0Prior(4) = 1;
% 
% [solPrior,~]=fmincon(objFcn,x0Prior,[],[],[],[],lb,ub,@mycon,options);

% Initial value based SVD
%------------------------------------------------------------
% construct the M matrix
M=[secM(1),secM(4),secM(5);....
    secM(4),secM(2),secM(6);...
    secM(5),secM(6),secM(3)];
[V,D]=eig(M);

% initialization via SVD
x0SVD(1)=real(V(1,3));
x0SVD(2)=real(V(2,3));
x0SVD(3)=real(V(3,3));
x0SVD(4)=1.5*real(D(3,3))-.5;

[solSVD,~]=fmincon(objFcn,x0SVD,[],[],[],[],lb,ub,@mycon,options);

% if fvalSVD <= fvalPrior
    mux=solSVD(1);
    muy=solSVD(2);
    muz=solSVD(3);
    rotMobil=solSVD(4);
    x0 = x0SVD;
% else
%     mux = solPrior(1);
%     muy = solPrior(2);
%     muz = solPrior(3);
%     rotMobil = solPrior(4);
%     x0 = x0Prior;
% end

% mux = -mux;

%% local functions

%     function out=symmCone2secM(z)
%
%             z=reshape(z,1,3);
%             %muz=sqrt(1-(z(:,1).^2+z(:,2).^2));
%             muxx=z(:,3).*z(:,1).^2+(1-z(:,3))/3;
%             muyy=z(:,3).*z(:,2).^2+(1-z(:,3))/3;
%             muzz=z(:,3).*muz.^2+(1-z(:,3))/3;
%             muxy=z(:,3).*z(:,1).*z(:,2);
%             muxz=z(:,3).*z(:,1).*muz;
%             muyz=z(:,3).*z(:,2).*muz;
%
%             out=[muxx,muyy,muzz,muxy,muxz,muyz]';
%     end

    function out=symmCone2secM(z)
        z=reshape(z,1,4);
        %             z(:,1) = -z(:,1); % flip x coordinate to match with imaging system coordinates
        muz_t=z(:,3);
        muxx=z(:,4).*z(:,1).^2+(1-z(:,4))/3;
        muyy=z(:,4).*z(:,2).^2+(1-z(:,4))/3;
        muzz=z(:,4).*muz_t.^2+(1-z(:,4))/3;
        muxy=z(:,4).*z(:,1).*z(:,2);
        muxz=z(:,4).*z(:,1).*muz_t;
        muyz=z(:,4).*z(:,2).*muz_t;
        out=[muxx,muyy,muzz,muxy,muxz,muyz]';
    end

    function [c,ceq] = mycon(z)
        c=[];
        ceq=sum(z(:,1:3).^2,2)-1;
    end

end



