% 201115 - pix_wei TD - Weighted by BFP intensity

% 201113 - pix TD - Pixel level optimization instead of Zernike
% polynomials.

% 201105 - TD - Calculate rmse of compensated x and y channel aberration
% against a phase mask

function [xrms]=rmse_abrr_2_cha_pix_wei(c,xAbrr,yAbrr,pmask,ind, Ex, Ey)

compPM = c;
compPM(~ind) = nan;

Ex(~ind) = nan;
Ey(~ind) = nan;
E = [reshape( Ex, [numel(Ex),1]); reshape( Ey, [numel(Ey),1]) ];
E(isnan(E)) = [];

xn = [reshape( (xAbrr+compPM) ,[numel(xAbrr),1]); reshape( (yAbrr+compPM) ,[numel(yAbrr),1])];
xn(isnan(xn)) = [];
tn = [reshape(pmask,[numel(pmask),1]); reshape(pmask,[numel(pmask),1])];
tn(isnan(tn)) = [];

xrms = sqrt(1/length(xn)*sum( E.*(xn-tn).^2 ));

% gxrms = 1./(length(xn).*xrms) .* sum( sum(Z,1),2 );
% gxrms = reshape( gxrms,[size(c,1),size(c,2)] );

end




