% 201105 - TD - Calculate rmse of compensated x and y channel aberration
% against a phase mask

function [xrms]=rmse_abrr_2_cha(c,Z,xAbrr,yAbrr,pmask,radius,ind)

compPM = sum(bsxfun(@times, Z, reshape(c, 1, 1, size(c, 1))), 3);
compPM = compPM( (size(compPM,1)/2-radius+1) : (size(compPM,1)/2+radius),...
    (size(compPM,1)/2-radius+1) : (size(compPM,1)/2+radius) );
compPM(~ind) = nan;

xn = [reshape( (xAbrr+compPM) ,[numel(xAbrr),1]); reshape( (yAbrr+compPM) ,[numel(yAbrr),1])];
xn(isnan(xn)) = [];
tn = [reshape(pmask,[numel(pmask),1]); reshape(pmask,[numel(pmask),1])];
tn(isnan(tn)) = [];

xrms = sqrt(1/length(xn)*sum( (xn-tn).^2 ));

gxrms = 1./(length(xn).*xrms) .* sum( sum(Z,1),2 );
gxrms = reshape( gxrms,[size(c,1),size(c,2)] );

end




