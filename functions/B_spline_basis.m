% Tianben Ding 181017
% B-spline function
% INPUTS
% t: coordinate
% q: order of B-spline basis function

function BOut = B_spline_basis(t,q)
if q == 1
    % BOut = size(t);
    ifSmall = t<1;
    ifLarge = t>=0;
    ifWithin = (double(ifSmall) + double(ifLarge)) == 2;
    
    BOut = double(ifWithin);
else
    BOut = t/(q-1).*B_spline_basis(t,q-1) + (q-t)/(q-1).*B_spline_basis(t-1,q-1);
end


end