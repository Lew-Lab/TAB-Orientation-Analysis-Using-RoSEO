% Tianben Ding 181016
% 2 dimensional Gaussian function with rotation
% INPUTS
% x: x(1) - amplitude, x(2) - mean in x direction, x(3) - std in x
% direction, x(4) - mean in y direction, x(5) - std in y direction, x(6) -
% clockwise rotation angle
% coord(:,:,1): x coordinates
% coord(:,:,2): y coordinates

% ref: https://en.wikipedia.org/wiki/Gaussian_function

function f = twoDGaussFun_withRot(x,coord)

a = cos(x(6))^2/(2*x(3)^2)+sin(x(6))^2/(2*x(5)^2);
b = -sin(2*x(6))/(4*x(3)^2)+sin(2*x(6))/(4*x(5)^2);
c = sin(x(6))^2/(2*x(3)^2)+cos(x(6))^2/(2*x(5)^2);

f = x(1)*exp(   -( a*(coord(:,:,1)-x(2)).^2 + 2*b*(coord(:,:,1)-x(2)).*(coord(:,:,2)-x(4)) + c*(coord(:,:,2)-x(4)).^2 )    ); 
end