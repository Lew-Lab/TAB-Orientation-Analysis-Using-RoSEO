function map = cold2hot(m)
%   COLD2HOT(M) returns an M-by-3 matrix containing a colormap. 
%   The colors begin with dark blue, range through light blue, gray, light orange, and 
%   end with dark red.
%


if nargin < 1
   f = get(groot,'CurrentFigure');
   if isempty(f)
      m = size(get(groot,'DefaultFigureColormap'),1);
   else
      m = size(f.Colormap,1);
   end
end
T = [76, 40, 130     %// violet
    49, 104, 174       %// light blue
     176, 180, 185        %// gray
     230, 148, 96       %// light orange
    255, 215, 4]./255; %// gold

x = 1:size(T,1);
 
 map = interp1(x,T,linspace(1,size(T,1),m));

% P = size(values,1);
% map = interp1(1:size(values,1), values, linspace(1,P,m), 'linear');
