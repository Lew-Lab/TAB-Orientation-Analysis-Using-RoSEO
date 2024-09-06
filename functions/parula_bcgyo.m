function map = parula_bcgyo(m)
%PARULA Blue-Cyan-Green-Yellow-Orange color map
%   PARULA_BCGYO(M) returns an M-by-3 matrix containing a colormap. 
%   The colors begin with blue, range through cyan, green and yellow, and 
%   end with bright orange. This is a modified parula colormap.
%
%   PARULA_BCGYO returns a colormap with the same number of colors as the 
%   current figure's colormap. If no figure exists, MATLAB uses the length 
%   of the default colormap.
%
%   EXAMPLE
%
%   This example shows how to reset the colormap of the current figure.
%
%       colormap(parula_bcgyo)
%
%   Modified in Nov. 2017 by Tianben Ding for better rendering of linear 
%   dichroism map.

if nargin < 1
   f = get(groot,'CurrentFigure');
   if isempty(f)
      m = size(get(groot,'DefaultFigureColormap'),1);
   else
      m = size(f.Colormap,1);
   end
end

values = [0	0.498039216	1
0	0.53150326	0.997385621
0	0.564967334	0.994771242
0	0.598431349	0.992156863
0	0.631895423	0.989542484
0	0.665359497	0.986928105
0	0.698823512	0.984313726
0	0.732287586	0.981699347
0	0.76575166	0.979084969
0	0.799215674	0.97647059
0	0.832679749	0.973856211
0	0.866143763	0.971241832
0	0.899607837	0.968627453
0	0.933071911	0.966013074
0	0.966535926	0.963398695
0	1	0.960784316
0	1	0.906617641
0	1	0.852450967
0	1	0.798284292
0	1	0.744117677
0	1	0.689951003
0	1	0.635784328
0	1	0.581617653
0	1	0.527450979
0	1	0.473284304
0	1	0.419117659
0	1	0.364950985
0	1	0.31078431
0	1	0.256617635
0	1	0.202450976
0	1	0.148284316
0	1	0.094117649
0.0625	1	0.088235296
0.125	1	0.082352944
0.1875	1	0.076470591
0.25	1	0.070588239
0.3125	1	0.064705886
0.375	1	0.05882353
0.4375	1	0.052941177
0.5	1	0.047058824
0.5625	1	0.041176472
0.625	1	0.035294119
0.6875	1	0.029411765
0.75	1	0.023529412
0.8125	1	0.01764706
0.875	1	0.011764706
0.9375	1	0.005882353
1	1	0
1	0.962499976	0
1	0.925000012	0
1	0.887499988	0
1	0.850000024	0
1	0.8125	0
1	0.774999976	0
1	0.737500012	0
1	0.699999988	0
1	0.662500024	0
1	0.625	0
1	0.587499976	0
1	0.550000012	0
1	0.512499988	0
1	0.474999994	0
1	0.4375	0
1	0.400000006	0
];

P = size(values,1);
map = interp1(1:size(values,1), values, linspace(1,P,m), 'linear');
