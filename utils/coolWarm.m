function map = coolWarm(m)
%coolWarm color map
%   COOLWARM(M) returns an M-by-3 matrix containing a colormap. 
%   The colors begin with blue, range through white, and end with red. 
%
%   COOLWARM returns a colormap with the same number of colors as the current
%   figure's colormap. If no figure exists, MATLAB uses the length of the
%   default colormap.
%
%   EXAMPLE
%
%   This example shows how to reset the colormap of the current figure.
%
%       colormap(coolWarm)
%
%   Modified in Oct. 2019 by Tianben Ding to better render two-end results
%   like bias characterization.

if nargin < 1
   f = get(groot,'CurrentFigure');
   if isempty(f)
      m = size(get(groot,'DefaultFigureColormap'),1);
   else
      m = size(f.Colormap,1);
   end
end

values = [
     0.2275    0.2980    0.7529
    0.2472    0.3157    0.7565
    0.2669    0.3333    0.7600
    0.2866    0.3510    0.7636
    0.3064    0.3686    0.7672
    0.3261    0.3863    0.7707
    0.3458    0.4039    0.7743
    0.3656    0.4216    0.7778
    0.3853    0.4392    0.7814
    0.4050    0.4569    0.7849
    0.4248    0.4745    0.7885
    0.4445    0.4922    0.7920
    0.4642    0.5098    0.7956
    0.4839    0.5275    0.7991
    0.5037    0.5451    0.8027
    0.5234    0.5627    0.8063
    0.5431    0.5804    0.8098
    0.5629    0.5980    0.8134
    0.5826    0.6157    0.8169
    0.6023    0.6333    0.8205
    0.6221    0.6510    0.8240
    0.6418    0.6686    0.8276
    0.6615    0.6863    0.8311
    0.6813    0.7039    0.8347
    0.7010    0.7216    0.8382
    0.7207    0.7392    0.8418
    0.7404    0.7569    0.8453
    0.7602    0.7745    0.8489
    0.7799    0.7922    0.8525
    0.7996    0.8098    0.8560
    0.8194    0.8275    0.8596
    0.8391    0.8451    0.8631
    0.8588    0.8627    0.8667
    0.8543    0.8372    0.8452
    0.8497    0.8116    0.8237
    0.8452    0.7861    0.8022
    0.8406    0.7605    0.7806
    0.8361    0.7350    0.7591
    0.8315    0.7094    0.7376
    0.8269    0.6839    0.7161
    0.8224    0.6583    0.6946
    0.8178    0.6328    0.6731
    0.8133    0.6072    0.6516
    0.8087    0.5817    0.6301
    0.8042    0.5561    0.6086
    0.7996    0.5306    0.5871
    0.7951    0.5050    0.5656
    0.7905    0.4794    0.5441
    0.7860    0.4539    0.5226
    0.7814    0.4283    0.5011
    0.7769    0.4028    0.4796
    0.7723    0.3772    0.4581
    0.7677    0.3517    0.4366
    0.7632    0.3261    0.4151
    0.7586    0.3006    0.3935
    0.7541    0.2750    0.3720
    0.7495    0.2495    0.3505
    0.7450    0.2239    0.3290
    0.7404    0.1984    0.3075
    0.7359    0.1728    0.2860
    0.7313    0.1472    0.2645
    0.7268    0.1217    0.2430
    0.7222    0.0961    0.2215
    0.7176    0.0706    0.2000
];

P = size(values,1);
map = interp1(1:size(values,1), values, linspace(1,P,m), 'linear');
