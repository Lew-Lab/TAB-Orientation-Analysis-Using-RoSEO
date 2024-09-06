function [img, FourierPupilPlane, FourierGradientTerm] = ComputeCommonTerms(Ey, Z, pmask, c, varargin)

% computing gradients with respect to Zernike coefficients
% and image formed on the camera

%inputs:

%E:    electric feild at pupil plane
%Z:     Zernike basis
%pmask: phase mask mounted on SLM
%c:     Zenike mode (Null index) coefficients

%output:

%img:            image formed on camera
%FourierPlane:   Fourier transform of the feild at pupil plane
%FourierGradientTerm: fourier transform of the gradient of the pupil with
%respect to coefficients
s = opt2struct(varargin);


N = size(pmask, 1);

if size(c, 1) == 1

    c = c';
end

if isempty(Z)

    PhaseAberr = ones(N);
    Z = zeros(N);
else

    PhaseAberr = exp(1i*sum(bsxfun(@times, Z, reshape(c, 1, 1, size(c, 1))), 3));
end


FourierPupilPlane = ifftshift(fftshift(fft2(bsxfun(@times, Ey, pmask .* PhaseAberr))), 3);


% computing the terms related to the gradients

FourierGradientTerm = ifftshift(ifftshift(fftshift(fft2(bsxfun(@times, bsxfun(@times, Ey, pmask .* PhaseAberr), ...
    1i * reshape(Z, N, N, 1, length(c))))), 3), 4);

% apply transform due to channels
if isfield(s, 'xchannel') && s.xchannel

    FourierPupilPlane = permute(FourierPupilPlane, ...
        [2, 1]);

    FourierPupilPlane = rot90(FourierPupilPlane, 2);

    FourierGradientTerm = permute(FourierGradientTerm, ...
        [2, 1, 3, 4]);

    FourierGradientTerm = rot90(FourierGradientTerm, 2);

elseif isfield(s, 'ychannel') && s.ychannel

    FourierPupilPlane = rot90(FourierPupilPlane, 2);

    FourierGradientTerm = rot90(FourierGradientTerm, 2);
end

img = abs(FourierPupilPlane).^2; % image formed on the camera

end