function [nll, gnll] = neg_log_like_Pois_abrr(c, Ex, Ey, Ez, Z, Ihat, Y, B, pmask, varargin)
%neg_log_like_Pois_abrr computes the negative log likelihood and its gradient
%with respect to  aberration coefficients
%of a set of independent Poisson random variables describing....
%the image signal consists of K frames of images of isolated beads.
%
%->---
%Input
%->---
% c:          array-(n,1)- Zernike coefficients (using Noll index)
% Ex,y,z:     array (N,N)- Electric fields at the pupil plane (X,Y and Z basis)
% Z:          array (N,N,n)- Zernike bases
% Ihat:       array (K,1)- Brightness estimates for emitters in image stack (K frames)
% Y:          array (m,m,K)- Noisy images
% B:          array (m,m,K)- Estimates of local background in  image stack (K frames)
% pmask:      array (N,N)- Phase mask from which images were acquired (in complex form)
%--->-
%output
%--->-
%nll          double- Values of the negative log-likelihood
%gnll:        array (n,1)- Gradients of the negative log-likelihood with
%respect to  aberration coefficients

s = opt2struct(varargin);

m = size(Y, 1); % image size
N = size(pmask, 1);
num_coeff = numel(c);
K = size(Y, 3); % number of frames


% common terms for each frame and each basis, X,Y,Z.
[qx, FourierPupilPlanex, FourierGradientTermx] = ComputeCommonTerms(Ex, Z, pmask, c, varargin{:});
[qy, FourierPupilPlaney, FourierGradientTermy] = ComputeCommonTerms(Ey, Z, pmask, c, varargin{:});
[qz, FourierPupilPlanez, FourierGradientTermz] = ComputeCommonTerms(Ez, Z, pmask, c, varargin{:});


% unnormalized gradients

Gx = 2 * real(bsxfun(@times, FourierGradientTermx, conj(FourierPupilPlanex)));
Gy = 2 * real(bsxfun(@times, FourierGradientTermy, conj(FourierPupilPlaney)));
Gz = 2 * real(bsxfun(@times, FourierGradientTermz, conj(FourierPupilPlanez)));


% normalization
q = (qx + qy + qz);
q = q(-(m - 1)/2+N/2+1:1:(m - 1)/2+N/2+1, -(m - 1)/2+N/2+1:1:(m - 1)/2+N/2+1, :);
sumnorm = sum(sum(q(:, :, 1)));
q = q / sumnorm;

if isfield(s, 'gausssmoothsigmapixel') && s.gausssmoothsigmapixel > 0

    q = imgaussfilt(q, s.gausssmoothsigmapixel);
end


% match to the image size m
Gx = Gx(-(m - 1)/2+N/2+1:1:(m - 1)/2+N/2+1, -(m - 1)/2+N/2+1:1:(m - 1)/2+N/2+1, :, :) / sumnorm;
Gy = Gy(-(m - 1)/2+N/2+1:1:(m - 1)/2+N/2+1, -(m - 1)/2+N/2+1:1:(m - 1)/2+N/2+1, :, :) / sumnorm;
Gz = Gz(-(m - 1)/2+N/2+1:1:(m - 1)/2+N/2+1, -(m - 1)/2+N/2+1:1:(m - 1)/2+N/2+1, :, :) / sumnorm;


Gx = repmat(Gx, 1, 1, K);
Gy = repmat(Gy, 1, 1, K);
Gz = repmat(Gz, 1, 1, K);

GIsoTropic = Gx + Gy + Gz;

if isfield(s, 'gausssmoothsigmapixel') && s.gausssmoothsigmapixel > 0

    GIsoTropic = imgaussfilt(GIsoTropic, s.gausssmoothsigmapixel);
end

q = repmat(reshape(q, m^2, 1), 1, K);

Y = reshape(Y(:, :, 1:K), m^2, K);
B = reshape(B(:, :, 1:K), m^2, K);
Ihat = reshape(Ihat, 1, K);
% Poisson negative log-likelihood

nll = sum(sum(bsxfun(@times, q, Ihat))) - sum(sum(Y .* log(bsxfun(@times, q, Ihat) + B)));
G = reshape(GIsoTropic, m^2, K, num_coeff);

% computing gradients
gnll = reshape(sum(sum(bsxfun(@times, G, Ihat), 2), 1), num_coeff, 1) - ...
    reshape((sum(sum(bsxfun(@times, Y, (bsxfun(@rdivide, bsxfun(@times, G, Ihat), ...
    bsxfun(@times, q, Ihat) + B))), 1), 2)), num_coeff, 1);
end
