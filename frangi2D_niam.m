function [outIm,whatScale,Direction] = frangi2D_niam(I, options, angle)
% This function FRANGIFILTER2D uses the eigenvectors of the Hessian to
% compute the likeliness of an image region to vessels, according
% to the method described by Frangi:2001 (Chapter 2).
%
% [J,Scale,Direction] = FrangiFilter2D(I, Options)
%
% inputs,
%   I : The input image (vessel image)
%   Options : Struct with input options,
%       .FrangiScaleRange : The range of sigmas used, default [1 8]
%       .FrangiScaleRatio : Step size between sigmas, default 2
%       .FrangiBetaOne : Frangi correction constant, default 0.5
%       .FrangiBetaTwo : Frangi correction constant, default 15
%       .BlackWhite : Detect black ridges (default) set to true, for
%                       white ridges set to false.
%       .verbose : Show debug information, default true
%
% outputs,
%   J : The vessel enhanced image (pixel is the maximum found in all scales)
%   Scale : Matrix with the scales on which the maximum intensity 
%           of every pixel is found
%   Direction : Matrix with directions (angles) of pixels (from minor eigenvector)   
%
% Example,
%   I=double(imread ('vessel.png'));
%   Ivessel=FrangiFilter2D(I);
%   figure,
%   subplot(1,2,1), imshow(I,[]);
%   subplot(1,2,2), imshow(Ivessel,[0 0.25]);
%
% Written by Marc Schrijver, 2/11/2001
% Re-Written by D.Kroon University of Twente (May 2009)


sigmas=options.FrangiScaleRange(1):options.FrangiScaleRatio:options.FrangiScaleRange(2);
sigmas = sort(sigmas, 'ascend');


% Make matrices to store all filterd images
ALLfiltered=zeros([size(I) length(sigmas)]);
ALLangles=zeros([size(I) length(sigmas)]);

% Frangi filter for all sigmas
for i = 1:length(sigmas),

    % Make 2D hessian
    [Dxx,Dxy,Dyy] = Hessian2D(I,sigmas(i));
    
    % Correct for scale
    Dxx = (sigmas(i)^2)*Dxx;
    Dxy = (sigmas(i)^2)*Dxy;
    Dyy = (sigmas(i)^2)*Dyy;
   
    % Calculate (abs sorted) eigenvalues and vectors
    [Lambda2,Lambda1,Ix,Iy]=eig2image(Dxx,Dxy,Dyy);

    % Compute the direction of the minor eigenvector
    angles = atan2(Ix,Iy);

    % Compute some similarity measures
    Lambda1(Lambda1==0) = eps;
    Rb = (Lambda2./Lambda1).^2; %lambda1 has higher value compared to lambda2, lambda2<lambda1. Maximum at non-tube places where lambda1 approx = lambda2 
    S2 = Lambda1.^2 + Lambda2.^2;
   
    % Compute the output image
    Ifiltered = exp(-Rb) .*(ones(size(I))-exp(-S2));
 
    mask_angle = angle / 180 * pi;  % we remove all lines that are 25 degree from vertical line
    
    Ifiltered(Lambda1>0)=0;
    Ifiltered( ( abs(angles) < mask_angle )  )=0;
    Ifiltered( ( abs(angles) > (pi - mask_angle) ) )=0;
    % store the results in 3D matrices
    ALLfiltered(:,:,i) = normalize(Ifiltered); % niam : added normalization in here, so similar contribution from different scales
    ALLangles(:,:,i) = angles;
end
% imagesc (ALLangles(:,:,1))
% Return for every pixel the value of the scale(sigma) with the maximum 
% output pixel value
if length(sigmas) > 1,
    [outIm,whatScale] = max(ALLfiltered,[],3);
    outIm = reshape(outIm,size(I));
    if(nargout>1)
        whatScale = reshape(whatScale,size(I));
    end
    if(nargout>2)
        Direction = reshape(ALLangles((1:numel(I))'+(whatScale(:)-1)*numel(I)),size(I));
    end
else
    outIm = reshape(ALLfiltered,size(I));
    if(nargout>1)
            whatScale = ones(size(I));
    end
    if(nargout>2)
        Direction = reshape(ALLangles,size(I));
    end
end