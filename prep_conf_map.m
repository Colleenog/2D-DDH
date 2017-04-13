function prep_im = prep_conf_map (im, ps)


% % prepare im for confidence map estimation
% Gaussian filtering. Note: structured anisotropic filtering for future

scale = round( size(im,2)/40 *3/2); % size of kernel
deviation = round (scale/6);

h = fspecial('gaussian', [scale scale], deviation);
im2 = imfilter(im, h, 'replicate');
ps = imfilter(ps, h, 'replicate');

% %% more systematic version
% [r, c] =size(im);
% % calculate weight of added structure, by correlation coefficient C
% ps = ps .* (ps>0);
% temp_Ivessel = reshape(sps1,[r*c 1]);
% [r2,c2,temp_Ivessel1] = find(temp_Ivessel);
% 
% temp_img = (sps1>0) .* im2; 
% temp_img2 = reshape(temp_img,[r*c 1]);
% [r1,c1,temp_img3] = find(temp_img2) ;
% temp_Ivessel2 = temp_Ivessel1 (1: length(temp_img3) );
% 
% % scale, to estimate the maximum of sps to be added to im
% scale_prep = max(temp_Ivessel2) ./ max(temp_Ivessel1);
% 
% C = corr2 (temp_img3, temp_Ivessel2);
% prep_im =   normalize (im2) + (1-C)/2 .* scale_prep .* normalize (ps); % prep_im is prepared image
% 
% prep_im = normalize (prep_im);

% %% shorter version
% prep_im =   normalize (im2) + normalize (ps);
% prep_im = normalize (prep_im);

%% shortest version % seems to give the best result
prep_im = im2;