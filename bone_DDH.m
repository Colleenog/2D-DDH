function [bone, probability, map] = bone_DDH (im, sps1, sps2)

%im = img;
[r, c] = size(im);


% sps correction
PS_correction = ones (r,c);
for j= 1:r
    PS_correction (j,:) = log (j);
end
sps_new = sps1 .* PS_correction;
sps_new = normalize (sps_new);





%prepare im for confidence map extraction
prep_im = prep_conf_map (im, sps_new);

% Confidence map and feature extraction
[ map, karamalis_bone, combined_attenuation] = bone_feature_extraction(prep_im);

% Search space for bone detection constrained using image-map bone
% indicator- karamalis thesis page 84
search_1_bone = sps_new .* (karamalis_bone>0.0);

% probablity of bone surface
probability = normalize( (sps_new + karamalis_bone + 2* combined_attenuation) .* (search_1_bone>0) );
%probability = normalize( (sps1 + karamalis_bone + combined_attenuation));


% Search space for bone detection constrained using confidence map, mean
% and 2*std, accounts for about 95% of the population
[map_min, map_max] = map_min_max (probability, map);
search_2_map = map .* (map>map_min) .* (map<map_max);
search_2_bone = search_1_bone .* (search_2_map>0);

% Enhancement of bone based on probability 
%[note: this is for visualization, whereas probability is used for classification
bone = search_2_bone .* probability;


probability = normalize( (sps_new + karamalis_bone + 2*combined_attenuation)) ; 
end