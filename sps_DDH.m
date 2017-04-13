function [sps1,sps2] = sps_DDH(im, angle)

OPS= phasesymmono_niam(im, 4, 5, 2.1, 0.15, 2.0, 1, -1); 
options.FrangiScaleRatio=2; 
options.FrangiScaleRange=[2 4];
%[Ivessel1, scale, direction]=frangi2D_niam(OPS, options, 0); 
Ivessel1 = frangi2D_niam(OPS, options, 0); 
sps1 = normalize (Ivessel1);
options.FrangiScaleRange=[2 4];
sps2=normalize(frangi2D_niam(OPS, options, angle));

end