% we take only the first 2mm points

function [out] = ray_cast_DDH (in)

im = in;
%im=bone;

[r,c] = size(im);
out = zeros(r,c);
num_points = round( size(im,2)/40 *2); % 2mm

for i=1:c
    tmp = find(im(:,i));
    tmp2 = length(tmp);
    if tmp2>0
        index1 = max(1, tmp2-num_points);
        index2 = tmp2;
        index3 = tmp(index1:index2);
        out(index3,i) = im(index3,i);
    else
    end
    
end
