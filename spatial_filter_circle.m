function filtered = spatial_filter_circle(im, radius, center_x, center_y)

imageSizeX = size(im,2);
imageSizeY = size(im,1);
[columnsInImage rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
% Next create the circle in the image.

%center_x = 120;
%center_y = 100;
%radius = 50;

center_y = round(center_y);
center_x = round(center_x);

circlePixels = (rowsInImage - center_y).^2 ...
    + (columnsInImage - center_x).^2 <= radius.^2;


filtered = im .*  circlePixels;
