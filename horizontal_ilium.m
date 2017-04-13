function [hor_pared_ilium, hor_length, ver_length, major_angle] = horizontal_ilium (ilium, size_mm) 

% Version 1: only left half of US image is analyzed. 
% assumed ilium in cranial direction, and FH and labrum and cartilage in
% caudal


% hor_pared_ilium
% ilium perpendicular to ultrasound beam, convert to bw

[rows, columns] = size(ilium); 
ilium_bw = zeros(rows, columns);
ilium_bw(round(rows/3): round(rows/4*3), :) = (ilium(round(rows/3): round(rows/4*3), :)>0);




% Pare it down to single pixel strands
imsk(:,1:columns/2) = bwmorph( ilium_bw(:,1:columns/2),'thin',Inf);

% Obtain the length (count of pixels for each segment, i.e the Area)
lineStats = regionprops(imsk, {'Area','PixelIdxList'});
[l, index] = max([lineStats.Area]); % l is the number of non-zeros in the maximum line

% Create a new image with just the longest line
if isempty(index)
    longestLine = zeros(size(imsk));
else
    longestLine = zeros(size(imsk));
    longestLine(lineStats(index).PixelIdxList)=1;
end

imsk2 = zeros(rows, columns);
imsk2(:,1:columns/2) = longestLine;


if nnz(imsk2)==0
    imsk2(floor(rows/2),floor(columns/2))=1; 
    imsk2(floor(rows/2),floor(columns/2)+1)=1; 
else
end
%imagesc(ilium)
hor_pared_ilium = imsk2;
%imagesc(imsk2)

%hor_length, ver_length
temp = sum(imsk2);
temp1 = max(find(temp));
temp2 = min(find(temp)); 
hor_length = (temp1-temp2)/size_mm; %hor_length is the horizontal length of the longest line

temp3 = min( find( imsk2(:,temp1) ));
temp4 = min( find( imsk2(:,temp2) ));
ver_length = abs(temp3-temp4)/size_mm;


    major_angle = atan(ver_length/hor_length) *180/pi;


% quality_ilium_binary = (hor_length > .33*columns/2) .* (major_angle < 10);
% quality_ilium_horizontality = (10-major_angle).*quality_ilium_binary; %value from 0 to 10

end
