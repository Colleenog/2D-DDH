function [center_FH_new, radius_FH_new, quality_FH, quality_trirad] = femoral_head_feb_miccai (sps1,probability,map,radius_FH,center_FH,center_alpha_beta,alpha,beta,size_mm,im) 



[r,c]=size(map);
 % extract confidence map based search area
 


    map_tmp = spatial_filter_circle(map, radius_FH, center_FH(1), center_FH(2));
   % map_tmp(:, round(center_FH(1)):end) = 0;% keep only left half, right half has hip
    map_tmp1 = map_tmp;
    map_tmp1 (r,c) = 1; % to keep colormap same
    %figure, imagesc(map_tmp1)
         
%    if (max(map_tmp1)>0)
%     else
%         center_FH_new=0;
%         radius_FH_new=0;
%         quality_FH=0;
%         quality_trirad=0;
%         return
%    end

   
    temporary = unique(map_tmp(:)); 
    if temporary(1,1)>0
        map_min = temporary(2); %map_min
    else
        map_min = 0.05;
    end
    map_max = max(max(map_tmp));
    
    
    map_BW = (map>=(map_min)).*(map<=map_max);
  %  figure, imagesc(map_BW), colormap gray, axis equal     


     
  %  figure, imagesc(map_tmp), axis equal     
    
    
%sps and gradient
    sps1_map = sps1 .* map_BW .* (sps1>0.1) ; 
    imagesc(sps1_map)
    % sps1_map=sps1; % 18 feb, try without map_BW
    sps1_map_cartilage = sps1_map .* (probability>.5);
    sps1_map_cartilage(:,1:center_alpha_beta(1,1))=0;
    img_temp = (sps1_map_cartilage>0); 
    tmp = round (center_FH(1,2));
    img4 = img_temp(1:tmp,:);
    h = [1;-1];
    img5 = imfilter(img4,h); 

    img4 = img_temp(tmp+1:end,:);
    h = [-1;1];
    img6 = imfilter(img4,h);
    img7 = [img5;img6]; 
    
    
    img7(:,1:center_alpha_beta(1,1))=0;
    
% % ROI
%     index1 = min(center_alpha_beta(1,1)+size_mm*30, c);
%     sps1_map_zoomed = img7 (round(r*1/5):round(r*4/5), center_alpha_beta(1,1):index1); % horizontal diameter max is 18mm, we take upto 22mm
%          % to remove the greater trochanter, just remove the sps that have vertical eigenvectors... do later
%     img7 = sps1_map_zoomed;
%     figure, imagesc(img7), colormap gray, axis equal
   
% hough based circle fit
 
    
    %18 feb attempt
        radii = size_mm*6:1:size_mm*9;
    h = circle_hough(img7, radii, 'same');  h1=h;
    
    [h_r, h_c]=size(h(:,:,1));
    region_bw_default = ones(h_r,h_c);
    for i = 1:size(h,3)
        region_bw=region_bw_default;
        region_bw(:,1:center_alpha_beta(1,1)+radii(i)+size_mm)=0;
        region_bw(:,center_alpha_beta(1,1)+round(radii(i)*1.5+size_mm):end)=0;
        region_bw(1:center_alpha_beta(1,2)-round(radii(i)*1.5+size_mm),:)=0;
        region_bw(center_alpha_beta(1,2)+round(radii(i)*1.5+size_mm):end,:)=0;
        region_bw = imresize(region_bw, [h_r, h_c]);
        h_new(:,:,i) = h(:,:,i) .* region_bw;
        %h_new enforces the circle to be close to the initial estimate obtained from angle locations/measurements
    end
    h = h_new;
    peaks = circle_houghpeaks(h, radii, 'nhoodxy', 15, 'nhoodr', 21, 'npeaks', 1);
    
%     h_number = peaks(3,1)-size_mm*6+1;
%     figure, imagesc(h(:,:,h_number)), colormap gray, axis equal
    

    
    
    
% check for hough circles        
    tanA = tan(alpha*pi/180);
    tanB = tan(beta*pi/180);

    FHC_check = tanA / (tanA+tanB);


           
    %hold on;
    for i = 1:size(peaks,2)    
        center = [peaks(1,i), peaks(2,i)];
        radius = peaks(3,i);
       % viscircles(center, radius,'EdgeColor','b');
        FHC_1 = center(1,2)+radius+round(r*1/5)-center_alpha_beta(1,2);
        FHC_2 = radius * 2;
        FHC = FHC_1/FHC_2;

        if abs(FHC-FHC_check)/FHC_check<0.4
            
        else
        end   

    end
    %hold off
    
% validation femoral head http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2903126/ Correlation of femoral head coverage and Graf ? angle in infants being screened for developmental dysplasia of the hip    
A = exist('center');
if (A)
else
    center_FH_new=0;
    radius_FH_new=0;
    quality_FH=0;
    quality_trirad=0;
    return
end

% return values
    center_FH_new = [center(1,1), center(1,2)];
    radius_FH_new = radius;
    
    
    quality_FH=zeros(1,7);
    map_tmp = spatial_filter_circle(map, radius_FH_new/4, center_FH_new(1), center_FH_new(2));
    quality_FH(:,1)=max(max(map_tmp));
    
    map_tmp = spatial_filter_circle(map, radius_FH_new/3, center_FH_new(1), center_FH_new(2));
    quality_FH(:,2)=max(max(map_tmp));
    
    map_tmp = spatial_filter_circle(map, radius_FH_new/2, center_FH_new(1), center_FH_new(2));
    quality_FH(:,3)=max(max(map_tmp));
    
    map_tmp = spatial_filter_circle(map, radius_FH_new/1.5, center_FH_new(1), center_FH_new(2));
    quality_FH(:,4)=max(max(map_tmp));
    
     map_tmp = spatial_filter_circle(map, radius_FH_new/1.25, center_FH_new(1), center_FH_new(2));
    quality_FH(:,5)=max(max(map_tmp));
    
            map_tmp = spatial_filter_circle(map, radius_FH_new*1.25, center_FH_new(1), center_FH_new(2));
    quality_FH(:,6)=max(max(map_tmp));
    
                map_tmp = spatial_filter_circle(map, radius_FH_new*1.5, center_FH_new(1), center_FH_new(2));
    quality_FH(:,7)=max(max(map_tmp));
    
    
% quality: triradiate cartilage
    sps1_tri_car = sps1;
    index1 = min(center(1,2)+radius_FH_new, size(sps1,1));
    index2 = max(center(1,1)-radius_FH_new,1);
    index3 = min(center(1,1)+radius_FH_new, size(sps1,2));
    sps1_tri_car(index1:end, index2 : index3);
        trirad_hist = reshape(sps1_tri_car,[size(sps1_tri_car,1)*size(sps1_tri_car,2),1]);
        [quality_trirad,n]= hist(trirad_hist); 
    
        
        
        
        
     %   figure, imagesc(map_BW.*im), colormap gray, axis equal   