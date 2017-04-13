
close all
clear all

%% adequate 
cd ..
cd('adequate_numbered')

for k = 1:258

    pngFileName = strcat('a ('...
        , num2str(k), ')','.png');
    img = imread(pngFileName);
    img = double (img (:,:,1) );     
    
   
    row_size = round (size(img,1) / ( size(img,2)/200 ) ) ;
    img = imresize(img, [row_size, 200]); %imagesc(img)
                    %remove the extra zero rows
                    temp =sum(img'); temp2 = max(find(temp>0));
                    img = img(1:temp2,:);
                    im2 = img;
                    %reduce speckle using gaussian, optional.. for new dataset, too grainy
%                     h = fspecial('gaussian', [4 4], 1);
%                     img = imfilter(img, h, 'replicate');
%figure, imagesc(img)

    size_mm = round(size(img,2)/40); %number of pixels per mm

    [r, c]=size(img);
                    ML1=r; ML2=c; ML3=r+c; ML4=r*c;
    % criterions
    angle = 75; % angle for ilium criterion

    % sps computation
    %[sps1, sps2] = sps_DDH_feb15 (img, angle);
    [sps1, sps2] = sps_DDH (img, angle);
                    sps1_hist = reshape(sps1,[r*c,1]);[ML5,n]= hist(sps1_hist); ML5=normalize(ML5);
                    sps2_hist = reshape(sps2,[r*c,1]);[ML6,n]= hist(sps2_hist); ML6=normalize(ML6);
    % bone contour segmentation
    [bone, probability, map] = bone_DDH (img, sps1, sps2);
                    probability_hist = reshape(probability,[r*c,1]);[ML7,n]= hist(probability_hist); ML7=normalize(ML7);
                    img_hist = reshape(img,[r*c,1]);[ML8,n]= hist(img_hist); ML8=normalize(ML8);
    %figure, imagesc(map), axis equal
    bone = bone .* (probability>0.5); % this threshold can be user controlled
                    bone_hist = reshape(bone,[r*c,1]);[ML9,n]= hist(bone_hist); ML9=normalize(ML9);
  
    % ray-cast to avoid the G muscles, may look similar to horizontal ilium    
    bone_rc = ray_cast_DDH (bone);
                    bone_rc_hist = reshape(bone_rc,[r*c,1]);[ML10,n]= hist(bone_rc_hist); ML10=normalize(ML10);
    ilium = bone_rc .* (sps2>0);
                    ilium_hist = reshape(ilium,[r*c,1]);[ML11,n]= hist(ilium_hist); ML11=normalize(ML11);
        
        if nnz(ilium)==0
            ML(k,:)=ML_default;
            continue
        else
        end
% horizontal ilium criteria, length and major angle
    [hor_pared_ilium, hor_length, ver_length, major_angle] = horizontal_ilium (ilium, size_mm) ; 
  %  figure, imagesc(hor_pared_ilium)
    
                    ML12 = nnz(hor_pared_ilium); ML13 = hor_length; ML14 = ver_length; ML15 = major_angle;
 
                    

   
    
%end % end of first test, ilium extraction  
    femoral_head_radius = size_mm * 10; % 10mm at start, default
    
    %for good visualization
            sps1 = sps1 .* (sps1>.1);
            sps1 = normalize(sps1);
    
    [R, radon_strength_1, radon_strength_2, radon_strength_3, ilium_angle, acetab_angle, labrum_angle,...
    alpha, beta, bisector1, bisector2, center_FH, radius_FH, row_acetab, col_acetab]= ilium_acetabulum_labrum2 (hor_pared_ilium,sps1, bone_rc, femoral_head_radius) ;
                    R_hist = reshape(R,[size(R,1)*size(R,2),1]);[ML16,n]= hist(R_hist); ML16=normalize(ML16);
                    ML17=radon_strength_1; ML18=radon_strength_2; ML19=radon_strength_3; ML20=ilium_angle; ML21=acetab_angle;
                    ML22=labrum_angle; ML23=alpha; ML24=beta; ML25=bisector1; ML26=bisector2; ML27=center_FH/size_mm; 
                    ML28=radius_FH/size_mm; ML29=row_acetab/size_mm; ML30=col_acetab/size_mm;


%intermediate visualization
    
%     figure, 
%     tmp = [bisector1;bisector2];
%     show_straight_line(tmp,im2) 
%     viscircles(center_FH, radius_FH,'EdgeColor','b');   


center_alpha_beta = [col_acetab, row_acetab];
    
% Femoral head segmentation
    [center_FH_new, radius_FH_new, quality_FH, quality_trirad] = femoral_head_feb_miccai (sps1,probability,map,radius_FH,center_FH,center_alpha_beta,alpha,beta,size_mm,img);
                    ML31=center_FH_new/size_mm; ML32=radius_FH_new/size_mm; ML33=quality_FH; ML49=quality_trirad;

    if radius_FH_new==0
            ML(k,:)=zeros(1,142);
            continue
        else
        end
    
% angle values refinement
femoral_head_radius = radius_FH_new; 

% visualization
   % figure, imagesc(img), colormap gray;axis equal
   % viscircles(center_FH_new, radius_FH_new,'EdgeColor','b');axis equal

    
    [R, radon_strength_1, radon_strength_2, radon_strength_3, ilium_angle, acetab_angle, labrum_angle,...
    alpha, beta, bisector1, bisector2, center_FH, radius_FH, row_acetab, col_acetab]= ilium_acetabulum_labrum2 (hor_pared_ilium,sps1, bone_rc, femoral_head_radius) ;

                    R_hist = reshape(R,[size(R,1)*size(R,2),1]);[ML34,n]= hist(R_hist); ML34=normalize(ML34);
                    ML35=radon_strength_1; ML36=radon_strength_2; ML37=radon_strength_3; ML38=ilium_angle; ML39=acetab_angle;
                    ML40=labrum_angle; ML41=alpha; ML42=beta; ML43=bisector1; ML44=bisector2; ML45=center_FH/size_mm; 
                    ML46=radius_FH/size_mm; ML47=row_acetab/size_mm; ML48=col_acetab/size_mm;
alpha;
beta;
    
alpha_1(k)=alpha;
beta_1(k)=beta;
    
% method 2 based angle calculation
    
    f_angles = [45 90 135];
    minW=[20 20 20];
    % 2D OPS
    OPS= optimized_phasesym(img,5,3,f_angles,minW,1.25, 1.5,.1,1);%k=.5 here, change to 3 for higher noise threshold
    [R, radon_strength_1, radon_strength_2, radon_strength_3, ilium_angle, acetab_angle, labrum_angle,...
    alpha, beta, bisector1, bisector2, center_FH, radius_FH, row_acetab, col_acetab]= ilium_acetabulum_labrum2 (hor_pared_ilium,OPS, bone_rc, femoral_head_radius) ;

    alpha;
    beta;

alpha_2(k)=alpha;
beta_2(k)=beta;
    
       % figure, imagesc(im2), colormap gray;axis equal
  %  viscircles(center_FH_new, radius_FH_new,'EdgeColor','b','LineWidth',2);axis equal
    
    
ML(k,:)=[ML1,ML2,ML3,ML4,ML5,ML6,ML7,ML8,ML9,ML10,ML11,ML12,ML13,ML14,ML15,ML16,ML17,ML18,ML19,ML20,ML21,ML22,ML23,ML24...
    ,ML25,ML26,ML27,ML28,ML29,ML30,ML31,ML32,ML33,ML34,ML35,ML36,ML37,ML38,ML39,ML40,ML41,ML42,ML43,ML44,ML45,ML46,ML47,ML48];    
 
% % intermediate visualization
%   
%     figure, subplot (2,1,1), imagesc(bone), axis equal, colormap gray
%    
%     
% str=sprintf('Adequate, iteration number = %.0f and feature strength = %.0f', k, ML12);
% title(str)
% 
% 
%     subplot (2,1,2)
%     imagesc(hor_pared_ilium), colormap hot, axis equal
%     pause
%     
%     
adequate_ML12(k) = ML12;
end
    
adequate_ML = ML;
cd ..
cd('classification')
save ('adequate1', 'adequate_ML')
cd ..
%% inadequate - coronal 
clear ML
cd('inadequate_coronal_numbered')

for k = 1:137

    pngFileName = strcat('z ('...
        , num2str(k), ')','.png');
    img = imread(pngFileName);
    img = double (img (:,:,1) );     
    
   
    row_size = round (size(img,1) / ( size(img,2)/200 ) ) ;
    img = imresize(img, [row_size, 200]); %imagesc(img)
                    %remove the extra zero rows
                    temp =sum(img'); temp2 = max(find(temp>0));
                    img = img(1:temp2,:);
                    im2 = img;
                    %reduce speckle using gaussian, optional.. for new dataset, too grainy
%                     h = fspecial('gaussian', [4 4], 1);
%                     img = imfilter(img, h, 'replicate');
%figure, imagesc(img)

    size_mm = round(size(img,2)/40); %number of pixels per mm

    [r, c]=size(img);
                    ML1=r; ML2=c; ML3=r+c; ML4=r*c;
    % criterions
    angle = 75; % angle for ilium criterion

    % sps computation
    %[sps1, sps2] = sps_DDH_feb15 (img, angle);
    [sps1, sps2] = sps_DDH (img, angle);
                    sps1_hist = reshape(sps1,[r*c,1]);[ML5,n]= hist(sps1_hist); ML5=normalize(ML5);
                    sps2_hist = reshape(sps2,[r*c,1]);[ML6,n]= hist(sps2_hist); ML6=normalize(ML6);
    % bone contour segmentation
    [bone, probability, map] = bone_DDH (img, sps1, sps2);
                    probability_hist = reshape(probability,[r*c,1]);[ML7,n]= hist(probability_hist); ML7=normalize(ML7);
                    img_hist = reshape(img,[r*c,1]);[ML8,n]= hist(img_hist); ML8=normalize(ML8);
    %figure, imagesc(map), axis equal
    bone = bone .* (probability>0.5); % this threshold can be user controlled
                    bone_hist = reshape(bone,[r*c,1]);[ML9,n]= hist(bone_hist); ML9=normalize(ML9);
  
    % ray-cast to avoid the G muscles, may look similar to horizontal ilium    
    bone_rc = ray_cast_DDH (bone);
                    bone_rc_hist = reshape(bone_rc,[r*c,1]);[ML10,n]= hist(bone_rc_hist); ML10=normalize(ML10);
    ilium = bone_rc .* (sps2>0);
                    ilium_hist = reshape(ilium,[r*c,1]);[ML11,n]= hist(ilium_hist); ML11=normalize(ML11);
        
        if nnz(ilium)==0
            ML(k,:)=ML_default;
            continue
        else
        end
% horizontal ilium criteria, length and major angle
    [hor_pared_ilium, hor_length, ver_length, major_angle] = horizontal_ilium (ilium, size_mm) ; 
  %  figure, imagesc(hor_pared_ilium)
    
                    ML12 = nnz(hor_pared_ilium); ML13 = hor_length; ML14 = ver_length; ML15 = major_angle;
 
                    

   
    
%end % end of first test, ilium extraction  
    femoral_head_radius = size_mm * 10; % 10mm at start, default
    
    %for good visualization
            sps1 = sps1 .* (sps1>.1);
            sps1 = normalize(sps1);
    
    [R, radon_strength_1, radon_strength_2, radon_strength_3, ilium_angle, acetab_angle, labrum_angle,...
    alpha, beta, bisector1, bisector2, center_FH, radius_FH, row_acetab, col_acetab]= ilium_acetabulum_labrum2 (hor_pared_ilium,sps1, bone_rc, femoral_head_radius) ;
                    R_hist = reshape(R,[size(R,1)*size(R,2),1]);[ML16,n]= hist(R_hist); ML16=normalize(ML16);
                    ML17=radon_strength_1; ML18=radon_strength_2; ML19=radon_strength_3; ML20=ilium_angle; ML21=acetab_angle;
                    ML22=labrum_angle; ML23=alpha; ML24=beta; ML25=bisector1; ML26=bisector2; ML27=center_FH/size_mm; 
                    ML28=radius_FH/size_mm; ML29=row_acetab/size_mm; ML30=col_acetab/size_mm;



center_alpha_beta = [col_acetab, row_acetab];
    
% Femoral head segmentation
    [center_FH_new, radius_FH_new, quality_FH, quality_trirad] = femoral_head_feb_miccai (sps1,probability,map,radius_FH,center_FH,center_alpha_beta,alpha,beta,size_mm,img);
                    ML31=center_FH_new/size_mm; ML32=radius_FH_new/size_mm; ML33=quality_FH; ML49=quality_trirad;

    if radius_FH_new==0
            ML(k,:)=zeros(1,142);
            continue
        else
        end
    
% angle values refinement
femoral_head_radius = radius_FH_new; 
    
    [R, radon_strength_1, radon_strength_2, radon_strength_3, ilium_angle, acetab_angle, labrum_angle,...
    alpha, beta, bisector1, bisector2, center_FH, radius_FH, row_acetab, col_acetab]= ilium_acetabulum_labrum2 (hor_pared_ilium,sps1, bone_rc, femoral_head_radius) ;

                    R_hist = reshape(R,[size(R,1)*size(R,2),1]);[ML34,n]= hist(R_hist); ML34=normalize(ML34);
                    ML35=radon_strength_1; ML36=radon_strength_2; ML37=radon_strength_3; ML38=ilium_angle; ML39=acetab_angle;
                    ML40=labrum_angle; ML41=alpha; ML42=beta; ML43=bisector1; ML44=bisector2; ML45=center_FH/size_mm; 
                    ML46=radius_FH/size_mm; ML47=row_acetab/size_mm; ML48=col_acetab/size_mm;
alpha;
beta;
    
alpha_1(k)=alpha;
beta_1(k)=beta;
    
% method 2 based angle calculation
    
    f_angles = [45 90 135];
    minW=[20 20 20];
    % 2D OPS
    OPS= optimized_phasesym(img,5,3,f_angles,minW,1.25, 1.5,.1,1);%k=.5 here, change to 3 for higher noise threshold
    [R, radon_strength_1, radon_strength_2, radon_strength_3, ilium_angle, acetab_angle, labrum_angle,...
    alpha, beta, bisector1, bisector2, center_FH, radius_FH, row_acetab, col_acetab]= ilium_acetabulum_labrum2 (hor_pared_ilium,OPS, bone_rc, femoral_head_radius) ;

    alpha;
    beta;

alpha_2(k)=alpha;
beta_2(k)=beta;
    
ML(k,:)=[ML1,ML2,ML3,ML4,ML5,ML6,ML7,ML8,ML9,ML10,ML11,ML12,ML13,ML14,ML15,ML16,ML17,ML18,ML19,ML20,ML21,ML22,ML23,ML24...
    ,ML25,ML26,ML27,ML28,ML29,ML30,ML31,ML32,ML33,ML34,ML35,ML36,ML37,ML38,ML39,ML40,ML41,ML42,ML43,ML44,ML45,ML46,ML47,ML48];    
 
% % intermediate visualization
%   
%     figure, subplot (2,1,1), imagesc(bone), axis equal, colormap gray
%    
%     
% str=sprintf('Inadequate_C, iteration number = %.0f and feature strength = %.0f', k, ML12);
% title(str)
% 
% 
%     subplot (2,1,2)
%     imagesc(hor_pared_ilium), colormap hot, axis equal
    
inadequate_coronal_ML12(k) = ML12;

end

inadequate_coronal_ML = ML;
cd ..
cd('classification')
save ('inadequate_coronal1', 'inadequate_coronal_ML')
cd ..

%% inadequate - transverse 
clear ML
cd('inadequate_transverse_numbered')

for k = 1:298

    pngFileName = strcat('z ('...
        , num2str(k), ')','.png');
    img = imread(pngFileName);
    img = double (img (:,:,1) );     
    
   
    row_size = round (size(img,1) / ( size(img,2)/200 ) ) ;
    img = imresize(img, [row_size, 200]); %imagesc(img)
                    %remove the extra zero rows
                    temp =sum(img'); temp2 = max(find(temp>0));
                    img = img(1:temp2,:);
                    im2 = img;
                    %reduce speckle using gaussian, optional.. for new dataset, too grainy
%                     h = fspecial('gaussian', [4 4], 1);
%                     img = imfilter(img, h, 'replicate');
%figure, imagesc(img)

    size_mm = round(size(img,2)/40); %number of pixels per mm

    [r, c]=size(img);
                    ML1=r; ML2=c; ML3=r+c; ML4=r*c;
    % criterions
    angle = 75; % angle for ilium criterion

    % sps computation
    %[sps1, sps2] = sps_DDH_feb15 (img, angle);
    [sps1, sps2] = sps_DDH (img, angle);
                    sps1_hist = reshape(sps1,[r*c,1]);[ML5,n]= hist(sps1_hist); ML5=normalize(ML5);
                    sps2_hist = reshape(sps2,[r*c,1]);[ML6,n]= hist(sps2_hist); ML6=normalize(ML6);
    % bone contour segmentation
    [bone, probability, map] = bone_DDH (img, sps1, sps2);
                    probability_hist = reshape(probability,[r*c,1]);[ML7,n]= hist(probability_hist); ML7=normalize(ML7);
                    img_hist = reshape(img,[r*c,1]);[ML8,n]= hist(img_hist); ML8=normalize(ML8);
    %figure, imagesc(map), axis equal
    bone = bone .* (probability>0.5); % this threshold can be user controlled
                    bone_hist = reshape(bone,[r*c,1]);[ML9,n]= hist(bone_hist); ML9=normalize(ML9);
  
    % ray-cast to avoid the G muscles, may look similar to horizontal ilium    
    bone_rc = ray_cast_DDH (bone);
                    bone_rc_hist = reshape(bone_rc,[r*c,1]);[ML10,n]= hist(bone_rc_hist); ML10=normalize(ML10);
    ilium = bone_rc .* (sps2>0);
                    ilium_hist = reshape(ilium,[r*c,1]);[ML11,n]= hist(ilium_hist); ML11=normalize(ML11);
        
        if nnz(ilium)==0
            ML(k,:)=ML_default;
            continue
        else
        end
% horizontal ilium criteria, length and major angle
    [hor_pared_ilium, hor_length, ver_length, major_angle] = horizontal_ilium (ilium, size_mm) ; 
  %  figure, imagesc(hor_pared_ilium)
    
                    ML12 = nnz(hor_pared_ilium); ML13 = hor_length; ML14 = ver_length; ML15 = major_angle;
 
                    

   
    
%end % end of first test, ilium extraction  
    femoral_head_radius = size_mm * 10; % 10mm at start, default
    
    %for good visualization
            sps1 = sps1 .* (sps1>.1);
            sps1 = normalize(sps1);
    
    [R, radon_strength_1, radon_strength_2, radon_strength_3, ilium_angle, acetab_angle, labrum_angle,...
    alpha, beta, bisector1, bisector2, center_FH, radius_FH, row_acetab, col_acetab]= ilium_acetabulum_labrum2 (hor_pared_ilium,sps1, bone_rc, femoral_head_radius) ;
                    R_hist = reshape(R,[size(R,1)*size(R,2),1]);[ML16,n]= hist(R_hist); ML16=normalize(ML16);
                    ML17=radon_strength_1; ML18=radon_strength_2; ML19=radon_strength_3; ML20=ilium_angle; ML21=acetab_angle;
                    ML22=labrum_angle; ML23=alpha; ML24=beta; ML25=bisector1; ML26=bisector2; ML27=center_FH/size_mm; 
                    ML28=radius_FH/size_mm; ML29=row_acetab/size_mm; ML30=col_acetab/size_mm;


%intermediate visualization
    
%     figure, 
%     tmp = [bisector1;bisector2];
%     show_straight_line(tmp,im2) 
%     viscircles(center_FH, radius_FH,'EdgeColor','b');   


center_alpha_beta = [col_acetab, row_acetab];
    
% Femoral head segmentation
    [center_FH_new, radius_FH_new, quality_FH, quality_trirad] = femoral_head_feb_miccai (sps1,probability,map,radius_FH,center_FH,center_alpha_beta,alpha,beta,size_mm,img);
                    ML31=center_FH_new/size_mm; ML32=radius_FH_new/size_mm; ML33=quality_FH; ML49=quality_trirad;

    if radius_FH_new==0
            ML(k,:)=zeros(1,142);
            continue
        else
        end
    
% angle values refinement
femoral_head_radius = radius_FH_new; 

% visualization
   % figure, imagesc(img), colormap gray;axis equal
   % viscircles(center_FH_new, radius_FH_new,'EdgeColor','b');axis equal

    
    [R, radon_strength_1, radon_strength_2, radon_strength_3, ilium_angle, acetab_angle, labrum_angle,...
    alpha, beta, bisector1, bisector2, center_FH, radius_FH, row_acetab, col_acetab]= ilium_acetabulum_labrum2 (hor_pared_ilium,sps1, bone_rc, femoral_head_radius) ;

                    R_hist = reshape(R,[size(R,1)*size(R,2),1]);[ML34,n]= hist(R_hist); ML34=normalize(ML34);
                    ML35=radon_strength_1; ML36=radon_strength_2; ML37=radon_strength_3; ML38=ilium_angle; ML39=acetab_angle;
                    ML40=labrum_angle; ML41=alpha; ML42=beta; ML43=bisector1; ML44=bisector2; ML45=center_FH/size_mm; 
                    ML46=radius_FH/size_mm; ML47=row_acetab/size_mm; ML48=col_acetab/size_mm;
alpha;
beta;
    
alpha_1(k)=alpha;
beta_1(k)=beta;
    
% method 2 based angle calculation
    
    f_angles = [45 90 135];
    minW=[20 20 20];
    % 2D OPS
    OPS= optimized_phasesym(img,5,3,f_angles,minW,1.25, 1.5,.1,1);%k=.5 here, change to 3 for higher noise threshold
    [R, radon_strength_1, radon_strength_2, radon_strength_3, ilium_angle, acetab_angle, labrum_angle,...
    alpha, beta, bisector1, bisector2, center_FH, radius_FH, row_acetab, col_acetab]= ilium_acetabulum_labrum2 (hor_pared_ilium,OPS, bone_rc, femoral_head_radius) ;

    alpha;
    beta;

alpha_2(k)=alpha;
beta_2(k)=beta;
    
       % figure, imagesc(im2), colormap gray;axis equal
  %  viscircles(center_FH_new, radius_FH_new,'EdgeColor','b','LineWidth',2);axis equal
    
    
ML(k,:)=[ML1,ML2,ML3,ML4,ML5,ML6,ML7,ML8,ML9,ML10,ML11,ML12,ML13,ML14,ML15,ML16,ML17,ML18,ML19,ML20,ML21,ML22,ML23,ML24...
    ,ML25,ML26,ML27,ML28,ML29,ML30,ML31,ML32,ML33,ML34,ML35,ML36,ML37,ML38,ML39,ML40,ML41,ML42,ML43,ML44,ML45,ML46,ML47,ML48];    
 
% % intermediate visualization
%   
%     figure, subplot (2,1,1), imagesc(bone), axis equal, colormap gray
%    
%     
% str=sprintf('Inadequate_T, iteration number = %.0f and feature strength = %.0f', k, ML12);
% title(str)
% 
% 
%     subplot (2,1,2)
%     imagesc(hor_pared_ilium), colormap hot, axis equal
    
inadequate_transverse_ML12(k) = ML12;
end

inadequate_transverse_ML = ML;
cd ..
cd('classification')
save ('inadequate_transverse1', 'inadequate_transverse_ML')
cd ..

%% adequate - holdout 
clear ML
cd('adequate_numbered')

for k = 1:68

    pngFileName = strcat('holdout ('...
        , num2str(k), ')','.png');
    img = imread(pngFileName);
    img = double (img (:,:,1) );     
    
   
    row_size = round (size(img,1) / ( size(img,2)/200 ) ) ;
    img = imresize(img, [row_size, 200]); %imagesc(img)
                    %remove the extra zero rows
                    temp =sum(img'); temp2 = max(find(temp>0));
                    img = img(1:temp2,:);
                    im2 = img;
                    %reduce speckle using gaussian, optional.. for new dataset, too grainy
%                     h = fspecial('gaussian', [4 4], 1);
%                     img = imfilter(img, h, 'replicate');
%figure, imagesc(img)

    size_mm = round(size(img,2)/40); %number of pixels per mm

    [r, c]=size(img);
                    ML1=r; ML2=c; ML3=r+c; ML4=r*c;
    % criterions
    angle = 75; % angle for ilium criterion

    % sps computation
    %[sps1, sps2] = sps_DDH_feb15 (img, angle);
    [sps1, sps2] = sps_DDH (img, angle);
                    sps1_hist = reshape(sps1,[r*c,1]);[ML5,n]= hist(sps1_hist); ML5=normalize(ML5);
                    sps2_hist = reshape(sps2,[r*c,1]);[ML6,n]= hist(sps2_hist); ML6=normalize(ML6);
    % bone contour segmentation
    [bone, probability, map] = bone_DDH (img, sps1, sps2);
                    probability_hist = reshape(probability,[r*c,1]);[ML7,n]= hist(probability_hist); ML7=normalize(ML7);
                    img_hist = reshape(img,[r*c,1]);[ML8,n]= hist(img_hist); ML8=normalize(ML8);
    %figure, imagesc(map), axis equal
    bone = bone .* (probability>0.5); % this threshold can be user controlled
                    bone_hist = reshape(bone,[r*c,1]);[ML9,n]= hist(bone_hist); ML9=normalize(ML9);
  
    % ray-cast to avoid the G muscles, may look similar to horizontal ilium    
    bone_rc = ray_cast_DDH (bone);
                    bone_rc_hist = reshape(bone_rc,[r*c,1]);[ML10,n]= hist(bone_rc_hist); ML10=normalize(ML10);
    ilium = bone_rc .* (sps2>0);
                    ilium_hist = reshape(ilium,[r*c,1]);[ML11,n]= hist(ilium_hist); ML11=normalize(ML11);
        
        if nnz(ilium)==0
            ML(k,:)=ML_default;
            continue
        else
        end
% horizontal ilium criteria, length and major angle
    [hor_pared_ilium, hor_length, ver_length, major_angle] = horizontal_ilium (ilium, size_mm) ; 
  %  figure, imagesc(hor_pared_ilium)
    
                    ML12 = nnz(hor_pared_ilium); ML13 = hor_length; ML14 = ver_length; ML15 = major_angle;
 
                    

   
    
%end % end of first test, ilium extraction  
    femoral_head_radius = size_mm * 10; % 10mm at start, default
    
    %for good visualization
            sps1 = sps1 .* (sps1>.1);
            sps1 = normalize(sps1);
    
    [R, radon_strength_1, radon_strength_2, radon_strength_3, ilium_angle, acetab_angle, labrum_angle,...
    alpha, beta, bisector1, bisector2, center_FH, radius_FH, row_acetab, col_acetab]= ilium_acetabulum_labrum2 (hor_pared_ilium,sps1, bone_rc, femoral_head_radius) ;
                    R_hist = reshape(R,[size(R,1)*size(R,2),1]);[ML16,n]= hist(R_hist); ML16=normalize(ML16);
                    ML17=radon_strength_1; ML18=radon_strength_2; ML19=radon_strength_3; ML20=ilium_angle; ML21=acetab_angle;
                    ML22=labrum_angle; ML23=alpha; ML24=beta; ML25=bisector1; ML26=bisector2; ML27=center_FH/size_mm; 
                    ML28=radius_FH/size_mm; ML29=row_acetab/size_mm; ML30=col_acetab/size_mm;


%intermediate visualization
    
%     figure, 
%     tmp = [bisector1;bisector2];
%     show_straight_line(tmp,im2) 
%     viscircles(center_FH, radius_FH,'EdgeColor','b');   


center_alpha_beta = [col_acetab, row_acetab];
    
% Femoral head segmentation
    [center_FH_new, radius_FH_new, quality_FH, quality_trirad] = femoral_head_feb_miccai (sps1,probability,map,radius_FH,center_FH,center_alpha_beta,alpha,beta,size_mm,img);
                    ML31=center_FH_new/size_mm; ML32=radius_FH_new/size_mm; ML33=quality_FH; ML49=quality_trirad;

    if radius_FH_new==0
            ML(k,:)=zeros(1,142);
            continue
        else
        end
    
% angle values refinement
femoral_head_radius = radius_FH_new; 

% visualization
   % figure, imagesc(img), colormap gray;axis equal
   % viscircles(center_FH_new, radius_FH_new,'EdgeColor','b');axis equal

    
    [R, radon_strength_1, radon_strength_2, radon_strength_3, ilium_angle, acetab_angle, labrum_angle,...
    alpha, beta, bisector1, bisector2, center_FH, radius_FH, row_acetab, col_acetab]= ilium_acetabulum_labrum2 (hor_pared_ilium,sps1, bone_rc, femoral_head_radius) ;

                    R_hist = reshape(R,[size(R,1)*size(R,2),1]);[ML34,n]= hist(R_hist); ML34=normalize(ML34);
                    ML35=radon_strength_1; ML36=radon_strength_2; ML37=radon_strength_3; ML38=ilium_angle; ML39=acetab_angle;
                    ML40=labrum_angle; ML41=alpha; ML42=beta; ML43=bisector1; ML44=bisector2; ML45=center_FH/size_mm; 
                    ML46=radius_FH/size_mm; ML47=row_acetab/size_mm; ML48=col_acetab/size_mm;
alpha;
beta;
    
alpha_1(k)=alpha;
beta_1(k)=beta;
    
% method 2 based angle calculation
    
    f_angles = [45 90 135];
    minW=[20 20 20];
    % 2D OPS
    OPS= optimized_phasesym(img,5,3,f_angles,minW,1.25, 1.5,.1,1);%k=.5 here, change to 3 for higher noise threshold
    [R, radon_strength_1, radon_strength_2, radon_strength_3, ilium_angle, acetab_angle, labrum_angle,...
    alpha, beta, bisector1, bisector2, center_FH, radius_FH, row_acetab, col_acetab]= ilium_acetabulum_labrum2 (hor_pared_ilium,OPS, bone_rc, femoral_head_radius) ;

    alpha;
    beta;

alpha_2(k)=alpha;
beta_2(k)=beta;
    
       % figure, imagesc(im2), colormap gray;axis equal
  %  viscircles(center_FH_new, radius_FH_new,'EdgeColor','b','LineWidth',2);axis equal
    
    
ML(k,:)=[ML1,ML2,ML3,ML4,ML5,ML6,ML7,ML8,ML9,ML10,ML11,ML12,ML13,ML14,ML15,ML16,ML17,ML18,ML19,ML20,ML21,ML22,ML23,ML24...
    ,ML25,ML26,ML27,ML28,ML29,ML30,ML31,ML32,ML33,ML34,ML35,ML36,ML37,ML38,ML39,ML40,ML41,ML42,ML43,ML44,ML45,ML46,ML47,ML48];    
 
% % intermediate visualization
%   
%     figure, subplot (2,1,1), imagesc(bone), axis equal, colormap gray
%    
%     
% str=sprintf('Adequate, iteration number = %.0f and feature strength = %.0f', k, ML12);
% title(str)
% 
% 
%     subplot (2,1,2)
%     imagesc(hor_pared_ilium), colormap hot, axis equal
%     pause
%     
%     
adequate_ML12(k) = ML12;
end
    
adequate_holdout_ML = ML;
cd ..
cd('classification')
save ('adequate1_holdout', 'adequate_holdout_ML')
cd ..
%% inadequate - coronal - holdout
clear ML
cd('inadequate_coronal_numbered')

for k = 1:38

    pngFileName = strcat('holdout ('...
        , num2str(k), ')','.png');
    img = imread(pngFileName);
    img = double (img (:,:,1) );     
    
   
    row_size = round (size(img,1) / ( size(img,2)/200 ) ) ;
    img = imresize(img, [row_size, 200]); %imagesc(img)
                    %remove the extra zero rows
                    temp =sum(img'); temp2 = max(find(temp>0));
                    img = img(1:temp2,:);
                    im2 = img;
                    %reduce speckle using gaussian, optional.. for new dataset, too grainy
%                     h = fspecial('gaussian', [4 4], 1);
%                     img = imfilter(img, h, 'replicate');
%figure, imagesc(img)

    size_mm = round(size(img,2)/40); %number of pixels per mm

    [r, c]=size(img);
                    ML1=r; ML2=c; ML3=r+c; ML4=r*c;
    % criterions
    angle = 75; % angle for ilium criterion

    % sps computation
    %[sps1, sps2] = sps_DDH_feb15 (img, angle);
    [sps1, sps2] = sps_DDH (img, angle);
                    sps1_hist = reshape(sps1,[r*c,1]);[ML5,n]= hist(sps1_hist); ML5=normalize(ML5);
                    sps2_hist = reshape(sps2,[r*c,1]);[ML6,n]= hist(sps2_hist); ML6=normalize(ML6);
    % bone contour segmentation
    [bone, probability, map] = bone_DDH (img, sps1, sps2);
                    probability_hist = reshape(probability,[r*c,1]);[ML7,n]= hist(probability_hist); ML7=normalize(ML7);
                    img_hist = reshape(img,[r*c,1]);[ML8,n]= hist(img_hist); ML8=normalize(ML8);
    %figure, imagesc(map), axis equal
    bone = bone .* (probability>0.5); % this threshold can be user controlled
                    bone_hist = reshape(bone,[r*c,1]);[ML9,n]= hist(bone_hist); ML9=normalize(ML9);
  
    % ray-cast to avoid the G muscles, may look similar to horizontal ilium    
    bone_rc = ray_cast_DDH (bone);
                    bone_rc_hist = reshape(bone_rc,[r*c,1]);[ML10,n]= hist(bone_rc_hist); ML10=normalize(ML10);
    ilium = bone_rc .* (sps2>0);
                    ilium_hist = reshape(ilium,[r*c,1]);[ML11,n]= hist(ilium_hist); ML11=normalize(ML11);
        
        if nnz(ilium)==0
            ML(k,:)=ML_default;
            continue
        else
        end
% horizontal ilium criteria, length and major angle
    [hor_pared_ilium, hor_length, ver_length, major_angle] = horizontal_ilium (ilium, size_mm) ; 
  %  figure, imagesc(hor_pared_ilium)
    
                    ML12 = nnz(hor_pared_ilium); ML13 = hor_length; ML14 = ver_length; ML15 = major_angle;
 
                    

   
    
%end % end of first test, ilium extraction  
    femoral_head_radius = size_mm * 10; % 10mm at start, default
    
    %for good visualization
            sps1 = sps1 .* (sps1>.1);
            sps1 = normalize(sps1);
    
    [R, radon_strength_1, radon_strength_2, radon_strength_3, ilium_angle, acetab_angle, labrum_angle,...
    alpha, beta, bisector1, bisector2, center_FH, radius_FH, row_acetab, col_acetab]= ilium_acetabulum_labrum2 (hor_pared_ilium,sps1, bone_rc, femoral_head_radius) ;
                    R_hist = reshape(R,[size(R,1)*size(R,2),1]);[ML16,n]= hist(R_hist); ML16=normalize(ML16);
                    ML17=radon_strength_1; ML18=radon_strength_2; ML19=radon_strength_3; ML20=ilium_angle; ML21=acetab_angle;
                    ML22=labrum_angle; ML23=alpha; ML24=beta; ML25=bisector1; ML26=bisector2; ML27=center_FH/size_mm; 
                    ML28=radius_FH/size_mm; ML29=row_acetab/size_mm; ML30=col_acetab/size_mm;



center_alpha_beta = [col_acetab, row_acetab];
    
% Femoral head segmentation
    [center_FH_new, radius_FH_new, quality_FH, quality_trirad] = femoral_head_feb_miccai (sps1,probability,map,radius_FH,center_FH,center_alpha_beta,alpha,beta,size_mm,img);
                    ML31=center_FH_new/size_mm; ML32=radius_FH_new/size_mm; ML33=quality_FH; ML49=quality_trirad;

    if radius_FH_new==0
            ML(k,:)=zeros(1,142);
            continue
        else
        end
    
% angle values refinement
femoral_head_radius = radius_FH_new; 
    
    [R, radon_strength_1, radon_strength_2, radon_strength_3, ilium_angle, acetab_angle, labrum_angle,...
    alpha, beta, bisector1, bisector2, center_FH, radius_FH, row_acetab, col_acetab]= ilium_acetabulum_labrum2 (hor_pared_ilium,sps1, bone_rc, femoral_head_radius) ;

                    R_hist = reshape(R,[size(R,1)*size(R,2),1]);[ML34,n]= hist(R_hist); ML34=normalize(ML34);
                    ML35=radon_strength_1; ML36=radon_strength_2; ML37=radon_strength_3; ML38=ilium_angle; ML39=acetab_angle;
                    ML40=labrum_angle; ML41=alpha; ML42=beta; ML43=bisector1; ML44=bisector2; ML45=center_FH/size_mm; 
                    ML46=radius_FH/size_mm; ML47=row_acetab/size_mm; ML48=col_acetab/size_mm;
alpha;
beta;
    
alpha_1(k)=alpha;
beta_1(k)=beta;
    
% method 2 based angle calculation
    
    f_angles = [45 90 135];
    minW=[20 20 20];
    % 2D OPS
    OPS= optimized_phasesym(img,5,3,f_angles,minW,1.25, 1.5,.1,1);%k=.5 here, change to 3 for higher noise threshold
    [R, radon_strength_1, radon_strength_2, radon_strength_3, ilium_angle, acetab_angle, labrum_angle,...
    alpha, beta, bisector1, bisector2, center_FH, radius_FH, row_acetab, col_acetab]= ilium_acetabulum_labrum2 (hor_pared_ilium,OPS, bone_rc, femoral_head_radius) ;

    alpha;
    beta;

alpha_2(k)=alpha;
beta_2(k)=beta;
    
ML(k,:)=[ML1,ML2,ML3,ML4,ML5,ML6,ML7,ML8,ML9,ML10,ML11,ML12,ML13,ML14,ML15,ML16,ML17,ML18,ML19,ML20,ML21,ML22,ML23,ML24...
    ,ML25,ML26,ML27,ML28,ML29,ML30,ML31,ML32,ML33,ML34,ML35,ML36,ML37,ML38,ML39,ML40,ML41,ML42,ML43,ML44,ML45,ML46,ML47,ML48];    
 
% % intermediate visualization
%   
%     figure, subplot (2,1,1), imagesc(bone), axis equal, colormap gray
%    
%     
% str=sprintf('Inadequate_C, iteration number = %.0f and feature strength = %.0f', k, ML12);
% title(str)
% 
% 
%     subplot (2,1,2)
%     imagesc(hor_pared_ilium), colormap hot, axis equal
    
inadequate_coronal_ML12(k) = ML12;

end

inadequate_coronal_holdout_ML = ML;
cd ..
cd('classification')
save ('inadequate_coronal1_holdout', 'inadequate_coronal_holdout_ML')
cd ..

%% inadequate - transverse - holdout
clear ML
cd('inadequate_transverse_numbered')

for k = 1:68

    pngFileName = strcat('holdout ('...
        , num2str(k), ')','.png');
    img = imread(pngFileName);
    img = double (img (:,:,1) );     
    
   
    row_size = round (size(img,1) / ( size(img,2)/200 ) ) ;
    img = imresize(img, [row_size, 200]); %imagesc(img)
                    %remove the extra zero rows
                    temp =sum(img'); temp2 = max(find(temp>0));
                    img = img(1:temp2,:);
                    im2 = img;
                    %reduce speckle using gaussian, optional.. for new dataset, too grainy
%                     h = fspecial('gaussian', [4 4], 1);
%                     img = imfilter(img, h, 'replicate');
%figure, imagesc(img)

    size_mm = round(size(img,2)/40); %number of pixels per mm

    [r, c]=size(img);
                    ML1=r; ML2=c; ML3=r+c; ML4=r*c;
    % criterions
    angle = 75; % angle for ilium criterion

    % sps computation
    %[sps1, sps2] = sps_DDH_feb15 (img, angle);
    [sps1, sps2] = sps_DDH (img, angle);
                    sps1_hist = reshape(sps1,[r*c,1]);[ML5,n]= hist(sps1_hist); ML5=normalize(ML5);
                    sps2_hist = reshape(sps2,[r*c,1]);[ML6,n]= hist(sps2_hist); ML6=normalize(ML6);
    % bone contour segmentation
    [bone, probability, map] = bone_DDH (img, sps1, sps2);
                    probability_hist = reshape(probability,[r*c,1]);[ML7,n]= hist(probability_hist); ML7=normalize(ML7);
                    img_hist = reshape(img,[r*c,1]);[ML8,n]= hist(img_hist); ML8=normalize(ML8);
    %figure, imagesc(map), axis equal
    bone = bone .* (probability>0.5); % this threshold can be user controlled
                    bone_hist = reshape(bone,[r*c,1]);[ML9,n]= hist(bone_hist); ML9=normalize(ML9);
  
    % ray-cast to avoid the G muscles, may look similar to horizontal ilium    
    bone_rc = ray_cast_DDH (bone);
                    bone_rc_hist = reshape(bone_rc,[r*c,1]);[ML10,n]= hist(bone_rc_hist); ML10=normalize(ML10);
    ilium = bone_rc .* (sps2>0);
                    ilium_hist = reshape(ilium,[r*c,1]);[ML11,n]= hist(ilium_hist); ML11=normalize(ML11);
        
        if nnz(ilium)==0
            ML(k,:)=ML_default;
            continue
        else
        end
% horizontal ilium criteria, length and major angle
    [hor_pared_ilium, hor_length, ver_length, major_angle] = horizontal_ilium (ilium, size_mm) ; 
  %  figure, imagesc(hor_pared_ilium)
    
                    ML12 = nnz(hor_pared_ilium); ML13 = hor_length; ML14 = ver_length; ML15 = major_angle;
 
                    

   
    
%end % end of first test, ilium extraction  
    femoral_head_radius = size_mm * 10; % 10mm at start, default
    
    %for good visualization
            sps1 = sps1 .* (sps1>.1);
            sps1 = normalize(sps1);
    
    [R, radon_strength_1, radon_strength_2, radon_strength_3, ilium_angle, acetab_angle, labrum_angle,...
    alpha, beta, bisector1, bisector2, center_FH, radius_FH, row_acetab, col_acetab]= ilium_acetabulum_labrum2 (hor_pared_ilium,sps1, bone_rc, femoral_head_radius) ;
                    R_hist = reshape(R,[size(R,1)*size(R,2),1]);[ML16,n]= hist(R_hist); ML16=normalize(ML16);
                    ML17=radon_strength_1; ML18=radon_strength_2; ML19=radon_strength_3; ML20=ilium_angle; ML21=acetab_angle;
                    ML22=labrum_angle; ML23=alpha; ML24=beta; ML25=bisector1; ML26=bisector2; ML27=center_FH/size_mm; 
                    ML28=radius_FH/size_mm; ML29=row_acetab/size_mm; ML30=col_acetab/size_mm;


%intermediate visualization
    
%     figure, 
%     tmp = [bisector1;bisector2];
%     show_straight_line(tmp,im2) 
%     viscircles(center_FH, radius_FH,'EdgeColor','b');   


center_alpha_beta = [col_acetab, row_acetab];
    
% Femoral head segmentation
    [center_FH_new, radius_FH_new, quality_FH, quality_trirad] = femoral_head_feb_miccai (sps1,probability,map,radius_FH,center_FH,center_alpha_beta,alpha,beta,size_mm,img);
                    ML31=center_FH_new/size_mm; ML32=radius_FH_new/size_mm; ML33=quality_FH; ML49=quality_trirad;

    if radius_FH_new==0
            ML(k,:)=zeros(1,142);
            continue
        else
        end
    
% angle values refinement
femoral_head_radius = radius_FH_new; 

% visualization
   % figure, imagesc(img), colormap gray;axis equal
   % viscircles(center_FH_new, radius_FH_new,'EdgeColor','b');axis equal

    
    [R, radon_strength_1, radon_strength_2, radon_strength_3, ilium_angle, acetab_angle, labrum_angle,...
    alpha, beta, bisector1, bisector2, center_FH, radius_FH, row_acetab, col_acetab]= ilium_acetabulum_labrum2 (hor_pared_ilium,sps1, bone_rc, femoral_head_radius) ;

                    R_hist = reshape(R,[size(R,1)*size(R,2),1]);[ML34,n]= hist(R_hist); ML34=normalize(ML34);
                    ML35=radon_strength_1; ML36=radon_strength_2; ML37=radon_strength_3; ML38=ilium_angle; ML39=acetab_angle;
                    ML40=labrum_angle; ML41=alpha; ML42=beta; ML43=bisector1; ML44=bisector2; ML45=center_FH/size_mm; 
                    ML46=radius_FH/size_mm; ML47=row_acetab/size_mm; ML48=col_acetab/size_mm;
alpha;
beta;
    
alpha_1(k)=alpha;
beta_1(k)=beta;
    
% method 2 based angle calculation
    
    f_angles = [45 90 135];
    minW=[20 20 20];
    % 2D OPS
    OPS= optimized_phasesym(img,5,3,f_angles,minW,1.25, 1.5,.1,1);%k=.5 here, change to 3 for higher noise threshold
    [R, radon_strength_1, radon_strength_2, radon_strength_3, ilium_angle, acetab_angle, labrum_angle,...
    alpha, beta, bisector1, bisector2, center_FH, radius_FH, row_acetab, col_acetab]= ilium_acetabulum_labrum2 (hor_pared_ilium,OPS, bone_rc, femoral_head_radius) ;

    alpha;
    beta;

alpha_2(k)=alpha;
beta_2(k)=beta;
    
       % figure, imagesc(im2), colormap gray;axis equal
  %  viscircles(center_FH_new, radius_FH_new,'EdgeColor','b','LineWidth',2);axis equal
    
    
ML(k,:)=[ML1,ML2,ML3,ML4,ML5,ML6,ML7,ML8,ML9,ML10,ML11,ML12,ML13,ML14,ML15,ML16,ML17,ML18,ML19,ML20,ML21,ML22,ML23,ML24...
    ,ML25,ML26,ML27,ML28,ML29,ML30,ML31,ML32,ML33,ML34,ML35,ML36,ML37,ML38,ML39,ML40,ML41,ML42,ML43,ML44,ML45,ML46,ML47,ML48];    
 
% % intermediate visualization
%   
%     figure, subplot (2,1,1), imagesc(bone), axis equal, colormap gray
%    
%     
% str=sprintf('Inadequate_T, iteration number = %.0f and feature strength = %.0f', k, ML12);
% title(str)
% 
% 
%     subplot (2,1,2)
%     imagesc(hor_pared_ilium), colormap hot, axis equal
    
inadequate_transverse_ML12(k) = ML12;
end

inadequate_transverse_holdout_ML = ML;
cd ..
cd('classification')
save ('inadequate_transverse1_holdout', 'inadequate_transverse_holdout_ML')
cd ..

