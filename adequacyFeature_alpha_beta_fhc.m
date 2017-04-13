function [img, feature, alpha_1, beta_1, fhc, I,L,A,edge_r, edge_c,center_final, radius_final]=adequacyFeature_alpha_beta_fhc(img)


% initialize variables
feature = zeros(1,142);
alpha_1=0; beta_1=0; fhc=0; I=0; L=0; A=0; edge_r=0; edge_c=0; center_final=0; radius_final=0;

 h = waitbar(0,'Please wait...');

 
%Colleen- trying to add img2 to preserve the image size for my GUI. Should
%remove img2 from the function outputs above if this doesn't work. 
img2 = img;


% resize image, so 5 pixels approxximately equals 1mm
row_size = round (size(img,1) / ( size(img,2)/200 ) ) ;
img = imresize(img, [row_size, 200]); 
                    %remove the extra zero rows to increase speed
                    temp =sum(img'); temp2 = max(find(temp>0));
                    img = img(1:temp2,:);   
size_mm = round(size(img,2)/40); %number of pixels per mm

%Colleen- trying to add img2 to preserve the image size for my GUI


[r, c]=size(img);
imagesc(img)
% ML are the features
ML1=r; ML2=c; ML3=r+c; ML4=r*c;

% ilium criteria - we expect ilium to lie 15 degrees (90-15=75 degree) from
% the horizontal axis
angle = 75; % angle for ilium criterion

% sps computation
[sps1, sps2] = sps_DDH (img, angle);
                sps1_hist = reshape(sps1,[r*c,1]);[ML5,n]= hist(sps1_hist); ML5=normalize(ML5);
                sps2_hist = reshape(sps2,[r*c,1]);[ML6,n]= hist(sps2_hist); ML6=normalize(ML6);

% bone contour segmentation
[bone, probability, map] = bone_DDH (img, sps1, sps2);
                    probability_hist = reshape(probability,[r*c,1]);[ML7,n]= hist(probability_hist); ML7=normalize(ML7);
                    img_hist = reshape(img,[r*c,1]);[ML8,n]= hist(img_hist); ML8=normalize(ML8);
    
bone = bone .* (probability>0.5); % this threshold can be user controlled
                bone_hist = reshape(bone,[r*c,1]);[ML9,n]= hist(bone_hist); ML9=normalize(ML9);
  
% ray-cast to avoid the G muscles, may look similar to horizontal ilium    
bone_rc = ray_cast_DDH (bone);
                    bone_rc_hist = reshape(bone_rc,[r*c,1]);[ML10,n]= hist(bone_rc_hist); ML10=normalize(ML10);
ilium = bone_rc .* (sps2>0);
                    ilium_hist = reshape(ilium,[r*c,1]);[ML11,n]= hist(ilium_hist); ML11=normalize(ML11);

% alpha, beta, fhc = 0 if no ilium present                    
        if nnz(ilium)==0
            return        
        end
       
 %CO- added waitbar
waitbar(200/100,h)

        
        
% horizontal ilium criteria features - length and major angle
[hor_pared_ilium, hor_length, ver_length, major_angle] = horizontal_ilium (ilium, size_mm) ; 
        ML12 = nnz(hor_pared_ilium); ML13 = hor_length; ML14 = ver_length; ML15 = major_angle;
 
                    
% femoral head initialization
    femoral_head_radius = size_mm * 10; % 10mm at start, default
    

   
%% extract initial alpha, beta, fhc
[R, radon_strength_1, radon_strength_2, radon_strength_3, ilium_angle, acetab_angle, labrum_angle,...
    alpha, beta, bisector1, bisector2, center_FH, radius_FH, row_acetab, col_acetab]= ilium_acetabulum_labrum2 (hor_pared_ilium,sps1, bone_rc, femoral_head_radius) ;
                    R_hist = reshape(R,[size(R,1)*size(R,2),1]);[ML16,n]= hist(R_hist); ML16=normalize(ML16);
                    ML17=radon_strength_1; ML18=radon_strength_2; ML19=radon_strength_3; ML20=ilium_angle; ML21=acetab_angle;
                    ML22=labrum_angle; ML23=alpha; ML24=beta; ML25=bisector1; ML26=bisector2; ML27=center_FH/size_mm; 
                    ML28=radius_FH/size_mm; ML29=row_acetab/size_mm; ML30=col_acetab/size_mm;

%             %intermediate visualization
% 
%                 figure, 
%                 tmp = [bisector1;bisector2];
%                 show_straight_line(tmp,img) 
%                 viscircles(center_FH, radius_FH,'EdgeColor','b');   

% edge of ilium
center_alpha_beta = [col_acetab, row_acetab];
    

% line between acetabulum and ilium, for use in femoral head segmentation
angle_bet_ilium_acetabulum = (acetab_angle+labrum_angle) / 2;
x1 = col_acetab; y1=row_acetab;
m_F = tand(90-angle_bet_ilium_acetabulum);
c_F = -m_F * x1 + y1;



% Femoral head segmentation
%new 4th april 2017
[center_FH_new, radius_FH_new, quality_FH, quality_trirad] = femoral_head_feb_miccai2 (sps1,probability,map,radius_FH,center_FH,center_alpha_beta,alpha,beta,size_mm,img,m_F,c_F);


% % attempt using  fhc_ref = tand(alpha)/ (tand(alpha)+tand(beta)) *100;                    
% fhc_ref = tand(alpha)/ (tand(alpha)+tand(beta)) *100;  
% 
%  for i=1:30    
%     [center_FH_new, radius_FH_new, quality_FH, quality_trirad] = femoral_head_feb_miccai3 (sps1,probability,map,radius_FH,center_FH,center_alpha_beta,alpha,beta,size_mm,img);
%                     ML31=center_FH_new/size_mm; ML32=radius_FH_new/size_mm; ML33=quality_FH; ML49=quality_trirad;
% 
%     center_final = center_FH_new;
%     radius_final = radius_FH_new;
% 
%     %CO- added waitbar
%     waitbar(600/100,h)
% 
%     % fhc computation
%     edge_of_ilium = center_alpha_beta;
%     fhc = extract_fhc(ilium_angle, edge_of_ilium, radius_FH_new, center_FH_new); 
%     
%     if (abs(fhc-fhc_ref))<10
%         break
%     elseif (fhc-fhc_ref)<0
%         center_FH(1) = center_FH(1) +2;
%         center_FH(2) = center_FH(2) +2;
%         radius_FH = radius_FH - 5;
%     else % (fhc-fhc_ref)>0
%         center_FH(1) = center_FH(1) +2;
%         center_FH(2) = center_FH(2) -2;
%         radius_FH = radius_FH - 5;
%     end
%     
%  end

%% fhc computation
center_final = center_FH_new;
radius_final = radius_FH_new;

%CO- added waitbar
waitbar(600/100,h)



%%
ML31=center_FH_new/size_mm; ML32=radius_FH_new/size_mm; ML33=quality_FH; ML49=quality_trirad;

% % old 
% [center_FH_new, radius_FH_new, quality_FH, quality_trirad] = femoral_head_feb_miccai (sps1,probability,map,radius_FH,center_FH,center_alpha_beta,alpha,beta,size_mm,img);
%                     ML31=center_FH_new/size_mm; ML32=radius_FH_new/size_mm; ML33=quality_FH; ML49=quality_trirad;




waitbar(400/100,h)

% alpha, beta, fhc = 0 if no femoral head present 
    if radius_FH_new==0
            return        
        end

        
% angle values refinement
femoral_head_radius = radius_FH_new; 

        % visualization
           % figure, imagesc(img), colormap gray;axis equal
           % viscircles(center_FH_new, radius_FH_new,'EdgeColor','b');axis equal


           
           
[R, radon_strength_1, radon_strength_2, radon_strength_3, ilium_angle, acetab_angle, labrum_angle,...
    alpha, beta, bisector1, bisector2, center_FH_new, radius_FH_new, row_acetab, col_acetab]= ilium_acetabulum_labrum2 (hor_pared_ilium,sps1, bone_rc, femoral_head_radius) ;

                    R_hist = reshape(R,[size(R,1)*size(R,2),1]);[ML34,n]= hist(R_hist); ML34=normalize(ML34);
                    ML35=radon_strength_1; ML36=radon_strength_2; ML37=radon_strength_3; ML38=ilium_angle; ML39=acetab_angle;
                    ML40=labrum_angle; ML41=alpha; ML42=beta; ML43=bisector1; ML44=bisector2; ML45=center_FH/size_mm; 
                    ML46=radius_FH/size_mm; ML47=row_acetab/size_mm; ML48=col_acetab/size_mm;
alpha;
beta;
 
I = ilium_angle;
L = labrum_angle;
A=acetab_angle;
edge_r = row_acetab;
edge_c =col_acetab;

% combine all feature vectors    
feature(1,:)=[ML1,ML2,ML3,ML4,ML5,ML6,ML7,ML8,ML9,ML10,ML11,ML12,ML13,ML14,ML15,ML16,ML17,ML18,ML19,ML20,ML21,ML22,ML23,ML24...
    ,ML25,ML26,ML27,ML28,ML29,ML30,ML31,ML32,ML33,ML34,ML35,ML36,ML37,ML38,ML39,ML40,ML41,ML42,ML43,ML44,ML45,ML46,ML47,ML48];    
 


% fhc computation
edge_of_ilium = center_alpha_beta;

if length(center_FH_new)<2
            return        
end



fhc = extract_fhc(ilium_angle, edge_of_ilium, radius_FH_new, center_FH_new);







alpha_1=alpha;
beta_1=beta;
  


% % method 2 based angle calculation
%     
%     f_angles = [45 90 135];
%     minW=[20 20 20];
%     % 2D OPS
%     OPS= optimized_phasesym(img,5,3,f_angles,minW,1.25, 1.5,.1,1);%k=.5 here, change to 3 for higher noise threshold
%     [R, radon_strength_1, radon_strength_2, radon_strength_3, ilium_angle, acetab_angle, labrum_angle,...
%     alpha_OPS, beta_OPS, bisector1, bisector2, center_FH, radius_FH, row_acetab, col_acetab]= ilium_acetabulum_labrum2 (hor_pared_ilium,OPS, bone_rc, femoral_head_radius) ;
% 
%    
%    
% 
% alpha_2=alpha_OPS;
% beta_2=beta_OPS;
% %CO- added waitbar
% waitbar(800/100,h)



