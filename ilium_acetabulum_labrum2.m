function [R, radon_strength_1, radon_strength_2, radon_strength_3, ilium_angle, acetab_angle, labrum_angle,...
    alpha, beta, bisector1, bisector2, center_FH, radius_FH, row_acetab, col_acetab] ...
     = ilium_acetabulum_labrum2 (hor_pared_ilium,sps1, bone_rc, femoral_head_radius) 

% size_mm = round (size(sps1,2)/40); %number of pixels per mm
% femoral_head_radius = size_mm * 10 % 10mm at start, default

% find column number of rightmost point of ilium
    temp = sum(hor_pared_ilium);
    col_acetab = max(find(temp));

    temp = hor_pared_ilium(:,col_acetab);
    row_acetab = max(find(temp));
    
    [r, c] =size(sps1);

% ROI of region around edge of ilium 
    sps1_ROI = spatial_filter_circle (sps1, femoral_head_radius, col_acetab, row_acetab);
    
    added_rows = zeros(femoral_head_radius, c);        % add rows and columns of zeros so that index in zoomed_in doesn't go to negative
    added_columns = zeros(r+2*femoral_head_radius, femoral_head_radius);
    sps1_ROI_extended = vertcat (added_rows,sps1_ROI,added_rows);
    sps1_ROI_extended = horzcat (added_columns,sps1_ROI_extended,added_columns);
    
    zoomed_in = sps1_ROI_extended ( row_acetab:row_acetab+2*femoral_head_radius, col_acetab:col_acetab+2*femoral_head_radius); % to center it
    zoomed_in=log(1+zoomed_in*50);
    %figure, imagesc(zoomed_in)
    
    zoomed_in_ilium = zoomed_in;
    zoomed_in_ilium(:, femoral_head_radius:end)=0; %1st and 4th quadrant 0
   % figure, imagesc(zoomed_in_ilium)
     
    zoomed_in_acetab = zoomed_in;
    zoomed_in_acetab(:, 1:femoral_head_radius)=0; %1st 2nd and 3rd quadrant 0
    zoomed_in_acetab(1:femoral_head_radius, femoral_head_radius:end)=0;
%     figure, imagesc(zoomed_in_acetab)    
    
    zoomed_in_labrum = zoomed_in;
    zoomed_in_labrum(:, 1:femoral_head_radius-20)=0; % 2nd 3rd and 4th quadrant 0
    zoomed_in_labrum(femoral_head_radius:end, :)=0; 
   % zoomed_in_labrum(1:10, :)=0; 
    %zoomed_in_labrum(1:round(femoral_head_radius/3), femoral_head_radius:end)=0;
    
    
%ilium angle
    theta = 85:95;
    [R,xp] = radon(zoomed_in_ilium,theta);
    R = medfilt2(R, [5 5]); % to remove noise/interference from integration of straight lines across a wider ROI than needed

    lower_bound = round(length(xp)/2-femoral_head_radius/4); % all lines cross within a circle of radius 2mm
    upper_bound = round(length(xp)/2+femoral_head_radius/4);

    zoomed_in_radon = R(lower_bound:upper_bound, :);
    [radon_strength_1, angle_ilium] = max(max(zoomed_in_radon));
    ilium_angle = 85+angle_ilium;
   % figure, imagesc(zoomed_in_radon)

%                     imagesc(theta,xp,R);
%                     title('R_{\theta} (X\prime)', 'fontsize', 30);
%                     xlabel('\theta (degrees)', 'fontsize', 30);
%                     ylabel('X\prime', 'fontsize', 30);
%                     set(gca,'XTick',85:5:95, 'fontsize', 30);
%                     colormap(hot);
%                     colorbar ('fontsize', 30)
                    
                    
%acetabulum angle
    theta = [ilium_angle-80:ilium_angle-25];
    [R,xp] = radon(zoomed_in_acetab,theta);
    R = medfilt2(R, [5 5]); % to remove noise/interference from integration of straight lines across a wider ROI than needed
    
    R(:,1:10)=R(:,1:10).*1.7;
    R(:,11:15)=R(:,11:15).*1.6;
    R(:,16:20)=R(:,16:20).*1.5;
    R(:,21:30)=R(:,21:30)*1.35; %cant make poor hip look good, but can make good hip look poor
    R(:,31:40)=R(:,31:40)*1.2; 
    
    lower_bound = round(length(xp)/2-femoral_head_radius/2); % all lines cross within a circle of radius 4mm for acetabulum
    upper_bound = round(length(xp)/2+femoral_head_radius/2);

    zoomed_in_radon = R(lower_bound:upper_bound, :);
    [radon_strength_2, index] = max(max(zoomed_in_radon));
    acetab_angle = theta(index);
   % figure, imagesc(zoomed_in_radon)

%                     imagesc(theta,xp,R);
%                     title('R_{\theta} (X\prime)', 'fontsize', 30);
%                     xlabel('\theta (degrees)', 'fontsize', 30);
%                     ylabel('X\prime', 'fontsize', 30);
%                     set(gca,'XTick',10:30:70, 'fontsize', 30);
%                     colormap(hot);
%                     colorbar ('fontsize', 30)
                    
    
%labrum angle
    theta = ilium_angle+15:ilium_angle+80;
    [R,xp] = radon(zoomed_in_labrum,theta);
    R = medfilt2(R, [5 5]); % to remove noise/interference from integration of straight lines across a wider ROI than needed

    R(:,15:30)=R(:,15:30)*3.5;
    R(:,31:40)=R(:,31:40)*2.5;
    R(:,55:end)=R(:,55:end)*4.3;

    
    lower_bound = round(length(xp)/2-femoral_head_radius/2); % all lines cross within a circle of radius 2mm
    upper_bound = round(length(xp)/2+femoral_head_radius/2);

    zoomed_in_radon = R(lower_bound:upper_bound, :);
    [radon_strength_3, index] = max(max(zoomed_in_radon));
    labrum_angle = theta(index);
    %figure, imagesc(zoomed_in_radon)

%                     imagesc(theta,xp,R);
%                     title('R_{\theta} (X\prime)', 'fontsize', 30);
%                     xlabel('\theta (degrees)', 'fontsize', 30);
%                     ylabel('X\prime', 'fontsize', 30);
%                     set(gca,'XTick',110:30:170, 'fontsize', 30);
%                     colormap(hot);
%                     colorbar ('fontsize', 30)

% angle calculation
    alpha = ilium_angle - acetab_angle;
    beta = labrum_angle - ilium_angle;






%% bisectors

    [bisector1, bisector2, center_FH, radius_FH] = bisector (col_acetab,row_acetab,femoral_head_radius,acetab_angle,labrum_angle);

    tmp = [bisector1;bisector2];
   % figure, show_straight_line (tmp, sps1)
