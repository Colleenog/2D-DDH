function show_straight_line (straight_line, im, edge_c)

%im=zoomed_in;
%straight_line=lab;
%# read and display image
%im = imread('autumn.tif');
imagesc(im), axis equal, axis off, colormap gray
% 
% %# make sure the image doesn't disappear if we plot something else
 hold on
% 
% %# define points (in matrix coordinates)
% %p1 = [10,100];
% %p2 = [100,20];
% 
% %# plot the points.
% %# Note that depending on the definition of the points,
% %# you may have to swap x and y
% %plot([p1(2),p2(2)],[p1(1),p2(1)],'Color','r','LineWidth',2)


for iter = 1:size(straight_line,1)

    m =straight_line(iter,1); c=straight_line(iter,2);
    pp = [m c];%pp = [1 0]; % Polynomial as MATLAB likes it.
    %X = 1:size(im,1):size(im,1)+1;  % Where to plot the line...
    
    %Colleen- tried changed the end point of the line by ending where X is
    %evaluated. currently defined by row size of im, maybe we want the
    %column size?
    X = edge_c-40:edge_c+40;%2010:size(im,2)/2;  % Where to plot the line...
    % X = -200:200;  
    pv = polyval(pp,X);  % Evaluate the polynomial.
    plot(X,pv,'Color','y','LineWidth',3)
    
end
    
 %figure
 % plot(X,pv,'Color','g','LineWidth',3)
 %out=1;