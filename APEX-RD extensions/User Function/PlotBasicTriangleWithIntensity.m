function PlotBasicTriangleWithIntensity(RMats, I, options)
% PlotBasicTriangleWithIntensitiy: This function plots the location of the
% loading direction in the stereographic triangle for all input grain
% orientations and colors the points according to 'I' and a colormap
% 
% USAGE: PlotBasicTriangleWithIntensity(RMats, I, options)
%
% AUTHORS: Timothy Long, Robert Carson, and Mark O. (maybe Darren P.)
% 
% INPUTS:
%   RMats is 3 x 3 x n:
%       A list of rotation matricies (Rsc * c = s) describing the
%       orientation of n grains.
% 
%   I is n x 1:
%       A list of intensity values to color the points according to.
% 
%   options is a 1x1 structure:
%       Contains options that change how the code runs
% 
%       .colorMap is a 256 x 3
%           The color map to be used.  Default is jet.
% 
%       .scale is a 1x2
%           Gives the values of intensity that correspond to the minimum
%           and maximum colormap values.
% 
%       .monochrome is a logical
%           If set to true, 'I' is ignored and all points are colored red.
%           Default is false.
% 
%       .loadingDir is 1x3;
%           A vector giving the direction the external load is applied in.
%           The sign of the vector does not matter, both the positive and
%           negative directions will be used.  Default is [0 1 0].
%
% OUTPUTS:
%   N/A
% 
%
% NOTES:
%   Started on 2015_Aug_13
% 
%   This code only works for cubic materials
%
%   This code was written by T.L. utilizing code from PlotBasicTriangleQuat
%   by R.C. and M.O.  T.L. is solely responcible for any bugs or errors in
%   this code.


% set default options
cmap = jet(256);
monochrome = false; % if true, plot all points with the same color
scale = [min(I), max(I)];
loadingDir = [0 1 0];
plotGrNums = 0;


% overwrite default options if passed in
if(exist('options','var'))
    
    if isfield(options,'colorMap')
        cmap = options.colorMap;
    end
    
    if isfield(options,'scale')
        scale = options.scale;
    end
    
    if isfield(options,'monochrome')
        monochrome = options.monochrome;
    end
    
    if isfield(options,'loadingDir')
        loadingDir = options.loadingDir;
    end
    
end

numGrains = length(I);

% set colormap
c = colormap(cmap);
hold on;
axis equal off

% map from intensity values to colormap
if monochrome || isempty(I)
    colors = repmat([1 0 0],256,1);
else
    index4 = floor(255/(scale(2)-scale(1))*(I-scale(1)))+1;
    index4(index4>255) = 255;
    index4(index4<1) = 1;
    colors = c( index4,: );
end

% include the negative loading direction.  Also make the loading
% direction a column vector
if(size(loadingDir,1)==1)
    loadingDir = [loadingDir', -1*loadingDir'];
else
    loadingDir = [loadingDir, -1*loadingDir];
end


% Cubic crystal symmetries
RMatSym = RMatOfQuat(CubSymmetries);
nsym = size(RMatSym,3);



% find the direction of the loading direction in the crystal coorindate
% system for each grain.
for ii = 1:numGrains
    
    % calculate the loading direction in the crystal coordinate system
    loadingDirC = RMats(:,:,ii)' * loadingDir;
    
    % initialize storage variables
    lDirsC_1 = zeros(3,nsym);
    lDirsC_2 = zeros(3,nsym);
    
    % calculate all symmetrically equivilant loading directions
    for jj = 1:nsym
        lDirsC_1(:,jj) = RMatSym(:,:,jj) * loadingDirC(:,1);
        lDirsC_2(:,jj) = RMatSym(:,:,jj) * loadingDirC(:,2);     
    end
    
    
    % Test +LD
    x = lDirsC_1(1,:);
    y = lDirsC_1(2,:);
    z = lDirsC_1(3,:);
    Y1 = y./(1+x);
    Z1 = z./(1+x);
    
    % Find the point in the fundamental stereographic triangle
    index1 = find( ((Y1 > Z1) & (Y1 > 0) & (Y1 < sqrt(2)/(2+sqrt(2))) & (Z1 > 0)),1 );
    
    
    % Test -LD
    x = lDirsC_2(1,:);
    y = lDirsC_2(2,:);
    z = lDirsC_2(3,:);
    Y2 = y./(1+x);
    Z2 = z./(1+x);
    
    % Find the point in the fundamental stereographic triangle
    index2 = find( ((Y2 > Z2) & (Y2 > 0) & (Y2 < sqrt(2)/(2+sqrt(2))) & (Z2 > 0)),1 );
    
    
    % Plot the point
    if ~isempty(index1)
        inside=insideTriangle([Y1(index1),Z1(index1)],index1);
        plot(Y1(inside),Z1(inside),'o','Color','b','MarkerFaceColor',colors(ii,:),'MarkerSize',6)
        if plotGrNums == 1
            text(Y1(inside)+0.005,Z1(inside), num2str(ii), 'Color', 'k');
        end
    elseif ~isempty(index2)
        inside=insideTriangle([Y2(index2),Z2(index2)],index2);
        plot(Y2(inside),Z2(inside),'o','Color','b','MarkerFaceColor',colors(ii,:),'MarkerSize',6)
        if plotGrNums == 1
            text(Y2(inside)+0.005,Z2(inside), num2str(ii), 'Color', 'k');
        end
    end
    
    
end

% add colorbar with scale
set(gca,'CLim',scale)
colorbar('Location','EastOutside')
    

% set font type
set(gca,'FontSize',18,'FontName','Times')

% Specify theta and phi range
% theta is the azimuthal angle (rotate counterclockwise from x about z-axis)
% phi is the elevation angle from the x-y plane
thetaGrid = [0:0.1:45]*pi/180;
phiGrid   = [0:0.1:45]*pi/180;
[theta, phi] = meshgrid(thetaGrid, phiGrid);
theta_array  = reshape(theta, 1, length(thetaGrid)^2);
phi_array    = reshape(phi, 1, length(thetaGrid)^2);


% Plot outlines of basic orientation triangle
% 45 degree line
index3 = (abs(tan(phi_array) - sin(theta_array)) < 0.0006);
[x_line, y_line, z_line] = sph2cart(theta_array(index3), phi_array(index3), ones(1,length(theta_array(index3))));
Y_line = y_line./(1+x_line);
Z_line = z_line./(1+x_line);
plot(Y_line, Z_line, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 2)

% bottom, horizontal line
theta_line = [0:0.3:45]*pi/180;
phi_line   = zeros(1,length(theta_line));
[x_line, y_line, z_line] = sph2cart(theta_line, phi_line, ones(1,length(theta_line)));
Y_line = y_line./(1+x_line);
Z_line = z_line./(1+x_line);
plot(Y_line, Z_line, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 3)

% right, cemicircular line
[~,p,~] = cart2sph(1,1,1);
phi_line   = [0:0.1:p*180/pi]*pi/180;
theta_line = ones(1,length(phi_line))*45*pi/180;
[x_line, y_line, z_line] = sph2cart(theta_line, phi_line, ones(1,length(phi_line)));
Y_line = y_line./(1+x_line);
Z_line = z_line./(1+x_line);
plot(Y_line, Z_line, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 2)

hold off

end


% helper function to ensure the points are inside the triangle
function inside=insideTriangle(points,index)

    inside=[];

    [~,p,~] = cart2sph(1,1,1);
    phi_line   = [0:0.1:p*180/pi]*pi/180;
    theta_line = ones(1,length(phi_line))*45*pi/180;
    [x_line, y_line, z_line] = sph2cart(theta_line, phi_line, ones(1,length(phi_line)));
    Y_line = y_line./(1+x_line);
    Z_line = z_line./(1+x_line);

    Z_new=Z_line;
    Y_new=Y_line;
    
    for i=1:length(points(:,1))
        
        point=points(i,:);
        
        if(all(point(1)<Y_new))
            
            inside=[inside,index(i)];
            
        else
        
            diffY=Y_new-point(1);
            diffZ=Z_new-point(2);
            
            if(any(diffZ(diffY<=0)<0))
                continue;
            else
                inside=[inside,index(i)];
            end
            
        end
        
    end
end

