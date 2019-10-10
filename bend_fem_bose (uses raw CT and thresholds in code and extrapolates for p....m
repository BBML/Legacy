function [ ]  = bend_fem_Bose_preload();

tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE!!!! 9/27/19 THIS CODE CAN BE USED BY MEMBERS OF THE BBML
% IT READS IN RAW FILES FROM CT AND THRESHOLDS INTERNALLY
% 
% Code is for femur and now allows this to work for 3 or 4 point bending 
% which is hard coded.
%
% AGB adapted on 7/24/15 to not zero the load and displacement when you choose
% the start point due to problems with rolling during testing. Instead, the
% program will use this point, then perform a linear regression to take
% this back to 0,0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Written by Joey Wallace, July 2012 (originally r4tib.m)

% Modified in August 2013 to work with femurs tested in 3 or 4 point bending

% Edited by Max Hammond Sept. 2014 Changed the output from a csv 
% file to an xls spreadsheet that included a title row. Code written by
% Alycia Berman was added into the CTgeom section of the code to subtract
% out the scale bars that appear in some CT images. Used while loop to
% semi-batch process. Hard coded in initial values like bendtype, slice
% number, and voxel size because each will be held constant within a study.
% Added the option to smooth or not during Testing Configuration. Smoothed
% using a moving average with a span of 10. Added a menu in case points
% need to be reselected.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a comprehensive program that reads in geometric and mechanical
% information to calculate force/displacement and stress/strain from 
% BENDING MECHANICAL TESTS IN THE FEMUR ONLY (3 OR 4 POINT)

% The program reads raw image files from Skyscan CT BMP ROIs to calculate
% geometric parameters.  As written, individual raw BMP images will be read
% in simultaneously and parameters will be calculated from each.

% Geometric files should be named "ID#_RF_slice_#.bmp" or "ID#_LF_slice_#.bmp".
% An example is 107_RF_slice_1.bmp.

% This program reads raw mechanical data from the Bose system
% from a a file named "specimen name.csv". It assumes that mechanical specimen
% names are written as "ID#_RF" or "ID#_LF". For femora, the assumption is that
% bending was about the ML axis with the anterior surface in tension.
%  
% The program adjusts for system compliance and then uses beam bending
% theory to convert force-displacement data to theoretical stress-strain
% values.  Mechanical properties are calculated in both domains and output
% to a file "specimen name_mechanics.csv".  It also outputs a figure
% showing the load-displacement and stress-strain curves with significant
% points marked.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%close all figure windows and clears all variables
close all
clear all
dbstop if error

%*****************\TESTING CONFIGURATION/**********************************
%                                                                         *
%       Adjust these values to match the system setup                     *
%                                                                         *
L = 9.00;           %span between bottom load points (mm)                 *
a = 3.00;           %distance between outer and inner points (if 4pt; mm) *
bendtype = '4';     %enter '4' for 4pt and '3' for 3pt bending            *
compliance = 0;     %system compliance (microns/N)                        *
side = 'RF';        %input 'LF' for left and 'RF' for right feumr         *
slices = 7;         %number of slices from uCT to read in                 *
res = 8.40781;      %voxel resolution in um of uCT                        *
ang = 0.5;          %angle step in degrees                                *
smoothing = 1;      %enter 1 to smooth using moving average (span=10)     *
threshold = 95;     %enter grayscale threshold from CTAn                  *
%**************************************************************************

%Check common errors in testing configuration
if bendtype ~= '3' && bendtype ~= '4'
        error('Please enter 3 or 4 for bendtype as a string in the Testing Configuration')
end

if strcmp(side,'LF') == 0 && strcmp(side,'RF') == 0
        error('Please enter LF or RF for side as a string in the Testing Configuration')
end

if smoothing ~= 1 && smoothing ~= 0
        error('Please enter a 1 or 0 for smoothing in the Testing Configuration')
end


%create a while loop to quickly run through multiple files without running
%the program more than once
zzz=1;
ppp=1;

while zzz==1

%input bone number
number = input('Bone Number: ','s');

%create ID from bone and side inputs
side_type = [side];
specimen_name = [number '_' side];
ID = [specimen_name];

%setting up arrays for data output later
A = 360/ang+1;
ang_out = [45:ang:405];
peri_out = zeros(slices,A);
endo_out = zeros(slices,A);
profiles = zeros((2*slices+3),A);
geom_out = zeros(slices,2);
aveprof = zeros(3,A);
geom_ave = zeros(1,2);

%Now setting up a loop for the actual calculations;
for j=1:slices
    close all   
    
    %Read in CT ROI BMPs
    layer = num2str(j);
    section(:,:,j) = imread([number '_' side_type '_slice_' layer '.bmp']);
    
    slice=section(:,:,j);
    
    slice = im2bw(slice,(threshold-1)/255); % use threshold-1 to take everything greater than or equal to the input threshold like CTAn
    
    %----------------------------Alycia-----------------------------------
    %Find the connected components in the image
    cc=bwconncomp(slice); 

    %Find the number of pixels in each component
    numPixels = cellfun(@numel,cc.PixelIdxList);

    %Find the index containing the most number of pixels
    [~,idx] = max(numPixels);

    %Remove all other components
    for i=1:length(cc.PixelIdxList)
        if i==idx
            %Do nothing
        else
            slice(cc.PixelIdxList{i})=0;
        end
    end
    
    clear idx
    clear cc
    clear numPixels
    %---------------------------------------------------------------------

    %create a box around the bone.  The 4 comnponents of this are 1) the x
    %coordinate of the UL corner 2) the y coordinate of the UL corner 3) the x
    %width of the box and 4) the y width of the box
    bb = regionprops(slice,'BoundingBox');
    bb = cat(1,bb.BoundingBox);

    %One problem is that the UL corner points in bb1 will actually be at the top left
    %corner of the pixel they represent, so we need to add 0.5 to the first and 
    %second entries to move these points to the center of the pixel.
    %Then, the widths will be 1 pixel long so we need to subtract 1 from
    %entries 3 and 4.  This will now represent the first on and last on pixel
    %in each direction

    bb(1)=bb(1)+0.5;
    bb(2)=bb(2)+0.5;
    bb(3)=bb(3)-1;
    bb(4)=bb(4)-1;

    %create a vector containing the coordinates of the corners of the bounding
    %box, going from UL, UR, LR, LL
    xbb = [bb(1),bb(1)+bb(3),bb(1)+bb(3),bb(1)];
    ybb = [bb(2),bb(2),bb(2)+bb(4),bb(2)+bb(4)];

    %Manually calculate the centroid from pixel locations:
    [index_y,index_x] = find (slice == 1); %this finds the x and y locations of each "on" pixel
    Qx = sum(index_y); %since the area of each dA is 1 pixel, we have x1 + x2 +...+ and this is the integral of y_dA
    Qy = sum(index_x); %since the area of each dA is 1 pixel, this is the integral of x_dA 
    area = length(index_y);     %since the area of each pixel is one, this is the total number of pixels or the area
    xbar = Qy/area; %xbar = integral of y_da/A
    ybar = Qx/area; %ybar = integral of x_da/A


    %Manually calculate MOI from pixel locations.  We need integral of x 
    %squared dA.  since each x location is in index_x,and since the area of of 
    %each is one, we need x1x1 + x2x2, so multiply the vectors and then sum 
    %the entries in the resulting vector.
    Ix  = sum(index_y.^2);       %Ix  = integral(y^2*dA); with respect to y-axis
    Iy  = sum(index_x.^2);       %Iy  = integral(x^2*dA); with respect to x-axis
    Ixy = sum(index_x.*index_y);  %Ixy = integral(xy *dA); with respect to xy-axis

    Ixc  = Ix  - area*ybar^2;      %Ixc  = Ix with respect to centroid; parallel axis theorema
    Iyc  = Iy  - area*xbar^2;      %Ixc  = Ix with respect to centroid; parallel axis theorema
    Ixyc = Ixy - area*xbar*ybar;     %Ixyc = Ixy with respect to centroid; parallel axis theorema

    %PLOT THE IMAGE WITH THE CENTROID MARKED
    subplot(2,2,1), imagesc(slice) %shows the image in blue vs red with axes labeled
    axis equal
    axis tight
    hold on
    plot(xbar, ybar, 'r+')
    xlabel('original x-positions in voxels')
    ylabel('original y-positions in voxels')
    title(['Original Voxel-Based BMP'])

    %start getting line profiles at various degrees:
    inner_fiber = []; %creats a blank vector for the endocortical radii
    outer_fiber = []; %creates a blank vector for the periosteal radii
    thickness = [];  %creates a blank vector for cortical thicknesses

    image_size = size(slice);
    x_size = image_size(2);
    y_size = image_size(1);

    %IN QUADRANT 1 (45 to 134.99 degrees, M quadrant):
    for i = -45:ang:44.9
        angle = i * pi / 180;
        yi=[ybar,0];
        xi=[xbar,xbar-(tan(angle)*yi(1))];
        [cx,cy,c] = improfile(slice,xi,yi,10000); 
            %improfiel chooses an arbirtray number of points to look at so I will choose alot to be accurate.
            %cx and cy are the pixel locations along the line and c is the intensity at each point
        cort_on = find(c == 1); 
        thick = length(cort_on);
        on_1 = cort_on(1);
        on_end = cort_on(thick);
        rad_x = [cx(on_1),cx(on_end)];
        rad_y = [cy(on_1),cy(on_end)];
        points=[xbar,ybar;rad_x(1),rad_y(1);rad_x(2),rad_y(2)]; %centroid, endo and peri points along this line
        radii = pdist(points); %radii(1) is endo, radii(2) is peri and radii(3) is c_thickness
        inner_fiber = [inner_fiber radii(1)];
        outer_fiber = [outer_fiber radii(2)];
        thickness = [thickness radii(3)];
    end

    %IN QUADRANT 2 (135 to 224.99 degrees, P quadrant):
    for i = -45:ang:44.9
        angle = i * pi / 180;
        xi=[xbar,0];  
        yi=[ybar,ybar+(tan(angle)*xi(1))];
        [cx,cy,c] = improfile(slice,xi,yi,10000); 
        cort_on = find(c == 1); 
        thick = length(cort_on);
        on_1 = cort_on(1);
        on_end = cort_on(thick);
        rad_x = [cx(on_1),cx(on_end)];
        rad_y = [cy(on_1),cy(on_end)];
        points=[xbar,ybar;rad_x(1),rad_y(1);rad_x(2),rad_y(2)]; 
        radii = pdist(points); 
        inner_fiber = [inner_fiber radii(1)];
        outer_fiber = [outer_fiber radii(2)];
        thickness = [thickness radii(3)];
    end

    %IN QUADRANT 3 (225 to 314.99 degrees, L quadrant):
    for i = -45:ang:44.9
        angle = i * pi / 180;
        yi=[ybar,y_size];
        xi=[xbar,xbar+(tan(angle)*(yi(2)-yi(1)))];
        [cx,cy,c] = improfile(slice,xi,yi,10000);        
        cort_on = find(c == 1); 
        thick = length(cort_on);
        on_1 = cort_on(1);
        on_end = cort_on(thick);
        rad_x = [cx(on_1),cx(on_end)];
        rad_y = [cy(on_1),cy(on_end)];
        points=[xbar,ybar;rad_x(1),rad_y(1);rad_x(2),rad_y(2)]; 
        radii = pdist(points); 
        inner_fiber = [inner_fiber radii(1)];
        outer_fiber = [outer_fiber radii(2)];
        thickness = [thickness radii(3)];
    end

    %IN QUADRANT 4 (315 to 404.99 or 49.99, A quadrant):
    for i = -45:ang:44.9
        angle = i * pi / 180;
        xi=[xbar,x_size];  
        yi=[ybar,ybar-(tan(angle)*(xi(2)-xi(1)))];
        [cx,cy,c] = improfile(slice,xi,yi,10000);
        cort_on = find(c == 1); 
        thick = length(cort_on);
        on_1 = cort_on(1);
        on_end = cort_on(thick);
        rad_x = [cx(on_1),cx(on_end)];
        rad_y = [cy(on_1),cy(on_end)];
        points=[xbar,ybar;rad_x(1),rad_y(1);rad_x(2),rad_y(2)]; 
        radii = pdist(points);     
        inner_fiber = [inner_fiber radii(1)];
        outer_fiber = [outer_fiber radii(2)];
        thickness = [thickness radii(3)];
    end

    %To plot this in polar,you need to append the inner and outer vectors with
    %the vaule at 360 degrees (0 deg) to close 
    inner_fiber = [inner_fiber inner_fiber(1)];
    outer_fiber = [outer_fiber outer_fiber(1)];

    %plots as a subplot in polar
    angle_deg = [45:ang:405];
    angle_rad = angle_deg.*pi./180; %convert to radians
    subplot(2,2,2)
    polar(360,50); %set  the axes for the polr plot
    hold on
    polar(angle_rad,inner_fiber);
    polar(angle_rad,outer_fiber);
    title(['Polar Plot of Original Voxel-Based BMP'])

    %Convert the geometric data from angle and radius to x and y coordinates
    outer_fiber_x = outer_fiber.*cos(angle_rad);
    outer_fiber_y = outer_fiber.*sin(angle_rad);
    inner_fiber_x = inner_fiber.*cos(angle_rad);
    inner_fiber_y = inner_fiber.*sin(angle_rad);

    %Before shifting the origin from the centroid, calculate extreme fiber in 
    %each anatomic direction.  A and P are not dependent on whether this is
    %a right or left bone, but M and L are so be careful:
    anterior_extreme = abs(max(outer_fiber_x));
    posterior_extreme = abs(min(outer_fiber_x));
    
    if side == 'RF'
        medial_extreme = abs(max(outer_fiber_y));
        lateral_extreme = abs(min(outer_fiber_y));
    else
        medial_extreme = abs(min(outer_fiber_y));
        lateral_extreme = abs(max(outer_fiber_y));
    end
    
    %Shift the coordinate system from (0,0) at centroid to the (0,0) at LL
    %corner.  For geometric properties, the outer perimeter needs to go in the CW
    %direction.  Currently, it is CCW so it needs to be flipped.  This does
    %both and plots to verify
    x_data_min = abs(min(outer_fiber_x));
    y_data_min = abs(min(outer_fiber_y));
    outer_fiber_x = outer_fiber_x+x_data_min;
    outer_fiber_x = fliplr(outer_fiber_x);
    outer_fiber_y = outer_fiber_y+y_data_min;
    outer_fiber_y = fliplr(outer_fiber_y);
    inner_fiber_x = inner_fiber_x+x_data_min;
    inner_fiber_y = inner_fiber_y+y_data_min;
    x_perimeter = [outer_fiber_x inner_fiber_x];
    y_perimeter = [outer_fiber_y inner_fiber_y];
    subplot(2,2,3)
    plot(x_perimeter,y_perimeter)
    hold on
    axis equal
    axis tight
    plot(x_data_min,y_data_min,'g*')
    xlabel('x-position in voxels')
    ylabel('y-position in voxels')
    title(['Calculated Perimeter Points'])

   %Convert all pixel values to um and plot
    outer_fiber_x = outer_fiber_x*res;
    outer_fiber_y = outer_fiber_y*res;
    inner_fiber_x = inner_fiber_x*res;
    inner_fiber_y = inner_fiber_y*res;
    x_perimeter = x_perimeter*res;
    y_perimeter = y_perimeter*res;

    subplot(2,2,4)
    plot(x_perimeter,y_perimeter)
    hold on
    axis equal
    axis tight
    xlabel('x-position in um')
    ylabel('y-position in um')
    title(['Calculated Perimeter Points'])

    %Calculate extreme fiber in anterior direction in um
    anterior_extreme = anterior_extreme * res;
    
    %*****NOW INCORPORATE POLYGEOM TO GET MOI - ALL INPUTS IN UM*****
    % all of the ouptuts from this part are already in um
    clear x y
    x = x_perimeter;
    y = y_perimeter;

    % check if inputs are same size
    if ~isequal( size(x), size(y) ),
      error( 'X and Y must be the same size');
    end

    % number of vertices
    [ x, ns ] = shiftdim( x );
    [ y, ns ] = shiftdim( y );
    [ n, c ] = size( x );

    % temporarily shift data to mean of vertices for improved accuracy
    xm = mean(x);
    ym = mean(y);
    x = x - xm*ones(n,1);
    y = y - ym*ones(n,1);

    % delta x and delta y
    dx = x( [ 2:n 1 ] ) - x;
    dy = y( [ 2:n 1 ] ) - y;

    % summations for CW boundary integrals
    A = sum( y.*dx - x.*dy )/2; %cortical area
    Axc = sum( 6*x.*y.*dx -3*x.*x.*dy +3*y.*dx.*dx +dx.*dx.*dy )/12; %First moment about the y-axis (xc*A)
    Ayc = sum( 3*y.*y.*dx -6*x.*y.*dy -3*x.*dy.*dy -dx.*dy.*dy )/12; %First moment about the x-axis (yc*A)
    Ixx = sum( 2*y.*y.*y.*dx -6*x.*y.*y.*dy -6*x.*y.*dy.*dy ...
              -2*x.*dy.*dy.*dy -2*y.*dx.*dy.*dy -dx.*dy.*dy.*dy )/12;%second moment about x axis
    Iyy = sum( 6*x.*x.*y.*dx -2*x.*x.*x.*dy +6*x.*y.*dx.*dx ...
              +2*y.*dx.*dx.*dx +2*x.*dx.*dx.*dy +dx.*dx.*dx.*dy )/12;%second moment about y axis
    Ixy = sum( 6*x.*y.*y.*dx -6*x.*x.*y.*dy +3*y.*y.*dx.*dx ...
              -3*x.*x.*dy.*dy +2*y.*dx.*dx.*dy -2*x.*dx.*dy.*dy )/24;%product of inertia about x-y axes
    P = sum( sqrt( dx.*dx +dy.*dy ) ); %Perimeter?

    % check for CCW versus CW boundary
    if A < 0,
      A = -A;
      Axc = -Axc;
      Ayc = -Ayc;
      Ixx = -Ixx;
      Iyy = -Iyy;
      Ixy = -Ixy;
    end

    % centroidal moments
    xc = Axc / A; %centroidal location in x direction
    yc = Ayc / A; %centroidal location in y direction
    Iuu = Ixx - A*yc*yc; %centroidal MOI about x axis
    Ivv = Iyy - A*xc*xc; %centroidal MOI anout y axis
    Iuv = Ixy - A*xc*yc; %Product of inertia
    J = Iuu + Ivv; %Polar MOI

    % replace mean of vertices
    x_cen = xc + xm;
    y_cen = yc + ym;
    Ixx = Iuu + A*y_cen*y_cen;
    Iyy = Ivv + A*x_cen*x_cen;
    Ixy = Iuv + A*x_cen*y_cen;

    %plot the centoid output from Polygeom on last graph
    plot(x_cen,y_cen,'b+')

    %converts the centroid back to pixels and overlays on 3rd plot to show
    %difference from this measure and original measure
    cen_x_pix = x_cen/res;
    cen_y_pix = y_cen/res;
    subplot(2,2,3)
    plot(cen_x_pix,cen_y_pix,'b+')

    geometry = [Ivv anterior_extreme];
    geom_out(j,:) = geometry;
        
end

% Save last figure as an image
geom_name = [number '_' side '_geom'];
print ('-dpng', geom_name)
close

% Calculate average c and MOI from the slices and convert to proper units 
geom_ave (1,:) = mean(geom_out(1:slices,:));
Ivv = geom_ave(1);              %Centroidal MOI about the ML axis in um^4
Ivv = Ivv * 1e-12;             %Centroidal MOI about the ML axis in mm^4     
I = Ivv;                       %Centroidal Moment of Inertia about ML axis (mm^4)

anterior_extreme = geom_ave(2);   %in um
c = anterior_extreme;            %in microns 

%Set up loop to redo point selection if desired
yyy = 1;
while yyy==1    

%Read in raw mechanical testing data from Bose system
imported_data = csvread([ID '.csv'],5,0);
load = imported_data (:,4);         %in N
load = load*(-1);
position = imported_data (:,3);     %in mm
position = position * 10^3 * (-1);  %in microns

%Moving average smoothing with span of 10 if selected initially
if smoothing == 1
    load = smooth(load,10,'moving');
end

%Plot raw data and choose failure point
figure (1)
plot (position,load)
xlabel ('Displacement (microns)')
ylabel ('Force (N)')

title ('Pick failure point:')
[x,y]=ginput(1);
hold on
plot(x,y,'ro')
close

%Find the failure point index
i=1;
while position(i) < x
    i=i+1;
end

%Truncate the data at the failure point
position = position(1:i);
load = load(1:i);

%Plot truncated data and choose starting point
figure (1)
plot (position,load)
xlabel ('Displacement (microns)')
ylabel ('Force (N)')

title ('Pick beginning point:')
[x,y]=ginput(1);
hold on
plot(x,y,'ro')
close

%Find the start point index
j=1;
while position(j) < x
    j=j+1;
end

%Truncate the data at the start point, but do not set this to zero
position = position(j:i-1);
load = load(j:i-1);
displacement = position - load*compliance;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Convert the corrected load/displacement data to stress/strain
if bendtype == '3'
    stress = (load*L*c) / (4*I) * 10^-3;             %MPa
    strain = (12*c*displacement) / (L^2);             %microstrain
end

if bendtype == '4'
   stress = (load*a*c) / (2*I) * 10^-3;             %MPa
   strain = (6*c*displacement) / (a*(3*L - 4*a));   %microstrain
end

%Plot the adjusted stress-strain curve and pick points to define modulus
figure (2)
plot(strain,stress)
axis xy
xlabel('Strain (microstrain)')
ylabel('Stress (MPa)')
title('Pick points to define modulus:')
[x,y] = ginput(2);
hold on
plot(x,y,'ro')

if x(2) < x(1)
    error ('The 2nd selected point must have a larger strain than the 1st selected point')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alycia Berman added the load and displacement portion of the code below 
% on 7/24/15 so that the data does not have to be zeroed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Create stress and strain vectors spanning region for modulus determination
i=1;
while x(1) > strain(i)
    i = i+1;
end
linear_strain(1) = strain(i);
linear_stress(1) = stress(i);

linear_displacement(1) = displacement(i);
linear_load(1) = load(i);

j=2;
while x(2) > strain(i)
    linear_strain(j) = strain(i);
    linear_stress(j) = stress(i);
    
    linear_displacement(j) = displacement(i);
    linear_load(j) = load(i);
    
    i = i+1;
    j = j+1;
end

plot(linear_strain,linear_stress,'r')


%Determine modulus by linear regression of selected points
coeff = polyfit(linear_strain,linear_stress,1);
slope = coeff(1);
modulus = slope * 10^3;                                 % GPa

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if bendtype == '3'
    stiffness = modulus*48*I / (L^3) * 10^3;   % N/mm
end

if bendtype == '4'
   stiffness = modulus*12*I / (a^2 * (3*L -4*a)) * 10^3;   % N/mm
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p=polyfit(linear_displacement,linear_load,1);

x_shift=-p(2)/p(1);
displacement=displacement-x_shift;

disp_extension=0:0.01:displacement(1);
load_extension=disp_extension.*p(1);

displacement=[disp_extension'; displacement];
load=[load_extension';load];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if bendtype == '3'
    stress = (load*L*c) / (4*I) * 10^-3;             %MPa
    strain = (12*c*displacement) / (L^2);             %microstrain
end

if bendtype == '4'
   stress = (load*a*c) / (2*I) * 10^-3;             %MPa
   strain = (6*c*displacement) / (a*(3*L - 4*a));   %microstrain
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Plot linear regression on top of selected region
hold on
linear_stress = slope*linear_strain;
plot(linear_strain,linear_stress,'g')

% Create line with a .2% offset (2000 microstrain)
y_int = -slope*2000;        %y intercept
y_offset = slope*strain + y_int;    %y coordinates of offest line


%Find indeces where the line crosses the x-axis and the stres-strain curve.
%Then truncates offset line between those points
for j = 1 : length(y_offset)
    if y_offset(j) <= 0
        i=j+1;
    end
    if y_offset(j) >= stress(j)
        break
    end
end
x_offset = strain(i:j);
y_offset = y_offset(i:j);
plot(x_offset,y_offset, 'k')

%FAILURE POINT DATA
i = length(load);
fail_load = load(i);
disp_to_fail = displacement(i);
fail_stress = stress(i);
strain_to_fail = strain(i);

%ULTIMATE LOAD POINT DATA
[ultimate_load,i] = max(load);
disp_to_ult = displacement(i);
ultimate_stress = stress(i);
strain_to_ult = strain(i);
ultimate_index = i;

%YIELD POINT DATA
if j > ultimate_index
    j=ultimate_index;
end
yield_load = load(j);
disp_to_yield = displacement(j);
yield_stress = stress(j);
strain_to_yield = strain(j);
yield_index = j;

%Get postyield deformation/strain
postyield_disp = disp_to_fail - disp_to_yield;
postyield_strain = strain_to_fail - strain_to_yield;


%**************************************************************************
%Find pre and post yield energies and toughnesses
%Divide curves up into pre- and post-yield regions. 
strain1 = strain(1:yield_index);
stress1 = stress(1:yield_index);
load1 = load(1:yield_index);
displacement1 = displacement(1:yield_index);

%Calculate areas under curves
preyield_toughness = trapz(strain1,stress1) / 10^6;            % In MPa
total_toughness = trapz(strain,stress) / 10^6;
postyield_toughness = total_toughness - preyield_toughness;

preyield_work = trapz(displacement1,load1) / 10^3;             % In mJ
total_work = trapz(displacement,load) / 10^3;
postyield_work = total_work - preyield_work;


%***********************************************************************
%Plot final graphs of stress/strain
close
figure(3)

%Stress-strain plot
subplot(2,1,1)
plot(strain,stress)
axis xy
xlabel('Strain (microstrain)')
ylabel('Stress (MPa)')
hold on
%plot(linear_strain,linear_stress,'r')
plot(x_offset,y_offset, 'k')
plot(strain_to_yield, yield_stress, 'k+', strain_to_ult, ultimate_stress, 'k+', ...
     strain_to_fail, fail_stress, 'k+')
hold off

%Load-displacement plot
subplot(2,1,2)
plot(displacement,load)
axis xy
xlabel('Displacement (microns)')
ylabel('Force (N)')
hold on
plot(disp_to_yield, yield_load, 'k+', disp_to_ult, ultimate_load, 'k+', ...
     disp_to_fail, fail_load, 'k+')
hold off
 
yyy=menu('Would you like to reselect these points?','Yes','No');

end

%**************************** OUTPUT *********************************************

% Saves an image of figure 3 (summary of mechanical properties)
print ('-dpng', specimen_name) 

% Writes values for mechanical properties to analyze to a xls file with column headers. There
% will be an empty cell afer which outputs for a schematic
% representation of the f/d and stress/strain curves will appear.

headers = {'Specimen','I_ml (mm^4)','c_ant (µm)','Yield Force (N)','Ultimate Force (N)','Displacement to Yield (µm)','Postyield Displacement (µm)','Total Displacment (µm)','Stiffness (N/mm)','Work to Yield (mJ)','Postyield Work (mJ)','Total Work (mJ)','Yield Stress (MPa)','Ultimate Stress (MPa)','Strain to Yield (µ?)','Total Strain (µ?)','Modulus (GPa)','Resilience (MPa)','Toughness (MPa)',' ','Specimen','Yield Force (N)','Ultimate Force (N)','Failure Force (N)','Displacement to Yield (µm)','Ultimate Displacement (µm)','Total Displacment (µm)','Yield Stress (MPa)','Ultimate Stress (MPa)','Failure Stress (MPa)','Strain to Yield (µ?)','Ultimate Strain (µ?)','Total Strain (µ?)'};

resultsxls = [{specimen_name, num2str(I), num2str(c), num2str(yield_load), ...
        num2str(ultimate_load), num2str(disp_to_yield), num2str(postyield_disp), num2str(disp_to_fail), ...
        num2str(stiffness), num2str(preyield_work), num2str(postyield_work), ...
        num2str(total_work), num2str(yield_stress), num2str(ultimate_stress), ...
        num2str(strain_to_yield), num2str(strain_to_fail), num2str(modulus),  ...
        num2str(preyield_toughness), num2str(total_toughness), '', specimen_name, ...
        num2str(yield_load), num2str(ultimate_load), num2str(fail_load), ...
        num2str(disp_to_yield), num2str(disp_to_ult), num2str(disp_to_fail), ...
        num2str(yield_stress), num2str(ultimate_stress), num2str(fail_stress), ...
        num2str(strain_to_yield), num2str(strain_to_ult), num2str(strain_to_fail)}]; 

row=num2str(ppp+1);
rowcount=['A' row];

xlswrite('femur_mechanics.xls', resultsxls, 'Data',rowcount)
xlswrite('femur_mechanics.xls', headers, 'Data', 'A1')
warning off MATLAB:xlswrite:AddSheet 

ppp=ppp+1;
zzz=menu('Do you have more data to analyze?','Yes','No');
    
end

toc