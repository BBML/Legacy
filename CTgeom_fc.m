function [ ]  = CTgeom_fc()
%% Revision History

% Edited Feb 4 2020 by Rachel Kohler to add error explanation if folders
% are setup wrong.

% Edited Oct 3 2019 by Rachel Kohler to output c_ant (anterior extreme).
% Further edits done Oct 10 2019 to clean up code (deleting unused code).

% Edited May 2015 by Max Hammond to optimize code, record a diary,
% calculate TMD, read BMD equations from multiple sub-folders, reduce Excel
% write functions, removed xlswrite1 dependence, fill pores, apply
% greyscale threshold, remove naming dependence, alter parameter output
% order, and alter figure output. Outdated code is commented out. Note that
% once code is commented out it may not function as expected if re-inserted
% given subsequent changes.

% Edited Sept 2014 by Max Hammond to add code from Alycia Berman that
% excludes scale bars in the images at line 105. Changed outpout to xlsx
% format and added headers. Adjust profile output to go from 0 to 360.
% Added studynum so the loop doesn't have to be hard coded. Used xlswrite1
% which only opens Excel once per file to speed up program. This adds a
% dependancy to xlswrite1.m which must be in the directory to run the file

% Written by Joey Wallace, Sept 2011

%% Function Overview

% This program reads grayscale .bmp images from CT, applies a theshold,
% removes everything except the bone, fills in pores, and calculates
% relevant geometric properties. Each slice is calculated individually.
% Four Excel spreadsheets are output containing slice by slice or average
% geometric properties or profiles. A .png image is output for each bone
% within the respective folder showing the bone's profile and major/minor
% axes.

% If there an error is generated during analysis, the program will save the
% data already analyzed if any and display a warning with information on
% when the error occured along with the actual error message.

%% Proper Setup

% These files should be vertically aligned and oriented with anterior to the
% right.  Therefore, right limbs will be oriented medial up and left limbs
% will be oriented medial down.

% Each day of scanning should have its own folder in the parent directory.
% Each bone's ROI should be saved as a separate folder named with sample's
% ID# in the fodler corresponding to the day it was scanned. For example, if
% sample 50 was scanned on 7/3/15 and sample 11 was scanned on 11/14/14, in
% a folder named 'Standard Site' you would have a folder named '7-3-15' and a
% folder named '11-14-14'. Within '7-3-15' there would be a folder named
% '50' containing the slice .bmp files of the standard site and within
% '11-14-14' there would be a folder named '11' containing its .bmps.

% There is no longer a naming convention that you must follow. Each bone
% can have a different number of slices, but all bones must have the sample
% resolution, be the same bone, and be from the same side.

%% Initialization
diary off
close all % close all figures
clear all % clear all variables
format long % change format to long
warning('off','MATLAB:xlswrite:AddSheet'); % disable add sheet warning for Excel

%% Input parameters and check for common erros

% Select the folder containing all the standard site ROIs separated into
% folders
overall_folder=uigetdir;
cd(overall_folder)
overall_listing=dir(overall_folder);

% Create an array of all the sample folders in the dataset
folders = overall_listing(arrayfun(@(x) x.name(1), overall_listing) ~= '.');
folders = folders(arrayfun(@(x) x.isdir(1), folders) ~= 0);

% Create a loop to count how many total samples there are in all sub
% folders
total_count = 0;

for m=1:length(folders)
    sub_folder=[overall_folder '\' folders(m,1).name];
    subfolder_listing=dir(sub_folder);
    sub_folders = subfolder_listing(arrayfun(@(x) x.name(1), subfolder_listing) ~= '.');
    sub_folders = sub_folders(arrayfun(@(x) x.isdir(1), sub_folders) ~= 0);
    total_count = total_count + length(sub_folders); 
end

diary([overall_folder '\Diary.txt'])
diary on
res = input ('\n\n Voxel resolution in um: '); % input the isotropic voxel size from CTAn

% Input step used for extreme calculations (typically 0.5) and display a
% warning for an odd angular resolution
ang = input (' Angle step in degrees (factor of 360): ');
while mod(360,ang)
    fprintf(2,'\n Please enter a factor of 360 (e.g. 0.5). Note: this is NOT the angular step size used in the uCT scan. \n\n')
    ang = input ('Angle step in degrees (factor of 360): ');
end

% Input global greyscale threshold and display a warning for values outside
% of the specified range
threshold = input(' Input threshold value (0-255): ');
while threshold>255 || threshold<0
    fprintf(2,'\n Please enter a threshold between 0 and 255. \n\n')
    threshold = input('Input threshold value (0-255): ');
end

% Input whether left or right limbs will be analyzed, correct common
% mistakes, and display a warning for other values
side = input(' Enter "l" for a left limb and "r" for a right limb: ', 's');
side = strrep(side,'"','');
side = lower(side);
side = side(1:1);
side_log = ~strcmp(side,'l') + ~strcmp(side,'r');
while side_log~=1
    fprintf(2, ['\n You entered "' side '". \n\n'])
    side = input(' Enter "l" for a left limb and "r" for a right limb: ', 's');
    side = strrep(side,'"','');
    side = lower(side);
    side = side(1:1);
    side_log = ~strcmp(side,'l') + ~strcmp(side,'r');
end

% Input whether femora or tibiae will be analyzed, correct common mistakes,
% and display a warning for other values
bone = input(' Enter "f" for a femur and "t" for a tibia: ', 's');
bone = strrep(bone,'"','');
bone = lower(bone);
bone = bone(1:1);
bone_log = ~strcmp(bone,'f') + ~strcmp(bone,'t');
while bone_log~=1
    fprintf(2, ['\n You entered "' bone '". \n\n']) 
    bone = input(' Enter "f" for a femur and "t" for a tibia: ', 's');
    bone = strrep(bone,'"','');
    bone = lower(bone);
    bone = bone(1:1);
    bone_log = ~strcmp(bone,'f') + ~strcmp(bone,'t');
end

if length(folders)==0
    fprintf(2,'\n----------------------ERROR-----------------------\n');
    fprintf(2,'Your files are not properly organized.\nSee line 43 in the code and set up your folders according to the prescribed hierarchy.\n');
    return
end

eq_num = zeros(length(folders),1);
eq_denom = zeros(length(folders),1);

for m=1:length(folders)
    sub_folder=[overall_folder '\' folders(m,1).name];
    day = sub_folder(length(overall_folder)+2:end);
    fprintf(['\n Enter attenuation coefficients (0-0.11) in CTAn''s BMD equation  for ' day '. \n'])
    
    
    % Input HA calibration equation from CTAn and display a warning if
    % needed. CTAn uses an equation to calculate BMD for every pixel
    % according to BMD = (AC - AC_min)/AC_range
    eq_n = input(' Enter the minimum attenuation coefficient (numerator): ');
    while eq_n>.11 || eq_n<0
        fprintf(2,'\n Please enter a proper minimum attenuation coefficient (0-0.11). \n\n')
        eq_n = input(' Enter the minimum attenuation coefficient (numerator): ');
    end
    
    eq_num(m, 1) = eq_n;
    
    eq_d = input(' Enter the attenuation coefficient range (denominator): ');
    while eq_d>.11 || eq_d<0 || eq_d<eq_n
        fprintf(2,'\n Please enter a proper attenuation coefficient range (0-0.11, > minimum). \n\n')
        eq_d = input(' Enter the attenuation coefficient range (denominator): ');
    end
    
    eq_denom(m, 1) = eq_d;
    
end

diary off

tot = tic; % start the stop watch

%% Calculation

try % use a try/catch block to run the code and save anything that has already ran
    
    % Preallocate variabless
    A = 360/ang+1;
    prof_out_peri_cell = cell(total_count, A);
    prof_out_endo_cell = cell(total_count, A);
    geom_out_cell = cell(total_count, 19);
    sample_list = cell(1, total_count);
    ac_step = .11/255;
    offset = 0;
    
    % Create a loop to index all the sub_folders
    for m=1:length(folders)
        sub_folder=[overall_folder '\' folders(m,1).name];
        cd(sub_folder);
        
        subfolder_listing=dir(sub_folder);
        
        % Create an array of all the sample folders in the dataset
        sub_folders = subfolder_listing(arrayfun(@(x) x.name(1), subfolder_listing) ~= '.');
        sub_folders = sub_folders(arrayfun(@(x) x.isdir(1), sub_folders) ~= 0);
        
        % Create loop to analyze every sample's folder within the overall
        % folder
        if length(sub_folders)==0
            fprintf(2,'\n----------------------ERROR-----------------------\n');
            fprintf(2,'Your files are not properly organized.\nSee line 43 in the code and set up your folders according to the prescribed hierarchy.\n');
            return
        end
            
        for k=1:length(sub_folders)
            
            tic
            
            clearvars -except offset res tot total_count ang side bone eq_num eq_denom overall_folder sub_folder overall_listing subfolder_listing threshold  sub_folders folders A prof_out_peri_cell prof_out_endo_cell sample_list geom_out_cell m k ac_step
            
            % Store the .bmp filenames in the folder
            filename=[sub_folder '\' sub_folders(k,1).name];
            cd(filename);
            slices=dir([filename '\*.bmp']);
            
            % Store the folder name to be used as the ID during output
            folder = (sub_folders(k,1).name);
            
            % Pre-allocating arrays for data output later
            peri_out = zeros(length(slices),A);
            endo_out = zeros(length(slices),A);
            profiles = zeros((2*length(slices)+3),A);
            geom_out = zeros(length(slices)+2,19);
            tmd_gs = zeros(length(slices), 3);
            centroid_x = zeros(length(slices), 1);
            centroid_y = zeros(length(slices), 1);
            
            % Create loop to calculate parameters for each slice
            for j=1:length(slices)
                
                % Read in each slice as a grayscale image
                section = imread(slices(j,1).name);
                
                % Read in each slice as a BW image and allow the variable
                % to change size to accommodate different sized ROIs and/or
                % different number of slices
  
                slice = imbinarize(section,(threshold-1)/255); % use threshold-1 to take everything greater than or equal to the input threshold like CTAn
                
                % Remove all but the largest connected component (i.e.
                % remove scales or fibula if applicable)
                cc=bwconncomp(slice); % find the connected components in the image
                numPixels = cellfun(@numel,cc.PixelIdxList); % find the number of pixels in each component
                [~,idx] = max(numPixels); % find the index containing the most number of pixels
                for i=1:length(cc.PixelIdxList) % remove all other components
                    if i==idx
                        % do nothing
                    else
                        slice(cc.PixelIdxList{i})=0;
                    end
                end
                clear cc numPixels idx i
                
                stats = regionprops(slice,section,'MeanIntensity','PixelValues');
                tmd_gs(j,:) = [stats.MeanIntensity sum([stats.PixelValues]) length([stats.PixelValues])];
                tmd_ac = tmd_gs(1).*ac_step;
                tmd_HA = (tmd_ac - eq_num(m,1))/eq_denom(m,1);
                
                % Manually calculate the centroid from pixel locations:
                [index_y,index_x] = find (slice == 1); % this finds the x and y locations of each "on" pixel
                Qx = sum(index_y); % since the area of each dA is 1 pixel, we have x1 + x2 +...+ and this is the integral of y_dA
                Qy = sum(index_x); % since the area of each dA is 1 pixel, this is the integral of x_dA
                area = length(index_y); % since the area of each pixel is one, this is the total number of pixels or the area
                xbar = Qy/area; % xbar = integral of y_da/A
                ybar = Qx/area; % ybar = integral of x_da/A
                
                % PLOT THE IMAGE WITH THE CENTROID MARKED
                % Start getting line profiles at various degrees:
                inner_fiber = zeros(1, A); % creats a zero vector for the endocortical radii
                outer_fiber = zeros(1, A); % creates a zero vector for the periosteal radii
                thickness = zeros(1, A-1);  % creates a zero vector for cortical thicknesses
                index = 0; % intialize index
                
                image_size = size(slice);
                x_size = image_size(2);
                y_size = image_size(1);
                
                % IN QUADRANT 1 (45 to 134.99 degrees, top quadrant):
                for i = -45:ang:44.9
                    index = index + 1;
                    angle = i * pi / 180;
                    yi=[ybar,0];
                    xi=[xbar,xbar-(tan(angle)*yi(1))];
                    [cx,cy,c] = improfile(slice,xi,yi,10000);
                    % improfile chooses an arbirtray number of points to
                    % look at so I will choose alot to be accurate. cx and
                    % cy are the pixel locations along the line and c is
                    % the intensity at each point
                    cort_on = find(c == 1);
                    thick = length(cort_on);
                    on_1 = cort_on(1);
                    on_end = cort_on(thick);
                    rad_x = [cx(on_1),cx(on_end)];
                    rad_y = [cy(on_1),cy(on_end)];
                    points=[xbar,ybar;rad_x(1),rad_y(1);rad_x(2),rad_y(2)]; % centroid, endo and peri points along this line
                    radii = pdist(points); % radii(1) is endo, radii(2) is peri and radii(3) is c_thickness
                    inner_fiber(1, index) = radii(1);
                    outer_fiber(1, index) = radii(2);
                    thickness(1, index) = radii(3);
                end
                
                
                % IN QUADRANT 2 (135 to 224.99 degrees, P quadrant):
                for i = -45:ang:44.9
                    index = index + 1;
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
                    inner_fiber(1, index) = radii(1);
                    outer_fiber(1, index) = radii(2);
                    thickness(1, index) = radii(3);
                end
                
                % IN QUADRANT 3 (225 to 314.99 degrees, bottom quadrant):
                for i = -45:ang:44.9
                    index = index + 1;
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
                    inner_fiber(1, index) = radii(1);
                    outer_fiber(1, index) = radii(2);
                    thickness(1, index) = radii(3);
                end
                
                % IN QUADRANT 4 (315 to 404.99 or 49.99, A quadrant):
                for i = -45:ang:44.9
                    index = index + 1;
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
                    inner_fiber(1, index) = radii(1);
                    outer_fiber(1, index) = radii(2);
                    thickness(1, index) = radii(3);
                end
                
                % To plot this in polar,you need to append the inner and
                % outer vectors with the vaule at 360 degrees (0 deg) to
                % close
                index = index + 1;
                inner_fiber(1, index) = inner_fiber(1);
                outer_fiber(1, index) = outer_fiber(1);
                
                % Setup angles for polar plot
                angle_deg = 45:ang:405;
                angle_rad = angle_deg.*pi./180; %convert to radians

                % Convert the geometric data from angle and radius to x and
                % y coordinates
                outer_fiber_x = outer_fiber.*cos(angle_rad);
                outer_fiber_y = outer_fiber.*sin(angle_rad);
                inner_fiber_x = inner_fiber.*cos(angle_rad);
                inner_fiber_y = inner_fiber.*sin(angle_rad);
                
                % Before shifting the origin from the centroid, calculate
                % extreme fiber in each anatomic direction.  A and P are
                % not dependent on wheter this is a right or left bone, but
                % M and L are:
                anterior_extreme = abs(max(outer_fiber_x));
                posterior_extreme = abs(min(outer_fiber_x));
                
                if side == 'r'
                    medial_extreme = abs(max(outer_fiber_y));
                    lateral_extreme = abs(min(outer_fiber_y));
                else
                    medial_extreme = abs(min(outer_fiber_y));
                    lateral_extreme = abs(max(outer_fiber_y));
                end
                
                % Shift the coordinate system from (0,0) at centroid to the
                % (0,0) at LL corner.  For geometric properties, the outer
                % perimeter needs to go in the CW direction.  Currently, it
                % is CCW so it needs to be flipped.  This does both and
                % plots to verify
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
                
                % USE THESE OUTPUTS PRIOR TO INCORPORATING POLYGEOM TO GET
                % SOME OF THE GEOMETRIC PROPERTIES OF INTEREST
                
                % Convert all pixel values to um and plot
                outer_fiber_x = outer_fiber_x*res;
                outer_fiber_y = outer_fiber_y*res;
                inner_fiber_x = inner_fiber_x*res;
                inner_fiber_y = inner_fiber_y*res;
                x_perimeter = x_perimeter*res;
                y_perimeter = y_perimeter*res;
                
                subplot(1,2,2)
                plot(x_perimeter,y_perimeter)
                hold on
                axis equal
                axis tight
                xlabel('x-position in um')
                ylabel('y-position in um')
                title('Calculated Profiles')
                
                % Calculate the average cortical thickness
                cort_thickness = mean(thickness);
                cort_thickness = cort_thickness*res;
                
                % Calculate extreme fiber in each anatomic direction in um
                anterior_extreme = anterior_extreme * res;
                medial_extreme = medial_extreme * res;
                posterior_extreme = posterior_extreme * res;
                lateral_extreme = lateral_extreme * res;
                
                % Calculate AP and ML diameters and AP/ML
                AP_width = anterior_extreme + posterior_extreme;
                ML_width = medial_extreme + lateral_extreme;
                APtoMLratio = AP_width/ML_width;
                
                %*****ADD POLYGEOM TO GET TOTAL CROSS SECTIONAL AREA AND
                %PERISOTEAL PERIMETER*****
                clear x y
                x = outer_fiber_x;
                y = outer_fiber_y;
                
                % Check if inputs are same size
                if ~isequal( size(x), size(y) )
                    error( 'X and Y must be the same size');
                end
                
                % Number of vertices
                [ x, ~ ] = shiftdim( x ); 
                [ y, ~ ] = shiftdim( y ); 
                [ n, ~ ] = size( x ); 
                
                % Temporarily shift data to mean of vertices for improved
                % accuracy
                xm = mean(x);
                ym = mean(y);
                x = x - xm*ones(n,1);
                y = y - ym*ones(n,1);
                
                % Delta x and delta y
                dx = x ( [ 2:n 1 ] ) - x;
                dy = y ( [ 2:n 1 ] ) - y;
                
                % Output of Total Cross Sectional Area
                total_cs_area = sum( y.*dx - x.*dy )/2;
                periosteal_BS = sum( sqrt( dx.*dx +dy.*dy ) );
                
                
                %*****ADD PART OF POLYGEOM TO GET CORTICAL AND MARROW
                %AREAS*****
                clear x y
                x = fliplr(inner_fiber_x);
                y = fliplr(inner_fiber_y);
                
                % check if inputs are same size
                if ~isequal( size(x), size(y) )
                    error( 'X and Y must be the same size');
                end
                
                % Number of vertices
                [ x, ~ ] = shiftdim( x ); 
                [ y, ~ ] = shiftdim( y ); 
                [ n, ~ ] = size( x ); 
                
                % Temporarily shift data to mean of vertices for improved
                % accuracy
                xm = mean(x);
                ym = mean(y);
                x = x - xm*ones(n,1);
                y = y - ym*ones(n,1);
                
                % Delta x and delta y
                dx = x ( [ 2:n 1 ] ) - x;
                dy = y ( [ 2:n 1 ] ) - y;
                
                % Output of Marrow Area and Cortical Area
                marrow_area = sum( y.*dx - x.*dy )/2;
                cortical_area = total_cs_area - marrow_area;
                endocortical_BS = sum( sqrt( dx.*dx +dy.*dy ) );
                
                %*****NOW INCORPORATE POLYGEOM TO GET OTHER PROPS*****
                clear x y
                x = x_perimeter;
                y = y_perimeter;
                
                % Check if inputs are same size
                if ~isequal( size(x), size(y) )
                    error( 'X and Y must be the same size');
                end
                
                % Number of vertices
                [ x, ~ ] = shiftdim( x ); 
                [ y, ~ ] = shiftdim( y ); 
                [ n, ~ ] = size( x ); 
                
                % Temporarily shift data to mean of vertices for improved
                % accuracy
                xm = mean(x);
                ym = mean(y);
                x = x - xm*ones(n,1);
                y = y - ym*ones(n,1);
                
                % Delta x and delta y
                dx = x( [ 2:n 1 ] ) - x;
                dy = y( [ 2:n 1 ] ) - y;
                
                % Summations for CW boundary integrals
                cA = sum( y.*dx - x.*dy )/2; % cortical area
                Axc = sum( 6*x.*y.*dx -3*x.*x.*dy +3*y.*dx.*dx +dx.*dx.*dy )/12; % first moment about the y-axis (xc*cA)
                Ayc = sum( 3*y.*y.*dx -6*x.*y.*dy -3*x.*dy.*dy -dx.*dy.*dy )/12; % first moment about the x-axis (yc*cA)
                Ixx = sum( 2*y.*y.*y.*dx -6*x.*y.*y.*dy -6*x.*y.*dy.*dy ...
                    -2*x.*dy.*dy.*dy -2*y.*dx.*dy.*dy -dx.*dy.*dy.*dy )/12;% second moment about x axis
                Iyy = sum( 6*x.*x.*y.*dx -2*x.*x.*x.*dy +6*x.*y.*dx.*dx ...
                    +2*y.*dx.*dx.*dx +2*x.*dx.*dx.*dy +dx.*dx.*dx.*dy )/12;% second moment about y axis
                Ixy = sum( 6*x.*y.*y.*dx -6*x.*x.*y.*dy +3*y.*y.*dx.*dx ...
                    -3*x.*x.*dy.*dy +2*y.*dx.*dx.*dy -2*x.*dx.*dy.*dy )/24;% product of inertia about x-y axes
                %P = sum( sqrt( dx.*dx +dy.*dy ) ); % perimeter?
                
                % Check for CCW versus CW boundary
                if cA < 0
                    cA = -cA;
                    Axc = -Axc;
                    Ayc = -Ayc;
                    Ixx = -Ixx;
                    Iyy = -Iyy;
                    Ixy = -Ixy;
                end
                
                % Centroidal moments
                xc = Axc / cA; % centroidal location in x direction
                yc = Ayc / cA; % centroidal location in y direction
                Iuu = Ixx - cA*yc*yc; % centroidal MOI about x axis
                Ivv = Iyy - cA*xc*xc; % centroidal MOI anout y axis
                Iuv = Ixy - cA*xc*yc; % product of inertia

                % Replace mean of vertices
                x_cen = xc + xm;
                y_cen = yc + ym;

                % Principal moments and orientation
                I = [ Iuu  -Iuv ;
                    -Iuv   Ivv ];
                [ eig_vec, eig_val ] = eig(I);
                I1 = eig_val(1,1); % principal MOI about 1 axis
                I2 = eig_val(2,2); % principal MOI about 2 axie
                ang1 = atan2( eig_vec(2,1), eig_vec(1,1) ); % orientation of 1 axis
                ang2 = atan2( eig_vec(2,2), eig_vec(1,2) ); % orientation of 2 axis
                
                % Plot the centoid output from Polygeom on last graph
                plot(x_cen,y_cen,'b+')
                
                % Preallocate memory to store the centroid for each slice
                % in a column vector and save the x and y coordinates
                
                centroid_x (j,:) = x_cen;
                centroid_y (j,:) = y_cen;
                
                % Section modulus is resistance to bending.  Here it is
                % Z=I/c where I is the centroidal MOI aboout the axis of
                % bending (the x axis) divided by the extreme fiber on the
                % failure surface (the medial surface is in tension) for a
                % tibia. A femur is tested about the mediolateral axis with
                % the anterior surface in tension
                
                if bone == 't'
                    section_mod = Iuu/medial_extreme;
                else
                    section_mod = Ivv/anterior_extreme;
                end
                
                %********************** SLICE OUTPUT***********************
                
                % For output purposes, we need the inner and outer fibers
                % converted to distance from pixels
                inner_fiber = inner_fiber * res;
                outer_fiber = outer_fiber * res;
                
                % Profiles can be dumbed into individual matirices which
                % will be added to duting this loop and then combined into
                % a single matrix after loop has ended for data output
                peri_out (j,:) = outer_fiber;
                endo_out (j,:) = inner_fiber;
                
                % Convert geometric props to proper units
                total_cs_area = total_cs_area * 1e-6;
                marrow_area = marrow_area * 1e-6;
                cortical_area = cortical_area * 1e-6;
                cort_thickness = cort_thickness * 1e-3;
                AP_width = AP_width * 1e-3;
                ML_width = ML_width * 1e-3;
                periosteal_BS = periosteal_BS * 1e-3;
                endocortical_BS = endocortical_BS * 1e-3;
                Iap = Iuu * 1e-12;
                Iml = Ivv * 1e-12;
                Imin = I1 * 1e-12;
                Imax = I2 * 1e-12;
                ang_min=ang1*180/pi;
                ang_max=ang2*180/pi;
                medial_extreme = medial_extreme * 1e-3;
                anterior_extreme = anterior_extreme * 1e-3;
                section_mod = section_mod * 1e-9;
                
                
                % Shift max and min angles into the 1st and 2nd Cartesian
                % quadrants
                if ang_min>=180
                    ang_min=ang_min-180;
                elseif ang_min<0
                    ang_min=ang_min+180;
                else
                    % do nothing
                end
                
                if ang_max>=180
                    ang_max=ang_max-180;
                elseif ang_max<0
                    ang_max=ang_max+180;
                else
                    % do nothing
                end
                
                % Store the output from each slice as a row in geom_out
                geometry = [total_cs_area marrow_area cortical_area cort_thickness ...
                    periosteal_BS endocortical_BS Imax ang_max Imin ang_min AP_width ML_width APtoMLratio  Iap ...
                    Iml section_mod medial_extreme anterior_extreme tmd_HA];
                geom_out(j,:) = geometry;
                
            end
            
            %********************** SAMPLE OUTPUT**************************
            
            % The matrice for angle, outer fiber, and inner fiber are
            % combined
            profiles (1,:) = 45:ang:405;
            profiles (2:length(slices)+1,:) = peri_out;
            profiles (length(slices)+2:2*length(slices)+1,:) = endo_out;
            profiles (2*length(slices)+2,:) = mean(profiles(2:length(slices)+1,:));
            profiles (2*length(slices)+3,:) = mean(profiles(length(slices)+2:2*length(slices)+1,:));
            
            % Overlay the average profile of all slices with the centroid
            % marked
            subplot(1, 2, 2)
            avg_peri = profiles (2*length(slices)+2,:);
            avg_endo = profiles (2*length(slices)+3,:);
            avg_peri_x = avg_peri.*cos(angle_rad);
            avg_peri_y = avg_peri.*sin(angle_rad);
            avg_endo_x = avg_endo.*cos(angle_rad);
            avg_endo_y = avg_endo.*sin(angle_rad);
            x_peri_min = abs(min(avg_peri_x));
            y_peri_min = abs(min(avg_peri_y));
            avg_peri_x = avg_peri_x+x_peri_min;
            avg_peri_x = fliplr(avg_peri_x);
            avg_peri_y = avg_peri_y+y_peri_min;
            avg_peri_y = fliplr(avg_peri_y);
            avg_endo_x = avg_endo_x+x_peri_min;
            avg_endo_y = avg_endo_y+y_peri_min;
            x_perimeter_avg = [avg_peri_x avg_endo_x];
            y_perimeter_avg = [avg_peri_y avg_endo_y];
            plot(x_perimeter_avg,y_perimeter_avg, 'r')
            hold on
            x_cen_avg = mean(centroid_x);
            y_cen_avg = mean(centroid_y);
            plot(x_cen_avg,y_cen_avg,'r+')
            axis equal
            axis tight
            xlabel('x-position in um')
            ylabel('y-position in um')
            title('Calculated Profiles')
            hold off
            
            % Prevent the polar subplot from shrinking (standard Matlab
            % error)
            h1 = subplot(1,2,1);  % save the handle of the subplot
            ax1=get(h1,'position'); % save the position as ax
            set(h1,'position',ax1);   % manually setting this holds the position
            
            % Create a polar plot of the profile with the Imax and Imin
            % axes
            subplot(1, 2, 1)
            polarplot(360,1000); % set  the axes for the polar plot
            hold on
            polarplot(angle_rad,avg_endo, '-b');
            polarplot(angle_rad,avg_peri, '-b');
            title('Polar Plot of Avg Profile')
            line_r = -max(avg_peri)+50:max(avg_peri)+50;
            ang_max_rad = mean(geom_out(1:length(slices),8)).*pi./180;
            ang_min_rad = mean(geom_out(1:length(slices),10)).*pi./180;
            ang_max_rad = ones(1,length(line_r)).*ang_max_rad;
            ang_min_rad = ones(1,length(line_r)).*ang_min_rad;
            pol_max = polarplot(ang_max_rad, line_r, '-r');
            pol_min = polarplot(ang_min_rad, line_r, '--r');
            legend([pol_max,pol_min],'\theta max','\theta min','Location','southoutside','Orientation','horizontal')
            hold off
            
            print ('-dpng', folder) % save the figure in the sample folder
            
            % The matrice for angle, outer fiber and inner fiber are
            % shifted to begin at 0 degrees
            Ao = 315/ang+1;
            ang_shift = 0:ang:360;
            prof_shift1 = profiles(1:2*length(slices)+3,Ao:end);
            prof_shift2 = profiles(1:2*length(slices)+3,2:Ao);
            prof_shift = horzcat(prof_shift1, prof_shift2);
            prof_shift(1, :) = ang_shift;
            prof_cell=num2cell(prof_shift);
 
            % Create a cell array for prof_avg
            prof_mean_peri = prof_shift(2*length(slices)+2,:);
            prof_cell_peri = num2cell(prof_mean_peri);
            prof_out_peri_cell(offset+k, :) = prof_cell_peri;
            prof_mean_endo = prof_shift(2*length(slices)+3,:);
            prof_cell_endo = num2cell(prof_mean_endo);
            prof_out_endo_cell(offset+k, :) = prof_cell_endo;
            
            % Calculate average TMD using total counts from each slice
            % instead of a simple average as before because the pixel count
            % will change between slices
            tot_gs = sum(tmd_gs(:,2), 1); % total greyscale intensity for all slices
            tot_px = sum(tmd_gs(:,3), 1); % total pixels for all slices
            tmd_gs_avg = tot_gs/tot_px; % average greyscale intensity for all slices
            tmd_ac_avg = tmd_gs_avg.*ac_step; % convert to attenuation coefficient (AC)
            tmd_HA_avg = (tmd_ac_avg - eq_num(m,1))/eq_denom(m,1); % convert AC to BMD using equation
            
            % Create a cell array for geom_avg
            geom_mean = mean(geom_out(1:length(slices),:));
            mean_cell=num2cell(geom_mean);
            geom_out_cell(offset+k, 1:18) = mean_cell(1: 18);
            geom_out_cell(offset+k, 19) = num2cell(tmd_HA_avg);
            sample_list{offset+k} = folder;
            
            timer = toc;
            out_msg = ['\n Sample ' num2str(folder) ' took ' num2str(timer) ' seconds.'];
            diary([overall_folder '\Diary.txt'])
            diary on
            fprintf(out_msg)
            diary off
        end
        
        offset = offset + length(sub_folders);
        
    end
    
catch ME
    
    if exist('k', 'var') && exist('folder', 'var') && exist('j', 'var')
        
        if k+offset>1 % only attempt to save data if this is not the first sample analyzed
            
            % Trim the variables to avoid dimension errors
            sample_list(:,all(cellfun(@isempty,sample_list),1)) = [];
            %geom_out_cell(all(cellfun(@isempty,geom_out_cell),2),:) = [];
            partial = length(sample_list);
            geom_out_cell(partial+1:end,:) = [];
            prof_out_peri_cell(partial+1:end,:) = [];
            prof_out_endo_cell(partial+1:end,:) = [];
            
            
            % Save the analysis that has already run
            peri_cell = ['Periosteal'; sample_list'];
            endo_cell = ['Endocortical'; sample_list'];
            col_prof_cell = [peri_cell; ' '; endo_cell];
            blank = cell(1, 360/ang +1);
            prof_out = [num2cell(0:ang:360); prof_out_peri_cell; blank; num2cell(0:ang:360); prof_out_endo_cell];
            prof_mean_out = horzcat(col_prof_cell, prof_out);
            xlswrite([overall_folder '\CTprof_avg_toerror.xlsx'], prof_mean_out, 'Raw Data', 'A1')
            
            geom_out_cell = horzcat(sample_list', geom_out_cell); % label the rows
            headers = {'Total CSA (mm^2)', 'Marrow Area (mm^2)', 'Cortical Area(mm^2)', 'Cortical Thickness (mm)', 'Periosteal BS (mm)', 'Endocortical BS (mm)',  'Imax (mm^4)', 'Theta max (deg)', 'Imin (mm^4)','Theta min (deg)', 'AP Width (mm)', 'ML Width (mm)', 'AP/ML',  'Iap (mm^4)', 'Iml (mm^4)',   'Section Modulus (mm^3)', 'Medial Extreme (mm)', 'Anterior Extreme (mm)','TMD (g/cm^3 HA)'};
            headers = horzcat(' ', headers);
            geom_mean_out = [headers; geom_out_cell]; % add column titles
            xlswrite([overall_folder '\CTgeom_avg_toerror.xlsx'], geom_mean_out, 'Raw Data', 'A1')
            
        end
        
        if isempty(j)==1
            
            diary([overall_folder '\Diary.txt'])
            diary on
            % Display the folder name where the error occured
            msg = ['\n Error occured at the beginning of loop ' num2str(k) '. No data could be read from folder ' num2str(folder) '.\n\n'];
            fprintf(2,msg) 
            
        else
            % Display the folder name where the error occured
            msg = ['\n Error occured during loop ' num2str(k) ' at folder ' num2str(folder) '  in slice ' num2str(j) '. \n\n'];
            fprintf(2,msg)
            
        end
        
    end
    
    tot_timer = toc(tot);
    tot_msg = ['  mn Total time for all samples to error was ' num2str(tot_timer) ' seconds.'];
    avg_msg = [' Average time per sample was ' num2str(tot_timer/total_count) ' seconds. \n\n'];
    fprintf(tot_msg)
    fprintf(avg_msg)
    
    rethrow(ME) % rethrow error for troubleshooting
    
end

%% Mean property and profile output for each bone

% Create cell array for output containing all column and row headers along
% with the profiles and save the data to an xlsx spreadsheet.
peri_cell = ['Periosteal'; sample_list'];
endo_cell = ['Endocortical'; sample_list'];
col_prof_cell = [peri_cell; ' '; endo_cell];
blank = cell(1, 360/ang +1);
prof_out = [prof_cell(1, :); prof_out_peri_cell; blank; prof_cell(1, :); prof_out_endo_cell];
prof_mean_out = horzcat(col_prof_cell, prof_out);
xlswrite([overall_folder '\CTprof_avg.xlsx'], prof_mean_out, 'Raw Data', 'A1')

% Create cell array for output containing all column and row headers along
% with the geometric properties and save the data to an xlsx spreadsheet.
geom_out_cell = horzcat(sample_list', geom_out_cell); % label the rows
headers = {'Total CSA (mm^2)', 'Marrow Area (mm^2)', 'Cortical Area(mm^2)', 'Cortical Thickness (mm)', 'Periosteal BS (mm)', 'Endocortical BS (mm)',  'Imax (mm^4)', 'Theta max (deg)', 'Imin (mm^4)','Theta min (deg)', 'AP Width (mm)', 'ML Width (mm)', 'AP/ML',  'Iap (mm^4)', 'Iml (mm^4)',   'Section Modulus (mm^3)', 'Medial Extreme (mm)', 'Anterior Extreme (mm)','TMD (g/cm^3 HA)'};
headers = horzcat(' ', headers);
geom_mean_out = [headers; geom_out_cell]; % add column titles
xlswrite([overall_folder '\CTgeom_avg.xlsx'], geom_mean_out, 'Raw Data', 'A1')

tot_timer = toc(tot);
tot_msg = ['\n\n Total time for all samples was ' num2str(tot_timer) ' seconds.'];
avg_msg = [' Average time per sample was ' num2str(tot_timer/total_count) ' seconds. \n\n'];
diary([overall_folder '\Diary.txt'])
diary on
fprintf(tot_msg)
fprintf(avg_msg)
diary off

cd(overall_folder);

end

