function [ ]  = bend_beams_bose()

tic
% Written by Joey Wallace,July 2012 to work with test resources system.
% On tibial testing data (originally r4tib.m)

% Modified in August 2013 to work with femurs tested in 3 or 4 point bending

% Modified in Sept 2014 to run for tests with anterior surface in tension

% Edited by Max Hammond Sept. 2014 Changed the output from a csv 
% file to an xls spreadsheet that included a title row. Used while loop to
% semi-batch process. Hard coded in initial values like bendtype, slice
% number, and voxel size because each will be held constant within a study.
% Added the option to smooth or not during Testing Configuration. Smoothed
% using a moving average with a span of 10. Added a menu in case points
% need to be reselected.

%Final modification by JW on 09/27/2019 to correct small errors 

% This program reads raw mechanical data from the Bose system
% from a a file named "specimen name.csv". 
%  
% The program adjusts for system compliance and then uses beam bending
% theory to convert force-displacement data to theoretical stress-strain
% values.  Mechanical properties are calculated in both domains and output
% to a file "specimen name_mechanics.csv".  It also outputs a figure
% showing the load-displacement and stress-strain curves with significant
% points marked.

%close all figure windows and clears all variables
close all
clear all

warning('off','MATLAB:xlswrite:AddSheet');

%*****************\TESTING CONFIGURATION/**********************************
%                                                                         *
%       Adjust these values to match the system setup                     *
%                                                                         *
L = 18.68;           %span between bottom load points (mm)                 *
a = 3.00;           %distance between outer and inner points (if 4pt; mm) *
bendtype = '3';     %enter '4' for 4pt and '3' for 3pt bending            *
compliance = 0;     %system compliance (microns/N)                        *
smoothing = 0;      %enter 1 to smooth using moving average (span=10)     *
%**************************************************************************

%Check common errors in testing configuration
if bendtype ~= '3' && bendtype ~= '4'
        error('Please enter 3 or 4 for bendtype as a string in the Testing Configuration')
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
ID = input('Input bone ID: ','s');

w = input('Input beam width (mm): ');
h = input('Input beam height (mm): ');

I = (1/12)*w*h^3; %MOI for rectangular beam in mm^4
c = h/2 * 10^3; %extreeme fiber in um

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

%Truncate the data at the start point
position = position(j:i);
load = load(j:i);

%Adjust load starting position equal to 0 Newtons.
load = load - load(1);

%Adjust position for compliance and set starting position equal to 0.
displacement = position - load*compliance;
displacement = displacement - displacement(1);


%Convert the corrected load/displacement data to stress/strain
if bendtype == '3'
    stress = (load*L*c) / (4*I) * 10^-3;             %MPa
    strain = (12*c*displacement) / (L^2);             %microstrain
end

if bendtype == '4'
   stress = (load*a*c) / (2*I) * 10^-3;             %MPa
   strain = (6*c*displacement) / (a*(3*L - 4*a));   %microstrain
end

    
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

%Plot the adjusted stress-strain curve and pick points to define modulus
figure(2)
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

%Create stress and strain vectors spanning region for modulus determination
i=1;
while x(1) > strain(i)
    i = i+1;
end
linear_strain(1) = strain(i);
linear_stress(1) = stress(i);

j=2;
while x(2) > strain(i)
    linear_strain(j) = strain(i);
    linear_stress(j) = stress(i);
    i = i+1;
    j = j+1;
end

plot(linear_strain,linear_stress,'r')

%Determine modulus by linear regression of selected points
coeff = polyfit(linear_strain,linear_stress,1);
slope = coeff(1);
modulus = slope * 10^3;                                 % GPa

if bendtype == '3'
    stiffness = modulus*48*I / (L^3) * 10^3;   % N/mm
end

if bendtype == '4'
   stiffness = modulus*12*I / (a^2 * (3*L -4*a)) * 10^3;   % N/mm
end

%Plot linear regression on top of selected region
hold on
linear_stress = slope*linear_strain + coeff(2);
plot(linear_strain,linear_stress,'g')

% Create line with a .2% offset (2000 microstrain)
y_int = coeff(2)-slope*2000;        %y intercept
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
plot(linear_strain,linear_stress,'r')
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
print ('-dpng', ID) 

% Writes values for mechanical properties to analyze to a xls file with column headers. There
% will be an empty cell afer which outputs for a schematic
% representation of the f/d and stress/strain curves will appear.

headers = {'Sample ID','I (mm^4)','c (µm)','Yield Force (N)','Ultimate Force (N)','Displacement to Yield (µm)','Postyield Displacement (µm)','Total Displacment (µm)','Stiffness (N/mm)','Work to Yield (mJ)','Postyield Work (mJ)','Total Work (mJ)','Yield Stress (MPa)','Ultimate Stress (MPa)','Strain to Yield (µ?)','Total Strain (µ?)','Modulus (GPa)','Resilience (MPa)','Toughness (MPa)',' ','Specimen','Yield Force (N)','Ultimate Force (N)','Failure Force (N)','Displacement to Yield (µm)','Ultimate Displacement (µm)','Total Displacment (µm)','Yield Stress (MPa)','Ultimate Stress (MPa)','Failure Stress (MPa)','Strain to Yield (µ?)','Ultimate Strain (µ?)','Total Strain (µ?)'};

resultsxls = [{ID, num2str(I), num2str(c), num2str(yield_load), ...
        num2str(ultimate_load), num2str(disp_to_yield), num2str(postyield_disp), num2str(disp_to_fail), ...
        num2str(stiffness), num2str(preyield_work), num2str(postyield_work), ...
        num2str(total_work), num2str(yield_stress), num2str(ultimate_stress), ...
        num2str(strain_to_yield), num2str(strain_to_fail), num2str(modulus),  ...
        num2str(preyield_toughness), num2str(total_toughness), '', ID, ...
        num2str(yield_load), num2str(ultimate_load), num2str(fail_load), ...
        num2str(disp_to_yield), num2str(disp_to_ult), num2str(disp_to_fail), ...
        num2str(yield_stress), num2str(ultimate_stress), num2str(fail_stress), ...
        num2str(strain_to_yield), num2str(strain_to_ult), num2str(strain_to_fail)}]; 

row=num2str(ppp+1);
rowcount=['A' row];

xlswrite('beam_mechanics.xls', resultsxls, 'Data',rowcount)
xlswrite('beam_mechanics.xls', headers, 'Data', 'A1')
warning off MATLAB:xlswrite:AddSheet 

ppp=ppp+1;
zzz=menu('Do you have more data to analyze?','Yes','No');
    
end

toc