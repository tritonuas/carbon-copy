%%This script is for completing HW5 for SE 142
%%I have this here as a test script. Please do not use inappropriately
%%Author: Andrew Fletcher
%%Email: afletche@ucsd.edu

%% General
close all;
clear;
clc;

format compact;
format shorte;

Name = "Andrew Fletcher";
PID = "A14509074";

inFile = "HW5_P1_InputFile.xlsx";
outFile = "HW5_P1_OutputResults.xlsx";
inSheet = readcell(inFile);

start_time = cputime;
%% Problem 1

%Material Properties
rowStart = "7";
rowEnd = "7";
colStart = "B";
colEnd = "M";
range = strcat(colStart,rowStart,":",colEnd,rowEnd);
values = readmatrix(inFile,'Range',range);
variables = ["E1", ...
       "E2", ...
       "G12", ...
       "Nu12", ...
       "thickness_per_ply", ...
       "cte1", ...
       "cte2", ...
       "sigma_1T", ...
       "sigma_1C",...
       "sigma_2T",...
       "sigma_2C",...
       "sigma_12"];      
for i = 1:length(values)
    eval(sprintf('%s = %d',variables(i),values(i)));
end

E1 = 2e7;
E2 = 2e7;
G12 = values(:,3)';
Nu12 = values(:,4)';
thickness_per_ply = values(:,5)';
cte1 = values(:,6)';
cte2 = values(:,7)';
sigma_1T = values(:,8)';
sigma_1C = values(:,9)';
sigma_2T = values(:,10)';
sigma_2C = values(:,11)';
sigma_12 = values(:,12)';

mat_props = [E1;E2;G12;Nu12];
cte_vec = [cte1; cte2; 0];
mat_strengths_t = [sigma_1T;sigma_2T;sigma_12];
mat_strengths_c = [sigma_1C;sigma_2C;sigma_12];

%Loading
rowStart = "12";
rowEnd = "12";
colStart = "A";
colEnd = "G";
range = strcat(colStart,rowStart,":",colEnd,rowEnd);
values = readmatrix(inFile,'Range',range);
variables = ["Nx", ...
       "Ny", ...
       "Nxy", ...
       "Mx", ...
       "My", ...
       "Mxy", ...
       "delta_T"];      
for i = 1:length(values)
    eval(sprintf('%s = %d',variables(i),values(i)));
end

mech_loading = [Nx; Ny; Nxy; Mx; My; Mxy];

%Laminate Layup
rowStart = "16";
rowEnd = "25";
colStart = "A";
colEnd = "D";
range = strcat(colStart,rowStart,":",colEnd,rowEnd);
lam_mats = readmatrix(inFile,'Range',range);

ply_num = lam_mats(:,1);
thetas_backwards = lam_mats(:,2);
mat_ID = lam_mats(:,3);
thicknesses = lam_mats(:,4);

thetas = zeros(length(thetas_backwards),1);
for i = 1:length(thetas_backwards)
    thetas(ply_num(i)) = thetas_backwards(i);
end
rad_or_deg = "deg";

%% Part C
[stresses_bot, stresses_top, z_all, ...
    mid_strains_and_curvatures, thermal_loading, ABD] = ...
    get_local_lamina_stresses_planar_ortho(mat_props, thetas, ...
    rad_or_deg, thicknesses, mech_loading, delta_T, cte_vec);

%Output A
rowStart = "19";
rowEnd = "21";
colStart = "A";
colEnd = "F";
range = strcat(colStart,rowStart,":",colEnd,rowEnd);
values = ABD(1:3, 1:6);
writematrix(values, outFile,'Range',range);

rowStart = "19";
rowEnd = "21";
colStart = "G";
colEnd = "I";
range = strcat(colStart,rowStart,":",colEnd,rowEnd);
values = ABD(4:6, 4:6);
writematrix(values, outFile,'Range',range);

%Output B
rowStart = "24";
rowEnd = "26";
colStart = "B";
colEnd = "B";
range = strcat(colStart,rowStart,":",colEnd,rowEnd);
values = thermal_loading(1:3);
writematrix(values, outFile,'Range',range);

rowStart = "24";
rowEnd = "26";
colStart = "F";
colEnd = "F";
range = strcat(colStart,rowStart,":",colEnd,rowEnd);
values = thermal_loading(4:6);
writematrix(values, outFile, 'Range',range);

%Output C.1
rowStart = "32";
rowEnd = "34";
colStart = "B";
colEnd = "B";
range = strcat(colStart,rowStart,":",colEnd,rowEnd);
values = mid_strains_and_curvatures(1:3);
writematrix(values, outFile,'Range',range);

rowStart = "32";
rowEnd = "34";
colStart = "F";
colEnd = "F";
range = strcat(colStart,rowStart,":",colEnd,rowEnd);
values = mid_strains_and_curvatures(4:6);
writematrix(values, outFile,'Range',range);

%Output C.2
rowStart = "38";
rowEnd = "40";
colStart = "B";
colEnd = "B";
range = strcat(colStart,rowStart,":",colEnd,rowEnd);
values = stresses_top(:,9);
writematrix(values, outFile, 'Range',range);


%% Failure Loading
fail_crit = "max_stress";
print_output = false;
SF = 1;

Nx = 0;
Ny = 0.5*Nx;
failed_plies = 0;
while(sum(failed_plies)==0)
    Nx = Nx + 1;
    Ny = 0.5*Nx;
    mech_loading = [Nx;Ny;Nxy;Mx;My;Mxy];
    
    [stresses_bot, stresses_top, z_all, ...
    mid_strains_and_curvatures, thermal_loading, ABD] = ...
    get_local_lamina_stresses_planar_ortho(mat_props, thetas, ...
    rad_or_deg, thicknesses, mech_loading, delta_T, cte_vec);

    [MS, failed_plies, failed_side, failed_z, fail_mode, fail_tcs] ...
    = report_ply_margins(stresses_bot, stresses_top, z_all, ...
    fail_crit, mat_strengths_t, mat_strengths_c, SF, print_output);
end

%Output D
rowStart = "46";
rowEnd = "47";
colStart = "B";
colEnd = "B";
range = strcat(colStart,rowStart,":",colEnd,rowEnd);
values = [Nx;Ny];
writematrix(values, outFile,'Range',range);

rowStart = "46";
rowEnd = "46";
colStart = "E";
colEnd = "E";
range = strcat(colStart,rowStart,":",colEnd,rowEnd);
values = failed_plies;
writematrix(values, outFile, 'Range',range);

rowStart = "47";
rowEnd = "47";
colStart = "E";
colEnd = "E";
range = strcat(colStart,rowStart,":",colEnd,rowEnd);
values = failed_side;
writematrix(values, outFile, 'Range',range);

rowStart = "48";
rowEnd = "48";
colStart = "E";
colEnd = "E";
range = strcat(colStart,rowStart,":",colEnd,rowEnd);
values = failed_z;
writematrix(values, outFile, 'Range',range);

rowStart = "46";
rowEnd = "47";
colStart = "H";
colEnd = "H";
range = strcat(colStart,rowStart,":",colEnd,rowEnd);
values = [fail_mode; fail_tcs];
writematrix(values, outFile, 'Range',range);