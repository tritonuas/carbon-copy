%% EXAMPLE CLPT ANALYSIS
% CLPT(classical laminate plate theory)

% calculate direction and number of plies 

% defining material properties - carbon fiber
E1 = 135e9; % Pa
E2 = 10e9; % Pa
G12 = 5e9; % Pa
Nu12 = 0.30;
thickness_per_ply = 0.005;
cte1 = -5e-7;
cte2 = 1.5e-5;
sigma_1T = 1500e6;
sigma_1C = 1200e6;
sigma_2T = 50e6;
sigma_2C = 250e6;
sigma_12 = 70e6;

mat_props = [E1;E2;G12;Nu12];
cte_vec = [cte1; cte2; 0];
mat_strengths_t = [sigma_1T;sigma_2T;sigma_12];
mat_strengths_c = [sigma_1C;sigma_2C;sigma_12];

fail_crit = "max_stress";
print_output = true;
SF = 10; % safety factor 

% defining mechanical loads 
Nx = 0; 
Ny = 0;
Nxy = 30000;
Mx = 0;
My = 0;
Mxy = 0;
delta_T = 0;
t_airfoil = .1;
thetas = [45 45 45 45 45];
rad_or_deg = "deg";
thickness_per_ply = 0.0003;
thicknesses = ones(length(thetas),1)*thickness_per_ply;

mech_loading = [Nx;Ny;Nxy;Mx;My;Mxy];

[stresses_bot, stresses_top, z_all, ...
mid_strains_and_curvatures, thermal_loading, ABD] = ...
get_local_lamina_stresses_planar_ortho(mat_props, thetas, ...
rad_or_deg, thicknesses, mech_loading, delta_T, cte_vec);

[MS, failed_plies, failed_side, failed_z, fail_mode, fail_tcs] ...
= report_ply_margins(stresses_bot, stresses_top, z_all, ...
fail_crit, mat_strengths_t, mat_strengths_c, SF, print_output)