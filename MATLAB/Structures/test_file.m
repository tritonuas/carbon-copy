E1 = 135e9; % Pa
E2 = 10e9; % Pa
G12 = 5e9; % Pa
Nu12 = 0.30;
thickness_per_ply = 0.005/2;
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
print_output = false;
SF = 100;

Nx = -4000;
Ny = 0;
Nxy = 0;
Mx = 0;
My = 0;
Mxy = 0;
mech_loading = [Nx,Ny,Nxy,Mx,My,Mxy]';
delta_T = 0;
thetas = [0,90];
rad_or_deg = "deg";
thicknesses = ones(length(thetas),1)*thickness_per_ply;

[stresses_bot, stresses_top,z_all,mid_strains_and_curvatures,...
     thermal_loading, ABD,MS, failed_plies, failed_side, failed_z,...
     fail_mode, fail_tcs] = structural_model(mat_props, thetas, ...
     rad_or_deg, thicknesses, mech_loading, delta_T, cte_vec,...
     mat_strengths_t,fail_crit,mat_strengths_c, SF, print_output)