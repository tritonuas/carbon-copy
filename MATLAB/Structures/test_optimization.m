%%Optimization of a composite layup using CLPT for analysis.
%%Author: Andrew Fletcher
%%Email: afletche@ucsd.edu

%% General
clear;
clc;

format shorte;
format compact;

%% Input Conditions

E1 = 2e7;
E2 = 1.5e6;
G12 = 1e6;
Nu12 = 0.29;
thickness_per_ply = 0.005;
cte1 = -5e-7;
cte2 = 1.5e-5;
sigma_1T = 3.1e5;
sigma_1C = -2e5;
sigma_2T = 9e3;
sigma_2C = -3e4;
sigma_12 = 1.5e4;

mat_props = [E1;E2;G12;Nu12];
cte_vec = [cte1; cte2; 0];
mat_strengths_t = [sigma_1T;sigma_2T;sigma_12];
mat_strengths_c = [sigma_1C;sigma_2C;sigma_12];

fail_crit = "max_stress";
print_output = false;
SF = 2;

Nx = 200;
Ny = 0;
Nxy = 0;
Mx = 0;
My = 0;
Mxy = 0;
delta_T = 0;

mech_loading = [Nx;Ny;Nxy;Mx;My;Mxy];

%% Design vars
thetas = symm([90]);
thicknesses = ones(1:length(thetas))*thickness_per_ply;
rad_or_deg = 'deg';
theta_h = 0.1;  %deg

% %% Objective
tol = 1e-12;

%% Optimization


thetas_hist = thetas;
obj_hist = 1;
grad_hist = ones(length(thetas)/2);
grad_norm_hist = 1;
hess_hist = eye(length(thetas)/2);

hess_init = eye(length(thetas)/2);
iter = 0;
while true
    iter = iter + 1;
    
    for i = 1:length(thetas)
        if thetas(i) > 90
            thetas(i) = thetas(i) - 180;
        elseif thetas(i) < -90
            thetas(i) = thetas(i) + 180;
        end
    end
    
    thetas_hist(iter,:) = thetas;  %keeping dv history
    
    % model / analysis valuation
    MS = CLPT(mat_props, thetas, rad_or_deg, thicknesses, mech_loading,...
                        delta_T, cte_vec, fail_crit, mat_strengths_t, ...
                        mat_strengths_c, SF, print_output);
    obj_hist(iter) = -min(MS(:));

    % preallocating for optimizer
    obj_grad = ones(1,length(thetas)/2);
    % calc derivatives using finite difference
    for i = 1:length(thetas)/2
        thetas(i) = thetas(i) + theta_h;
        thetas = symm(thetas(1:length(thetas)/2));
        MS_dist = CLPT(mat_props, thetas, rad_or_deg, thicknesses, ...
            mech_loading, delta_T, cte_vec, fail_crit, mat_strengths_t, ...
            mat_strengths_c, SF, print_output);
        d_MS = min(min(MS_dist)) - min(min(MS));
            
        thetas(i) = thetas(i) - theta_h;
        thetas = symm(thetas(1:length(thetas)/2));
        obj_grad(i) = -d_MS/theta_h
    end
    grad_hist(iter,:) = obj_grad;
    
    
    % calc hessian using secant method
    if iter > 1
%     if false
        d_obj_grad = obj_grad - grad_hist(iter-1);
        obj_hess = diag(d_obj_grad./(pk*alpha));
    else
        obj_hess = hess_init;
    end
    hess_hist(iter,:,:) = obj_hess;
    
    % break condition
    obj_grad_norm = norm(obj_grad, 1);
    grad_norm_hist(iter) = obj_grad_norm;
    if(obj_grad_norm < tol)
        break;
    end
    
    % search direction vector
    pk = -obj_hess\obj_grad;
    
    % line search
    alpha = 1;
    
    % update design variables
    thetas = thetas + pk*alpha
end

%% Output
figure(1);
plot(thetas_hist)
figure(2);
plot(grad_hist)
figure(3);
plot(grad_norm_hist)

%% Test
thetas = linspace(-90, 90, 5000);
for i = 1:length(thetas)
    MS = CLPT(mat_props, [thetas(i) thetas(i)], rad_or_deg, thicknesses, mech_loading,...
    delta_T, cte_vec, fail_crit, mat_strengths_t, mat_strengths_c, SF, ...
    print_output);
    obj(i) = -min(min(MS));
    
    thetas_dist = [thetas(i)+theta_h thetas(i)+theta_h];
    MS_dist = CLPT(mat_props, thetas_dist, rad_or_deg, thicknesses, ...
        mech_loading, delta_T, cte_vec, fail_crit, mat_strengths_t, ...
        mat_strengths_c, SF, print_output);
    d_MS = min(min(MS_dist)) - min(min(MS));

    obj_grad(i) = -d_MS/theta_h;
end

plot(obj)
plot(obj_grad)

%% CLPT Analysis Function Definition
function [MS, failed_plies, failed_side, failed_z, fail_mode, fail_tcs] ...
    = CLPT(mat_props, thetas, rad_or_deg, thicknesses, mech_loading,...
    delta_T, cte_vec, fail_crit, mat_strengths_t, mat_strengths_c, SF, ...
    print_output)
[stresses_bot, stresses_top, z_all] = ...
get_local_lamina_stresses_planar_ortho(mat_props, thetas, ...
rad_or_deg, thicknesses, mech_loading, delta_T, cte_vec);

[MS, failed_plies, failed_side, failed_z, fail_mode, fail_tcs] ...
= report_ply_margins(stresses_bot, stresses_top, z_all, ...
fail_crit, mat_strengths_t, mat_strengths_c, SF, print_output);
end

% %% Example CLPT Analysis
% [stresses_bot, stresses_top, z_all, ...
% mid_strains_and_curvatures, thermal_loading, ABD] = ...
% get_local_lamina_stresses_planar_ortho(mat_props, thetas, ...
% rad_or_deg, thicknesses, mech_loading, delta_T, cte_vec);
% 
% [MS, failed_plies, failed_side, failed_z, fail_mode, fail_tcs] ...
% = report_ply_margins(stresses_bot, stresses_top, z_all, ...
% fail_crit, mat_strengths_t, mat_strengths_c, SF, print_output);