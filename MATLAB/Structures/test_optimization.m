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

Nx = 100;
Ny = 0;
Nxy = 0;
Mx = 0;
My = 0.5;
Mxy = 0;
delta_T = 0;

mech_loading = [Nx;Ny;Nxy;Mx;My;Mxy];

%% Design vars
thetas = symm([45 45]);
thicknesses = ones(1:length(thetas))*thickness_per_ply;
rad_or_deg = 'deg';
theta_h = 0.1;  %deg

[stresses_bot, stresses_top, z_all, ...
    mid_strains_and_curvatures, thermal_loading, ABD] = ...
    get_local_lamina_stresses_planar_ortho(mat_props, thetas, ...
    rad_or_deg, thicknesses, mech_loading, delta_T, cte_vec);

% %% Objective
tol = 1e-6;

%% Optimization

thetas_hist = thetas;
obj_hist = 1;
grad_hist = ones(length(thetas)/2);
grad_norm_hist = 1;
hess_hist = eye(length(thetas)/2);

hess_init = eye(length(thetas)/2);
iter = 0;
start_time = cputime;
while true
    iter = iter + 1;
    
%     for i = 1:length(thetas)
%         if thetas(i) > 90
%             thetas(i) = thetas(i) - 180;
%         elseif thetas(i) < -90
%             thetas(i) = thetas(i) + 180;
%         end
%     end
    
    thetas_hist(iter, :) = thetas;  %keeping dv history
    
    % model / analysis evaluation
    [stresses_bot, stresses_top] = ...
        get_local_lamina_stresses_planar_ortho(mat_props, thetas, ...
        rad_or_deg, thicknesses, mech_loading, delta_T, cte_vec);
    obj_hist(iter) = norm(stresses_bot./(mat_props(1:3)*ones(1, length(thetas)))...
        + stresses_top./(mat_props(1:3)*ones(1, length(thetas))));

    % preallocating for optimizer
    obj_grad = ones(1,length(thetas)/2);
    % calc derivatives using finite difference
    for i = 1:length(thetas)/2
        thetas(i) = thetas(i) + theta_h;
        thetas = symm(thetas(1:length(thetas)/2));
        [stresses_bot, stresses_top] = ...
            get_local_lamina_stresses_planar_ortho(mat_props, thetas, ...
            rad_or_deg, thicknesses, mech_loading, delta_T, cte_vec);
        obj_dist = norm(stresses_bot./(mat_props(1:3)*ones(1, length(thetas)))...
            + stresses_top./(mat_props(1:3)*ones(1, length(thetas))));
            d_obj = obj_dist - obj_hist(iter);
            
        thetas(i) = thetas(i) - theta_h;
        thetas = symm(thetas(1:length(thetas)/2));
        obj_grad(i) = d_obj/theta_h;
    end
    obj_grad = obj_grad';
    grad_hist(iter,:) = obj_grad;
    
    
    % calc hessian using secant method
    if iter > 1
%     if false
        d_obj_grad = obj_grad - grad_hist(iter-1);
        obj_hess = abs(diag(d_obj_grad./(pk*alpha)));
    else
        obj_hess = hess_init;
    end
    hess_hist(:,:, iter) = obj_hess;
    
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
    thetas = thetas + symm(pk'*alpha);
end
end_time = cputime;
run_time = end_time - start_time
%% Output
thetas
figure(1);
plot(thetas_hist(:,1:end/2))
title("Layup Orientation Optimization History")
xlabel("Iteration Number"); ylabel("Theta (degrees)");
legend("Ply 1", "Ply 2")
figure(2);
plot(obj_hist)

%% Check MS
MS = CLPT(mat_props, thetas, rad_or_deg, thicknesses, mech_loading,...
    delta_T, cte_vec, fail_crit, mat_strengths_t, mat_strengths_c, SF, ...
    print_output)

%% Test
% thetas = linspace(-90, 90, 500);
% for i = 1:length(thetas)
%     MS = CLPT(mat_props, [thetas(i) thetas(i)], rad_or_deg, thicknesses, mech_loading,...
%     delta_T, cte_vec, fail_crit, mat_strengths_t, mat_strengths_c, SF, ...
%     print_output);
%     obj(i) = -min(min(MS));
%     
%     
%     thetas_dist = [thetas(i)+theta_h thetas(i)+theta_h];
%     MS_dist = CLPT(mat_props, thetas_dist, rad_or_deg, thicknesses, ...
%         mech_loading, delta_T, cte_vec, fail_crit, mat_strengths_t, ...
%         mat_strengths_c, SF, print_output);
%     d_MS = min(min(MS_dist)) - min(min(MS));
% 
%     obj_grad(i) = -d_MS/theta_h;
% %     obj_grad(i) = sum(MS_dist(:)) - sum(MS(:));
% end
% 
% figure(4)
% plot(thetas, obj)
% figure(5)
% plot(thetas, obj_grad)

%% Test2
% thetas = linspace(-90, 90, 1000);
% for i = 1:length(thetas)
%     [stresses_bot, stresses_top] = ...
%         get_local_lamina_stresses_planar_ortho(mat_props, [thetas(i) thetas(i)], ...
%         rad_or_deg, thicknesses, mech_loading, delta_T, cte_vec);
%     obj(i) = norm(stresses_bot./(mat_props(1:3)*ones(1, 2))...
%         + stresses_top./(mat_props(1:3)*ones(1, length(thetas(i)))));
%     
%     
%     thetas_dist = [thetas(i)+theta_h thetas(i)+theta_h];
%     [stresses_bot, stresses_top] = ...
%         get_local_lamina_stresses_planar_ortho(mat_props, thetas_dist, ...
%         rad_or_deg, thicknesses, mech_loading, delta_T, cte_vec);
%     obj_dist = norm(stresses_bot./(mat_props(1:3)*ones(1, 2))...
%         + stresses_top./(mat_props(1:3)*ones(1, length(thetas(i)))));
%     d_obj = obj_dist - obj(i);
% 
%     obj_grad(i) = d_obj/theta_h;
% %     obj_grad(i) = sum(MS_dist(:)) - sum(MS(:));
% end
% 
% figure(4)
% plot(thetas, obj)
% figure(5)
% plot(thetas, obj_grad)

%% Example CLPT Analysis
% [stresses_bot, stresses_top, z_all, ...
% mid_strains_and_curvatures, thermal_loading, ABD] = ...
% get_local_lamina_stresses_planar_ortho(mat_props, thetas, ...
% rad_or_deg, thicknesses, mech_loading, delta_T, cte_vec);
% 
% [MS, failed_plies, failed_side, failed_z, fail_mode, fail_tcs] ...
% = report_ply_margins(stresses_bot, stresses_top, z_all, ...
% fail_crit, mat_strengths_t, mat_strengths_c, SF, print_output);