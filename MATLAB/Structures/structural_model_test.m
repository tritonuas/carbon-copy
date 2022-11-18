% we need one ply, choose arbitrary direction, afterwards checks margin
function [MS, thetas] = structural_model_test(mat_props, thetas, ...
         rad_or_deg, thicknesses, mech_loading, delta_T, cte_vec,...
         mat_strengths_t,fail_crit,mat_strengths_c, SF, print_output);

        [stresses_bot, stresses_top,z_all,mid_strains_and_curvatures,...
         thermal_loading, ABD,MS, failed_plies, failed_side, failed_z,...
         fail_mode, fail_tcs] = structural_model(mat_props, thetas, ...
         rad_or_deg, thicknesses, mech_loading, delta_T, cte_vec,...
         mat_strengths_t,fail_crit,mat_strengths_c, SF, print_output)

%
wing_num_plies = 0;
thetas = [0 90];
MS = 0;

while MS < 0

    thetas = [0 90]

    [stresses_bot, stresses_top,z_all,mid_strains_and_curvatures,...
    thermal_loading, ABD, MS, failed_plies, failed_side, failed_z,...
    fail_mode, fail_tcs] = structural_model(mat_props, thetas, ...
    rad_or_deg, thicknesses, mech_loading, delta_T, cte_vec,...
    mat_strengths_t,fail_crit,mat_strengths_c, SF, print_output)
    
    MS_0_90 = MS;
    thetas = [45 45];

    [stresses_bot, stresses_top,z_all,mid_strains_and_curvatures,...
    thermal_loading, ABD, MS, failed_plies, failed_side, failed_z,...
    fail_mode, fail_tcs] = structural_model(mat_props, thetas, ...
    rad_or_deg, thicknesses, mech_loading, delta_T, cte_vec,...
    mat_strengths_t,fail_crit,mat_strengths_c, SF, print_output)

    MS_45_45 = MS;

%     is MS_1 > MS_2

    %determine which MS is higher and which one is greater than 0
    % if it is greater than 0 it wont fail
    % if it fails add another ply
    % if 

    % for num_plies = 1:i

% 
% if MS > 0
%     break
% end
% 
% end
    wing_num_plies = wing_num_plies + 1;
end