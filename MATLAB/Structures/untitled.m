% we need one ply, choose arbitrary direction, afterwards checks margin
function [MS, thetas] = structural_model_test(mat_props, thetas, ...
         rad_or_deg, thicknesses, mech_loading, delta_T, cte_vec,...
         mat_strengths_t,fail_crit,mat_strengths_c, SF, print_output);

        [stresses_bot, stresses_top,z_all,mid_strains_and_curvatures,...
         thermal_loading, ABD,MS, failed_plies, failed_side, failed_z,...
         fail_mode, fail_tcs] = structural_model(mat_props, thetas, ...
         rad_or_deg, thicknesses, mech_loading, delta_T, cte_vec,...
         mat_strengths_t,fail_crit,mat_strengths_c, SF, print_output)

wing_num_plies

% while MS < 0
% 
% num
% 
% if MS > 0
%     break
% end
% 
% end

end