function [stresses_bot, stresses_top,z_all,mid_strains_and_curvatures,...
         thermal_loading, ABD,MS, failed_plies, failed_side, failed_z,...
         fail_mode, fail_tcs] = structural_model(mat_props, thetas, ...
         rad_or_deg, thicknesses, mech_loading, delta_T, cte_vec,...
         mat_strengths_t,fail_crit,mat_strengths_c, SF, print_output)
     
      [stresses_bot, stresses_top, z_all, ...
       mid_strains_and_curvatures, thermal_loading, ABD] = ...
       get_local_lamina_stresses_planar_ortho(mat_props, thetas, ...
       rad_or_deg, thicknesses, mech_loading, delta_T, cte_vec);

       [MS, failed_plies, failed_side, failed_z, fail_mode, fail_tcs] ...
       = report_ply_margins(stresses_bot, stresses_top, z_all, ...
       fail_crit, mat_strengths_t, mat_strengths_c, SF, print_output)
end