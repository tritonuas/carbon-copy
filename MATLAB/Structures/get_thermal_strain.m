%Write comments later

function thermal_strains = get_thermal_strain(delta_T, cte_mat_local, ...
    thetas, rad_or_deg)

cte_mat_global = local_to_global_cte(cte_mat_local, thetas, rad_or_deg);
thermal_strains = cte_mat_global*delta_T;