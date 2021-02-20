%Write comments later

function global_cte_mat = local_to_global_cte(local_cte_mat, thetas, rad_or_deg)

for k = 1:length(thetas)
    T_sigma = get_T_planar(thetas(k), 'stress', rad_or_deg);
    global_cte_mat(:,k) = T_sigma'*local_cte_mat(:,k);
end