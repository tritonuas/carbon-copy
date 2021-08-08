
function gross = compute_weight_analytic(battery, payload, S, num_spar_wing,  spar_width_wing, density_balsa, t_divinycell, density_divinycell,...
    num_plies_wing, density_carbon_epoxy, t_tip, t_root, wingspan, fudge_factor,...
    num_spar_htail, spar_width_htail, num_plies_htail, s_h, density_blue_foam, t_htail_root, t_htail_tip, htail_span,...
    num_spar_vtail, spar_width_vtail, num_plies_vtail, s_v, t_vtail_root, t_vtail_tip, vtail_span,...
    num_plies_fuse, num_bulkheads, t_bulkhead, area_fraction_bulkhead, density_plywood, fuse_height, fuse_width, fuse_length,...
    wing_fudge_factor, htail_fudge_factor, vtail_fudge_factor, fuse_fudge_factor)
    

    w_wing = weight_wing_v2(S, num_spar_wing,  spar_width_wing, density_balsa, t_divinycell, density_divinycell, num_plies_wing, density_carbon_epoxy, t_tip, t_root, wingspan, wing_fudge_factor);
    w_htail = weight_htail_v2(num_spar_htail, spar_width_htail, num_plies_htail, density_carbon_epoxy, s_h, density_blue_foam, t_htail_root, t_htail_tip, htail_span, htail_fudge_factor);
    w_vtail = weight_vtail_v2(num_spar_vtail, spar_width_vtail, num_plies_vtail, density_carbon_epoxy, s_v, density_blue_foam, t_vtail_root, t_vtail_tip, vtail_span, vtail_fudge_factor);
    w_fuse = weight_fuse_v2(density_carbon_epoxy, num_plies_fuse, num_bulkheads, t_bulkhead, area_fraction_bulkhead, density_plywood, fuse_height, fuse_width, fuse_length, fuse_fudge_factor);
    
    gross = w_wing + w_htail + w_vtail + w_fuse + battery + payload;
    gross = gross*fudge_factor;
end