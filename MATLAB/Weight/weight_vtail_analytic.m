function w_vtail = weight_vtail_analytic(num_spar, spar_width, num_plies, density_carbon_epoxy, s_v, density_blue_foam, t_vtail_root, t_vtail_tip, vtail_span, fudgefactor)
    
    carbon_epoxy_area_density = density_carbon_epoxy*num_plies;
    
    spar_length_density = (t_vtail_tip+t_vtail_root/2)*num_spar*density_blue_foam*spar_width; 
    
    w_vtail = 2*s_v*carbon_epoxy_area_density + vtail_span*spar_length_density;
    
    w_vtail = 9.81*w_vtail; % convert from kg to N
    w_vtail = w_vtail*fudgefactor; 

end 