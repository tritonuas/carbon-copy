function w_htail = weight_htail_analytic(num_spar, spar_width, num_plies, density_carbon_epoxy, s_h, density_blue_foam, t_htail_root, t_htail_tip, htail_span, fudgefactor)
    
    carbon_epoxy_area_density = density_carbon_epoxy*num_plies;
    
    spar_length_density = (t_htail_tip+t_htail_root/2)*num_spar*density_blue_foam*spar_width; 
    
    w_htail = 2*s_h*carbon_epoxy_area_density + htail_span*spar_length_density;
    
    w_htail = 9.81*w_htail; % convert from kg to N
    w_htail = w_htail*fudgefactor; 

end 