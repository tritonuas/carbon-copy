function w_wing = weight_wing_analytic(S, num_spar,  spar_width, density_balsa, t_divinycell, density_divinycell, num_plies, density_carbon_epoxy, t_tip, t_root, wingspan, fudge_factor)
    
    carbon_epoxy_area_density = density_carbon_epoxy*num_plies; % kg/m^2 
    divinycell_area_density = t_divinycell*density_divinycell; % kg/m^2
    spar_length_density = (t_tip+t_root/2)*num_spar*density_balsa*spar_width; % kg/m
    
    
    w_wing = 2*S*(carbon_epoxy_area_density+divinycell_area_density) + wingspan*spar_length_density; % multiply area density by 2 to account for the top and bottom wing skins
    w_wing = 9.81*w_wing; % convert from kg to N
    w_wing = w_wing*fudge_factor;
end 
