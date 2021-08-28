%%Authors: Allyson Chen & Kevin Vo
%%Emails:  aechen@ucsd.edu    k1vo@ucsd.edu
%%Calculates the weight of the the vertical tail

%%Variable Legend
%%num_spar = number of spars
%%spar_width = width of the spar
%%num_plies = number of carbon epoxy layers
%%density_carbon_epoxy = density of the carbon epoxy composite
%%s_v = vertical tail area
%%density_blue_foam = density of the blue foam
%%t_vtail_tip = thickness of the vertical tail tip
%%t_vtail_root = thickness of the vertical tail root
%%vtail_span = vertical tail span
%%fudgefactor = fudgefactor

function w_vtail = weight_vtail_analytic(num_spar, spar_width, num_plies, density_carbon_epoxy, s_v, density_blue_foam, t_vtail_root, t_vtail_tip, vtail_span, fudgefactor)
    
    t_carbon_epoxy = 0.0003;
    carbon_epoxy_area_density = density_carbon_epoxy*num_plies*t_carbon_epoxy;
    
    spar_length_density = ((t_vtail_tip+t_vtail_root)/2)*num_spar*density_blue_foam*spar_width; 
    
    w_vtail = 2*s_v*carbon_epoxy_area_density + vtail_span*spar_length_density; % multiply area density by 2 to account for the top and bottom tail skins
    
    w_vtail = 9.81*w_vtail; % convert from kg to N
    w_vtail = w_vtail*fudgefactor; 

end 