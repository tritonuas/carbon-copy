%%Authors: Allyson Chen & Kevin Vo
%%Emails:  aechen@ucsd.edu    k1vo@ucsd.edu
%%Calculates the weight of the the horizontal tail

%%Variable Legend
%%w_htail = weight of horizontal tail
%%num_spar = number of spars
%%spar_width = width of the spar
%%num_plies = number of carbon epoxy layers
%%density_carbon_epoxy = density of the carbon epoxy composite
%%s_h = horizontal tail area
%%density_blue_foam = density of the blue foam
%%t_htail_tip = thickness of the horizontal tail tip
%%t_htail_root = thickness of the horizontal tail root
%%htail_span = horizontal tail span
%%fudgefactor = fudgefactor


function w_htail = weight_htail_analytic(num_spar, spar_width, num_plies, density_carbon_epoxy, s_h, density_blue_foam, t_htail_root, t_htail_tip, htail_span, fudgefactor)
    
    t_carbon_epoxy = 0.0003;
    carbon_epoxy_area_density = density_carbon_epoxy*num_plies*t_carbon_epoxy;
    
    spar_length_density = ((t_htail_tip+t_htail_root)/2)*num_spar*density_blue_foam*spar_width; 
    
    w_htail = 2*s_h*carbon_epoxy_area_density + htail_span*spar_length_density; % multiply area density by 2 to account for the two sides of tail skins
    skin = 2*s_h*carbon_epoxy_area_density;
    htail = htail_span*spar_length_density;
    
    w_htail = 9.81*w_htail; % convert from kg to N
    w_htail = w_htail*fudgefactor; 

end 